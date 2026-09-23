#!/usr/bin/env python3
"""
download_ncbi_fasta.py

Batch-downloads FASTA sequences from NCBI (via E-utilities' efetch) for a
list of accessions, with automatic retries for anything that fails or is
missing, so a single run finishes with everything found -- or clearly
reports what's still missing -- instead of needing a separate manual
retry pass.

Why this exists: a naive one-shot chunked download over thousands of
accessions predictably leaves some records missing (dropped requests,
NCBI rate-limit responses, truncated chunks, ...). Rather than writing
one script to download, another to diff against the original list, and
a third to retry the gaps, this script folds all of that into a single
run: after the first pass, whatever is still missing is retried in
progressively smaller chunks (down to one accession per request on the
final passes) until either everything is found or --max-retries is hit.

Requirements:
    pip install requests tqdm      # tqdm is optional; a progress bar is
                                    # shown if it's installed, skipped if not

Usage:
    python3 download_ncbi_fasta.py -i accessions.txt -o sequences.fasta
    python3 download_ncbi_fasta.py -i accessions.txt -o sequences.fasta \\
        --workers 4 --chunk-size 200 --max-retries 5

Input file format: one NCBI accession per line (e.g. NC_012920.1).
Output: a single FASTA file, in the same order as the input list.
If any accessions could not be retrieved after --max-retries passes,
they are written to "<output>.missing.txt" instead of silently dropped.
"""

import argparse
import sys
import threading
import time
from pathlib import Path

import requests

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI's commonly-cited rate limits for E-utilities requests per second:
# 3/sec without an API key, 10/sec with one. See NCBI's E-utilities
# documentation for the current, authoritative limits.
RATE_LIMIT_NO_KEY = 3.0
RATE_LIMIT_WITH_KEY = 10.0


class RateLimiter:
    """Thread-safe rate limiter: blocks the calling thread until at
    least `1 / max_per_second` seconds have passed since the last call,
    across ALL threads sharing this instance. This -- not the thread
    count -- is what keeps the total request rate within NCBI's limit
    when multiple worker threads are downloading concurrently."""

    def __init__(self, max_per_second: float):
        self._min_interval = 1.0 / max_per_second
        self._lock = threading.Lock()
        self._last_call = 0.0

    def wait(self):
        with self._lock:
            now = time.monotonic()
            sleep_for = self._min_interval - (now - self._last_call)
            if sleep_for > 0:
                time.sleep(sleep_for)
            self._last_call = time.monotonic()


def parse_fasta_records(fasta_text: str) -> dict:
    """Splits a multi-record FASTA string into {accession: record_text}.
    The accession is the first whitespace-delimited token after '>' on
    each header line -- NCBI echoes back exactly the accession you asked
    for there, so this can be compared directly against the request list."""
    records = {}
    current_id = None
    current_lines = []

    for line in fasta_text.splitlines(keepends=True):
        if line.startswith(">"):
            if current_id is not None:
                records[current_id] = "".join(current_lines)
            current_id = line[1:].split()[0].strip()
            current_lines = [line]
        else:
            current_lines.append(line)

    if current_id is not None:
        records[current_id] = "".join(current_lines)

    return records


def fetch_chunk(session: requests.Session, accessions: list, api_key: str,
                 rate_limiter: RateLimiter, timeout: int = 60) -> dict:
    """Fetches one chunk of accessions from NCBI efetch in a single
    request (POST, so large ID lists don't hit URL-length limits), and
    returns {accession: fasta_record_text} for whichever of them were
    actually present in the response."""
    data = {
        "db": "nucleotide",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
    }
    if api_key:
        data["api_key"] = api_key

    rate_limiter.wait()
    response = session.post(EFETCH_URL, data=data, timeout=timeout)
    response.raise_for_status()

    text = response.text
    if ">" not in text:
        raise ValueError("response did not contain any FASTA records")

    return parse_fasta_records(text)


def download_with_retries(accessions: list, session: requests.Session,
                           api_key: str, rate_limiter: RateLimiter,
                           chunk_size: int, max_retries: int,
                           progress=None):
    """
    Attempts to download every accession in `accessions`, retrying
    whatever is still missing after each pass with a smaller chunk size
    (shrinking toward one accession per request on the final passes), up
    to `max_retries` passes.

    Returns (results, still_missing):
      results: {accession: fasta_record_text} for everything found
      still_missing: list of accessions never successfully retrieved
    """
    results = {}
    remaining = list(accessions)
    current_chunk_size = max(1, chunk_size)

    for attempt in range(1, max_retries + 1):
        if not remaining:
            break

        chunks = [remaining[i:i + current_chunk_size]
                  for i in range(0, len(remaining), current_chunk_size)]

        newly_missing = []
        for chunk in chunks:
            try:
                found = fetch_chunk(session, chunk, api_key, rate_limiter)
            except Exception as exc:
                print(f"  [!] Pass {attempt}: chunk of {len(chunk)} failed "
                      f"({exc}); will retry individually next pass.")
                newly_missing.extend(chunk)
                time.sleep(1.0)  # brief extra backoff after any failure
                continue

            for acc in chunk:
                if acc in found:
                    results[acc] = found[acc]
                    if progress is not None:
                        progress.update(1)
                else:
                    newly_missing.append(acc)

        remaining = newly_missing

        # Shrink the chunk size each pass: a smaller chunk means a
        # dropped request or a rate-limit response affects fewer
        # accessions, and the last couple of passes fall back to one
        # accession per request so a single stubborn ID can't keep
        # blocking the rest of the batch.
        current_chunk_size = max(1, current_chunk_size // 4)

    return results, remaining


def read_accession_list(path: Path) -> list:
    with open(path) as f:
        # dict.fromkeys() dedupes while preserving the original order
        return list(dict.fromkeys(line.strip() for line in f if line.strip()))


def main():
    parser = argparse.ArgumentParser(
        description="Batch-download FASTA sequences from NCBI for a list "
                    "of accessions, automatically retrying whatever is "
                    "missing after the first pass."
    )
    parser.add_argument("-i", "--input", required=True,
                         help="Text file with one NCBI accession per line")
    parser.add_argument("-o", "--output", required=True,
                         help="Output FASTA file")
    parser.add_argument("--chunk-size", type=int, default=200,
                         help="Accessions per request on the first pass "
                              "(default: 200)")
    parser.add_argument("--workers", type=int, default=3,
                         help="Number of concurrent download threads "
                              "(default: 3). Requests are still throttled "
                              "to NCBI's overall rate limit regardless of "
                              "this value -- it controls concurrency, not "
                              "the total request rate.")
    parser.add_argument("--max-retries", type=int, default=5,
                         help="Maximum number of passes over any "
                              "still-missing accessions (default: 5)")
    parser.add_argument("--api-key", default="",
                         help="Optional NCBI API key, raises the rate "
                              "limit from 3 to 10 requests/second. Pass "
                              "your own key on the command line or via an "
                              "environment variable -- never hardcode a "
                              "real key in a script committed to a public "
                              "repo.")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"Error: input file not found: {input_path}")

    accessions = read_accession_list(input_path)
    if not accessions:
        sys.exit(f"Error: no accessions found in {input_path}")

    print(f"[*] {len(accessions)} unique accessions to download")

    rate_limit = RATE_LIMIT_WITH_KEY if args.api_key else RATE_LIMIT_NO_KEY
    rate_limiter = RateLimiter(rate_limit)

    session = requests.Session()
    session.headers.update({"User-Agent": "download_ncbi_fasta.py"})

    progress = (tqdm(total=len(accessions), unit="acc", desc="Downloading")
                if tqdm else None)

    # Split the accession list across worker threads up front; each
    # worker runs its own share through download_with_retries (including
    # its own retry passes). The shared RateLimiter -- not the thread
    # count -- is what keeps the combined request rate within NCBI's
    # limit.
    worker_count = max(1, args.workers)
    shares = [accessions[i::worker_count] for i in range(worker_count)]
    shares = [s for s in shares if s]

    all_results = {}
    all_missing = []
    results_lock = threading.Lock()

    def run_share(share):
        results, missing = download_with_retries(
            share, session, args.api_key, rate_limiter,
            args.chunk_size, args.max_retries, progress,
        )
        with results_lock:
            all_results.update(results)
            all_missing.extend(missing)

    threads = [threading.Thread(target=run_share, args=(share,))
               for share in shares]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    if progress is not None:
        progress.close()

    # Write output in the same order as the input list, skipping
    # whatever ended up permanently missing.
    output_path = Path(args.output)
    with open(output_path, "w") as out:
        for acc in accessions:
            if acc in all_results:
                out.write(all_results[acc])

    print(f"\n[Done] Wrote {len(all_results)} / {len(accessions)} "
          f"sequences to {output_path}")

    if all_missing:
        unique_missing = sorted(set(all_missing))
        missing_path = output_path.with_name(output_path.name + ".missing.txt")
        with open(missing_path, "w") as f:
            for acc in unique_missing:
                f.write(f"{acc}\n")
        print(f"[!] {len(unique_missing)} accessions could not be "
              f"retrieved after {args.max_retries} passes.")
        print(f"    See: {missing_path}")
        print(f"    You can re-run with just that file as --input to "
              f"try again (e.g. after NCBI traffic has settled down).")
    else:
        print("[OK] All accessions retrieved successfully.")


if __name__ == "__main__":
    main()
