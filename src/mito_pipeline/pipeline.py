"""
Batch pipeline for mitochondrial genome reconstruction from Nanopore reads:

  1. find_best_reference() -- maps a sample's reads against a reference
     database and selects the single best-matching reference sequence.
  2. run_medaka_iterative() -- iteratively polishes that reference against
     the sample's reads with Medaka, feeding a quality-filtered read set
     into each successive round, then depth-trims unreliable 5'/3' ends.
  3. main() -- runs both steps over every FASTQ file in a directory.

Output naming:
    {prefix}_ref.fasta               reference selected for this sample
    {prefix}.consensus_untrimmed.fasta   consensus before depth trimming
    {prefix}.final.fasta             the pipeline's actual result (trimmed)

Only "{prefix}.final.fasta" is meant to be used downstream (e.g. submitted
alongside the raw reads); the untrimmed consensus is kept alongside it for
reference/comparison only.

Usage (after `pip install runmitopipe`):
    runmitopipe \
        -i raw_fastq_dir \
        -r reference_db.fasta \
        -o output_base \
        --threads 64 --min_count 100 --max_iter 5
"""

import argparse
import subprocess
import shutil
import sys
import traceback
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


DEFAULT_MEDAKA_MODEL = "r1041_e82_400bps_sup_v5.2.0"  # medaka_consensus's own default

# External command-line programs this script (pipeline.py) actually shells
# out to. These are compiled binaries -- pip cannot install them, users need
# to install them themselves via conda or a system package manager.
# Note: if other scripts are added to this project later (e.g. ones that use
# blast or seqkit), each of those scripts should declare its own list of
# required tools and pass it into check_external_tools() -- there's no need
# to merge everything into one project-wide list here.
REQUIRED_EXTERNAL_TOOLS = ["minimap2", "samtools", "medaka_consensus"]

# Fallback thresholds for samples where mitochondrial reads are a small
# fraction of the input (e.g. reads pulled from whole-genome sequencing
# data rather than from mtDNA-enriched sequencing), where the default
# thresholds below can leave zero reads passing per-iteration filtering.
RELAXED_FILTER_PARAMS = dict(
    min_mapq=10,
    min_aligned_fraction=0.5,
    max_softclip_fraction=0.4,
    min_aligned_bases=100,
    min_identity=0.6,
)


class PipelineError(Exception):
    """Raised for an anticipated, sample-level failure (e.g. no reads pass
    filtering). Caught in main() so one failed sample does not abort the
    rest of the batch."""


def check_external_tools(tools: list[str] = None) -> list[str]:
    """Checks whether the given external tools are all on PATH, and returns
    the list of any that are missing.

    `tools` defaults to this script's own REQUIRED_EXTERNAL_TOOLS, but a
    different list can be passed in to reuse this function elsewhere --
    e.g. if other scripts using blast or seqkit are added to this project
    later, each can pass in its own required-tools list rather than every
    script's requirements being merged into one combined list (so a user
    who only installed minimap2/samtools/medaka, and not blast, won't be
    told to install blast when running a script that has nothing to do
    with it).
    """
    if tools is None:
        tools = REQUIRED_EXTERNAL_TOOLS
    return [tool for tool in tools if shutil.which(tool) is None]


class _Tee:
    """Writes to multiple streams at once (used to mirror printed output
    into a per-sample log file, in addition to the console)."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()


@contextmanager
def log_to_file(log_path):
    """Mirrors everything printed to stdout/stderr within this block into
    `log_path`, in addition to the console, so a sample's full processing
    log is preserved even though the FASTA outputs themselves stay clean."""
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_file:
        original_stdout, original_stderr = sys.stdout, sys.stderr
        sys.stdout = _Tee(original_stdout, log_file)
        sys.stderr = _Tee(original_stderr, log_file)
        try:
            yield
        finally:
            sys.stdout, sys.stderr = original_stdout, original_stderr


# ============================================================
# Step 1: select the best-matching reference for a sample
# ============================================================

def run_minimap2_to_paf_df(ref_fasta, query_fastq, threads=64):
    """Runs minimap2 and parses PAF output into a DataFrame (12 standard
    columns, any extra tags collected into a single 'Tags' column)."""
    core_columns = [
        'Qname', 'Qlen', 'Qstart', 'Qend', 'strand',
        'Tname', 'Tlen', 'Tstart', 'Tend', 'Nmatch',
        'Alen', 'MapQ'
    ]

    print(f"Running minimap2 with {threads} threads...")
    cmd = ["minimap2", "-t", str(threads), "-x", "map-ont", str(ref_fasta), str(query_fastq)]

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    paf_data = []
    for line in process.stdout:
        parts = line.strip().split('\t')
        if len(parts) < 12:
            continue
        core_data = parts[:12]
        tags = parts[12:] if len(parts) > 12 else []
        paf_data.append(core_data + [tags])

    process.stdout.close()
    process.wait()
    return pd.DataFrame(paf_data, columns=core_columns + ['Tags'])


def find_best_reference(fastq, ref_db, output_ref_fasta, threads=64, min_count=100):
    """
    Maps `fastq` against `ref_db`, finds the reference sequence with the most
    reads whose best hit lands on it, and writes that single winning
    reference sequence to `output_ref_fasta` if it clears `min_count` reads.

    Returns True if a winning reference was found and written, False
    otherwise (no output file is written in the False case).
    """
    df = run_minimap2_to_paf_df(ref_db, fastq, threads)

    if df.empty:
        print("No alignments found.")
        return False

    numeric_cols = ['Qlen', 'Qstart', 'Qend', 'Tlen', 'Tstart', 'Tend', 'Nmatch', 'Alen']
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Keep only reasonably long alignments, then keep each read's single
    # best hit so every read votes for exactly one reference sequence.
    df = df[df['Alen'] > 500]
    df = df.sort_values(by=['Qname', 'Alen'], ascending=[True, False])
    df_unique = df.drop_duplicates(subset=['Qname'])

    tname_counts = df_unique['Tname'].value_counts()
    if tname_counts.empty:
        print("No Tname passed the filter criteria.")
        return False

    best_tname = tname_counts.index[0]
    best_count = tname_counts.iloc[0]

    print("\n" + "=" * 40)
    print(f"{'Rank':<5} {'Reference Name':<30} {'Read Count':<10}")
    for i, (name, count) in enumerate(tname_counts.head(3).items(), 1):
        print(f"{i:<5} {name:<30} {count:<10}")
    print("=" * 40)

    if best_count < min_count:
        print(f"\n[!] Threshold not met: best reference '{best_tname}' only has {best_count} reads "
              f"(required: {min_count}). No output file generated.")
        return False

    print(f"\nWinning reference: {best_tname} (count: {best_count})")

    ref_dict = SeqIO.to_dict(SeqIO.parse(str(ref_db), "fasta"))
    if best_tname not in ref_dict:
        print(f"Error: {best_tname} not found in reference database.")
        return False

    Path(output_ref_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(ref_dict[best_tname], str(output_ref_fasta), "fasta")
    print(f"Result written to: {output_ref_fasta}")
    return True


# ============================================================
# Step 2: iterative Medaka polishing with per-read quality
# filtering feeding into the next iteration
# ============================================================

def read_fasta_as_set(fasta_path):
    """Reads a FASTA's sequence content as a set, used to detect when two
    iterations produced identical consensus sequences."""
    return {str(record.seq).upper() for record in SeqIO.parse(str(fasta_path), "fasta")}


def run_medaka(fastq, reference, outdir, threads=64, quiet=True, model=DEFAULT_MEDAKA_MODEL):
    """
    Runs medaka_consensus.

    -m must match the flow cell/basecaller chemistry used to generate
    `fastq` (e.g. an R9.4.1 model such as 'r941_min_hac_g507' for R9 data);
    medaka does not warn on a mismatched model, it silently produces a
    lower-quality result.

    -f forces medaka to start from a clean state rather than reusing any
    stale output already present in `outdir`, which is medaka's own default
    behavior when the output directory already exists.
    """
    outdir = Path(outdir)
    cmd = [
        "medaka_consensus",
        "-i", str(fastq),
        "-d", str(reference),
        "-r", "N",
        "-m", model,
        "-f",
        "-t", str(threads),
        "-o", str(outdir),
    ]

    print(f"[*] Running Medaka on {Path(reference).name} (model: {model}) ...")

    if quiet:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print("[X] medaka_consensus failed:")
            print("----- stdout -----")
            print(result.stdout)
            print("----- stderr -----")
            print(result.stderr)
            raise subprocess.CalledProcessError(result.returncode, cmd)
    else:
        subprocess.run(cmd, check=True)


def get_softclip_lengths(cigartuples):
    """Total soft-clipped length from a CIGAR (op 4 = soft-clip)."""
    if cigartuples is None:
        return 0
    return sum(length for op, length in cigartuples if op == 4)


def get_max_indel_length(cigartuples):
    """Largest single indel length in a CIGAR (op 1 = insertion, op 2 = deletion)."""
    if cigartuples is None:
        return 0
    return max([length for op, length in cigartuples if op in (1, 2)] + [0])


def passes_alignment_filters(aln, min_mapq=20, min_aligned_fraction=0.80,
                              max_softclip_fraction=0.15, min_aligned_bases=200,
                              min_identity=0.85, max_indel_size=None):
    """
    Judges whether a single (primary) alignment record is good enough to
    feed into the next polishing iteration. Identity is estimated from the
    NM tag as 1 - NM/aligned_length.
    """
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        return False
    if aln.mapping_quality < min_mapq:
        return False

    query_len = aln.query_length
    if not query_len:
        return False

    aligned_bases = aln.query_alignment_length
    if not aligned_bases or aligned_bases < min_aligned_bases:
        return False

    if (aligned_bases / query_len) < min_aligned_fraction:
        return False

    if (get_softclip_lengths(aln.cigartuples) / query_len) > max_softclip_fraction:
        return False

    if max_indel_size and get_max_indel_length(aln.cigartuples) > max_indel_size:
        return False

    if aln.has_tag("NM"):
        identity = 1.0 - (aln.get_tag("NM") / aligned_bases)
        if identity < min_identity:
            return False

    return True


def filter_bam_to_fastq(bam_path, output_fastq, **kwargs):
    """
    Writes reads from `bam_path` that pass passes_alignment_filters() to
    `output_fastq` (one FASTQ record per read, deduplicated by name).

    Returns the number of reads kept. Does not raise on zero reads kept --
    callers decide how to respond (e.g. retry with relaxed thresholds).
    """
    import pysam

    seen_reads = set()
    forward_count = 0
    reverse_count = 0
    kept = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam, open(output_fastq, "w") as out:
        for aln in bam.fetch(until_eof=True):
            if aln.query_name in seen_reads:
                continue
            if not passes_alignment_filters(aln, **kwargs):
                continue

            seq = aln.query_sequence
            qual_array = aln.query_qualities
            if seq is None or qual_array is None:
                continue

            qual = pysam.qualities_to_qualitystring(qual_array)
            out.write(f"@{aln.query_name}\n{seq}\n+\n{qual}\n")

            seen_reads.add(aln.query_name)
            kept += 1
            if aln.is_reverse:
                reverse_count += 1
            else:
                forward_count += 1

    print(f"[OK] Kept {kept} reads (Fwd: {forward_count}, Rev: {reverse_count}) -> {output_fastq}")
    return kept


def compute_depth_profile(bam_path, contig_name):
    """Per-position read depth for `contig_name`, via `samtools depth -a`
    (the -a flag reports every position, including zero-depth ones)."""
    result = subprocess.run(
        ["samtools", "depth", "-a", "-r", contig_name, str(bam_path)],
        capture_output=True, text=True, check=True,
    )
    depths = []
    for line in result.stdout.strip().splitlines():
        if line:
            depths.append(int(line.split("\t")[2]))
    return depths


def trim_low_depth_ends(consensus_fasta, bam_path, output_fasta, min_depth=5):
    """
    Trims bases from the 5' and 3' ends of a single-contig consensus wherever
    per-base read depth falls below min_depth, scanning inward from each end
    until the threshold is met. Low-depth ends are the least reliable part
    of a consensus (e.g. reads not fully spanning the linearization point of
    a circular genome) and are discarded rather than reported.

    Returns (n_trimmed_5prime, n_trimmed_3prime, final_length).
    """
    records = list(SeqIO.parse(str(consensus_fasta), "fasta"))
    if len(records) != 1:
        print(f"[!] Warning: expected exactly 1 contig in {consensus_fasta}, found {len(records)}. "
              f"Trimming only the first; others are dropped from the output.")
    record = records[0]
    seq = str(record.seq)
    seq_len = len(seq)

    depths = compute_depth_profile(bam_path, record.id)
    if len(depths) != seq_len:
        print(f"[!] Warning: depth profile length ({len(depths)}) does not match consensus "
              f"length ({seq_len}) for '{record.id}' -- skipping trimming, copying unchanged.")
        SeqIO.write(record, str(output_fasta), "fasta")
        return 0, 0, seq_len

    start_trim = 0
    while start_trim < seq_len and depths[start_trim] < min_depth:
        start_trim += 1

    end_trim = 0
    while end_trim < seq_len - start_trim and depths[seq_len - 1 - end_trim] < min_depth:
        end_trim += 1

    trimmed_seq = seq[start_trim: seq_len - end_trim]
    final_len = len(trimmed_seq)

    print(f"[*] Depth-based end trimming (min_depth={min_depth}x):")
    print(f"    5' end trimmed: {start_trim} bp | 3' end trimmed: {end_trim} bp")
    print(f"    Length: {seq_len} bp -> {final_len} bp")

    if final_len == 0:
        print(f"    [!] WARNING: entire sequence fell below {min_depth}x depth -- output is EMPTY.")
    elif final_len < 0.5 * seq_len:
        print(f"    [!] WARNING: trimmed sequence lost more than half its length -- "
              f"check this sample's coverage.")

    # Header/description are carried through unchanged -- trimming details
    # are printed above (and captured in the run log) rather than embedded
    # in the FASTA itself.
    trimmed_record = SeqRecord(Seq(trimmed_seq), id=record.id, description=record.description)
    SeqIO.write(trimmed_record, str(output_fasta), "fasta")
    return start_trim, end_trim, final_len


def align_reads_to_reference(fastq, reference, output_bam, threads=8):
    """
    Aligns `fastq` to `reference` with minimap2, producing a sorted, indexed
    BAM guaranteed to share `reference`'s own coordinate system.

    Used specifically for depth-based trimming: Medaka's own
    "calls_to_draft.bam" is aligned to the DRAFT reference fed INTO an
    iteration, not to the CONSENSUS it produced, and those two can differ by
    a base or more whenever Medaka's polishing edits the sequence -- which
    breaks the length match that depth-based trimming depends on.
    """
    reference = Path(reference)
    subprocess.run(["samtools", "faidx", str(reference)],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cmd = (
        f"minimap2 -ax map-ont -t {threads} {reference} {fastq} 2>/dev/null "
        f"| samtools sort -@ {threads} -o {output_bam} - 2>/dev/null"
    )
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(["samtools", "index", str(output_bam)],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _handle_zero_reads(filter_kwargs, interactive):
    """
    Called when filter_bam_to_fastq() keeps zero reads. This is most often
    caused by thresholds tuned for mtDNA-enriched sequencing being too
    strict for a sample where mitochondrial reads are a small fraction of
    the input (e.g. reads recovered from whole-genome sequencing data).

    Prints the current thresholds alongside a relaxed set known to work for
    such low-signal samples, and either prompts the user to retry with them
    (interactive=True) or raises PipelineError so the caller can skip this
    sample and continue the batch.

    Returns the filter_kwargs to retry with, or raises PipelineError.
    """
    print("\n[!] No reads passed the per-iteration alignment filters.")
    print("    Current thresholds:")
    for k, v in filter_kwargs.items():
        print(f"      {k} = {v}")
    print("    This commonly happens when mitochondrial reads are a small fraction of the "
          "input (e.g. whole-genome sequencing data). A relaxed threshold set that has "
          "worked for such samples:")
    for k, v in RELAXED_FILTER_PARAMS.items():
        print(f"      {k} = {v}")

    if interactive:
        answer = input("\n    Retry this sample with the relaxed thresholds above? [y/N]: ").strip().lower()
        if answer == "y":
            return {**filter_kwargs, **RELAXED_FILTER_PARAMS}

    raise PipelineError(
        "No reads passed alignment filtering, and no relaxed retry was applied. "
        "Consider re-running this sample with relaxed thresholds, e.g.: " +
        " ".join(f"--{k} {v}" for k, v in RELAXED_FILTER_PARAMS.items())
    )


def rename_fasta_header(fasta_path, new_id):
    """Overwrites a single-contig FASTA's header with `new_id` (both the id
    and description), so the sequence is identified only by the sample it
    came from -- provenance details (selected reference, trimming) live in
    the run's log file instead of being embedded in the FASTA itself."""
    record = next(SeqIO.parse(str(fasta_path), "fasta"))
    record.id = new_id
    record.description = new_id
    SeqIO.write(record, str(fasta_path), "fasta")


def run_medaka_iterative(fastq, init_ref, out_prefix, sample_name=None,
                          max_iter=10, threads=64,
                          quiet=True, model=DEFAULT_MEDAKA_MODEL, min_mapq=20,
                          min_aligned_fraction=0.80, max_softclip_fraction=0.15,
                          min_aligned_bases=200, min_identity=0.85, max_indel_size=None,
                          min_trim_depth=5, interactive=False):
    """
    Iteratively polishes `init_ref` against `fastq` with Medaka, up to
    `max_iter` times.

    Each iteration's alignments are quality-filtered (mapping quality,
    aligned length/fraction, soft-clip fraction, identity, indel size), and
    the resulting FASTQ becomes the input for the next iteration --
    progressively excluding poorly-aligned or chimeric reads from the
    polishing pool. If a round of filtering keeps zero reads, see
    _handle_zero_reads().

    Convergence is declared once the current iteration's consensus matches
    ANY previously-seen iteration's consensus (not just the immediately
    preceding one), since Medaka's consensus can oscillate at homopolymer
    positions rather than settle monotonically.

    Writes "{out_prefix}.consensus_untrimmed.fasta" (pre-trim),
    "{out_prefix}.final.fasta" (the pipeline's actual result, after
    depth-based end trimming), and "{out_prefix}.filtered_reads.fastq" (the
    read set that produced the winning iteration). If `sample_name` is
    given, both output FASTA headers are renamed to it (the selected
    reference's own header, e.g. an NCBI accession, is kept in the
    description instead).

    Returns the Path to "{out_prefix}.final.fasta" on success, or None if
    Medaka never produced a consensus.
    """
    current_fastq = Path(fastq)
    prev_ref = Path(init_ref)
    seen_seqs = []       # [(iteration_index, sequence_set), ...] -- full convergence history
    iteration_dirs = []
    temp_fastqs = []

    filter_kwargs = dict(
        min_mapq=min_mapq, min_aligned_fraction=min_aligned_fraction,
        max_softclip_fraction=max_softclip_fraction, min_aligned_bases=min_aligned_bases,
        min_identity=min_identity, max_indel_size=max_indel_size,
    )

    for i in range(max_iter):
        outdir = Path(f"{out_prefix}-{i}")
        iteration_dirs.append(outdir)
        print(f"\n========== Iteration {i} ==========")

        run_medaka(current_fastq, prev_ref, outdir, threads=threads, quiet=quiet, model=model)

        consensus = outdir / "consensus.fasta"
        bam_file = outdir / "calls_to_draft.bam"
        if not consensus.exists():
            raise PipelineError(f"consensus.fasta missing in {outdir}")

        curr_seq = read_fasta_as_set(consensus)

        match_iter = next((j for j, s in seen_seqs if s == curr_seq), None)
        if match_iter is not None:
            gap = i - match_iter
            note = "" if gap == 1 else f", {gap} iterations apart -- likely homopolymer oscillation"
            print(f"[OK] Consensus converged at iteration {i} (matches iteration {match_iter}{note})")
            prev_ref = consensus
            break

        seen_seqs.append((i, curr_seq))
        prev_ref = consensus

        next_fastq = Path(f"{out_prefix}.iter{i + 1}.filtered.fastq")
        kept = filter_bam_to_fastq(bam_file, next_fastq, **filter_kwargs)
        if kept == 0:
            filter_kwargs = _handle_zero_reads(filter_kwargs, interactive)
            kept = filter_bam_to_fastq(bam_file, next_fastq, **filter_kwargs)
            if kept == 0:
                raise PipelineError("No reads passed even the relaxed thresholds. Convergence impossible.")

        if current_fastq in temp_fastqs and current_fastq.exists():
            print(f"[*] Removing intermediate FASTQ: {current_fastq}")
            current_fastq.unlink()

        current_fastq = next_fastq
        temp_fastqs.append(current_fastq)

    untrimmed_fasta = Path(f"{out_prefix}.consensus_untrimmed.fasta")
    shutil.copyfile(prev_ref, untrimmed_fasta)

    # Medaka's consensus otherwise keeps the selected reference's own header
    # (e.g. an NCBI accession) rather than identifying the sample it was
    # reconstructed from. Renaming here means both the untrimmed and final
    # (trimmed) outputs carry the sample's own name.
    if sample_name:
        rename_fasta_header(untrimmed_fasta, sample_name)

    print(f"\n[*] Untrimmed consensus: {untrimmed_fasta}")

    # Realign the read set that produced this consensus specifically against
    # it (rather than reusing this iteration's calls_to_draft.bam, which is
    # aligned to the DRAFT fed in, not this OUTPUT) to get a depth profile
    # guaranteed to match its coordinate system.
    depth_bam = Path(f"{out_prefix}.depth_check.bam")
    align_reads_to_reference(current_fastq, untrimmed_fasta, depth_bam, threads=threads)

    final_fasta = Path(f"{out_prefix}.final.fasta")
    trim_low_depth_ends(untrimmed_fasta, depth_bam, final_fasta, min_depth=min_trim_depth)

    # current_fastq is always one of the per-iteration temp files by this
    # point (never the original raw input), so renaming it in place to a
    # predictable name -- rather than leaving it as "...iterN.filtered.fastq"
    # with N varying per sample -- is safe.
    filtered_reads_fastq = Path(f"{out_prefix}.filtered_reads.fastq")
    if current_fastq.exists():
        current_fastq.rename(filtered_reads_fastq)
        print(f"[*] Filtered reads used for final consensus: {filtered_reads_fastq}")
        current_fastq = filtered_reads_fastq

    for tf in temp_fastqs:
        if tf.exists() and tf != current_fastq:
            tf.unlink()

    if len(iteration_dirs) > 2:
        for d in iteration_dirs[:-2]:
            if d.exists():
                print(f"[*] Cleaning up old directory: {d}")
                shutil.rmtree(d)

    return final_fasta if final_fasta.exists() else None


# ============================================================
# Step 3: batch orchestration over a directory of FASTQ samples
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        prog="runmitopipe",
        description="Mitochondrial genome reconstruction: reference selection -> "
                    "iterative Medaka polishing -> depth trimming, batched over a "
                    "directory of FASTQ samples."
    )
    parser.add_argument("-i", "--input_dir", required=True, help="Input FASTQ directory")
    parser.add_argument("-r", "--reference_db", required=True, help="Path to reference database")
    parser.add_argument("-o", "--output_base", default="output", help="Base output directory")
    parser.add_argument("--threads", type=int, default=64,
                         help="Threads for minimap2 and Medaka (default: 64)")
    parser.add_argument("--min_count", type=int, default=100,
                         help="Minimum best-hit read count for a reference to be accepted (default: 100)")
    parser.add_argument("--max_iter", type=int, default=10,
                         help="Maximum Medaka polishing iterations per sample (default: 10)")
    parser.add_argument("--pattern", default="*.fastq*",
                         help="Glob pattern (relative to --input_dir) selecting which files to "
                              "process, e.g. 'Mmus*.fastq*' (default: '*.fastq*')")
    parser.add_argument("--show_medaka_log", action="store_true",
                         help="Show Medaka's own stdout/stderr instead of suppressing it")
    parser.add_argument("--interactive", action="store_true",
                         help="If a sample's iterative filtering keeps zero reads, prompt "
                              "whether to retry with relaxed thresholds instead of skipping "
                              "the sample outright (default: off, for unattended batch runs)")
    parser.add_argument("-m", "--medaka_model", default=DEFAULT_MEDAKA_MODEL,
                         help=f"Medaka basecaller/flowcell model (default: '{DEFAULT_MEDAKA_MODEL}'; "
                              "use an R9.4.1 model, e.g. r941_min_hac_g507, for R9 data)")

    # Per-read alignment filters applied each iteration to build the FASTQ
    # fed into the next round.
    parser.add_argument("--min_mapq", type=int, default=20)
    parser.add_argument("--min_aligned_fraction", type=float, default=0.80)
    parser.add_argument("--max_softclip_fraction", type=float, default=0.15)
    parser.add_argument("--min_aligned_bases", type=int, default=200)
    parser.add_argument("--min_identity", type=float, default=0.85)
    parser.add_argument("--max_indel_size", type=int, default=None)
    parser.add_argument("--min_trim_depth", type=int, default=5,
                         help="Trim 5'/3' consensus ends below this read depth (default: 5)")
    args = parser.parse_args()

    missing = check_external_tools()
    if missing:
        sys.exit(
            f"Error: could not find the following external tools, please install "
            f"them first: {', '.join(missing)}\n"
            f"e.g.: conda install -c bioconda minimap2 samtools medaka"
        )

    input_path = Path(args.input_dir).resolve()
    output_base = Path(args.output_base).resolve()
    ref_db = Path(args.reference_db).resolve()

    output_base.mkdir(parents=True, exist_ok=True)

    fastq_files = sorted(input_path.glob(args.pattern))
    if not fastq_files:
        print(f"No FASTQ files matching pattern '{args.pattern}' found in {input_path}")
        return

    print(f"[*] Found {len(fastq_files)} samples. Starting pipeline...")

    for fastq in fastq_files:
        prefix = fastq.name.split('.')[0]

        sample_ref_fasta = output_base / f"{prefix}_ref.fasta"
        sample_out_prefix = output_base / prefix
        sample_log_path = output_base / f"{prefix}.log"

        with log_to_file(sample_log_path):
            print("\n" + "=" * 50)
            print(f"Processing sample: {prefix}")
            print("=" * 50)

            try:
                found = find_best_reference(
                    fastq, ref_db, sample_ref_fasta,
                    threads=args.threads, min_count=args.min_count,
                )
                if not found or not sample_ref_fasta.exists():
                    print(f"[!] Skipped {prefix}: no reference cleared --min_count {args.min_count} reads.")
                    continue

                print(f"[+] Found reference: {sample_ref_fasta}. Proceeding to iterative Medaka polishing...")
                final_fasta = run_medaka_iterative(
                    fastq, sample_ref_fasta, sample_out_prefix,
                    sample_name=prefix,
                    max_iter=args.max_iter, threads=args.threads,
                    quiet=(not args.show_medaka_log),
                    model=args.medaka_model,
                    min_mapq=args.min_mapq,
                    min_aligned_fraction=args.min_aligned_fraction,
                    max_softclip_fraction=args.max_softclip_fraction,
                    min_aligned_bases=args.min_aligned_bases,
                    min_identity=args.min_identity,
                    max_indel_size=args.max_indel_size,
                    min_trim_depth=args.min_trim_depth,
                    interactive=args.interactive,
                )

                if final_fasta and Path(final_fasta).exists():
                    print(f"[+] Final consensus ready: {final_fasta}")
                else:
                    print(f"[!] Warning: Medaka polishing did not produce a final consensus for {prefix}.")

            except PipelineError as e:
                print(f"[!] Skipped {prefix}: {e}")
                continue
            except Exception:
                # Anything unanticipated still shouldn't abort the rest of the batch.
                print(f"[!] Skipped {prefix} due to an unexpected error:")
                traceback.print_exc(file=sys.stdout)
                continue

    print("\n[Done] Batch process finished.")


if __name__ == "__main__":
    main()
