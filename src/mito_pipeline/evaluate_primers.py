#!/usr/bin/env python3
"""
evaluate_primers.py

Runs `seqkit locate` for a panel of universal mitochondrial primers against
a reference multi-FASTA of complete animal mitochondrial genomes, and
tabulates in silico primer detection.

Detection criterion: a genome is considered in silico detectable when at
least two distinct primers are matched, with at least one primer-binding
site on each DNA strand (both thresholds are configurable).

seqkit locate note (verified against seqkit v2.8.2): by default it already
searches both strands of every reference sequence and reports the matched
strand in its "strand" column, so a single `seqkit locate -p SEQ` call per
primer is enough -- no need to run it twice.

Outputs (in --outdir):
    primer_occurrence_summary.tsv   per primer: how many genomes it
                                     matched (either strand)
    genome_primer_matrix.tsv        per genome: which primers hit it and
                                     on which strand(s), plus the overall
                                     detection call -- this is the table
                                     used to build Supplementary Table
                                     S1/S2 (see build_taxonomy_summary.py
                                     and merge_matrix_lineage.py)
    overall_detection_summary.txt   headline in silico detection rate
    raw_locate/<primer>.tsv         raw seqkit locate output per primer

genome_primer_matrix.tsv columns: genome_id, one column per primer (each
cell blank / "+" / "-" / "+/-"), n_primers_hit, has_plus_strand,
has_minus_strand, detectable (1 = meets the detection criterion above).

Usage:
    python evaluate_primers.py \
        --ref refseq_mito_sequences.fasta \
        --primers primers.tsv \
        --outdir results/ \
        --threads 8

primers.tsv format (tab-separated, NO header), one primer per line:
    Primer1<TAB>ACTGGGATTA<TAB>+
    Primer2<TAB>GCTGTTATCC<TAB>-
The third column (designed strand, from your primer design table) is
optional and for reporting only -- it does not affect the search, since
seqkit locate always searches both strands of the reference.
"""

import argparse
import subprocess
import sys
import os
import csv
from collections import defaultdict


def run_seqkit_locate(seqkit_bin, ref_fasta, primer_seq, out_path, threads):
    cmd = [
        seqkit_bin, "locate", ref_fasta,
        "-p", primer_seq,
        "-o", out_path,
        "-j", str(threads),
    ]
    subprocess.run(cmd, check=True)


def parse_locate_output(path):
    """Return dict: seqID -> set of strands ('+'/'-') matched by this primer."""
    hits = defaultdict(set)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            hits[row["seqID"]].add(row["strand"])
    return hits


def get_all_seq_ids(ref_fasta):
    """Read just the FASTA headers (fast, no need to load sequences)."""
    ids = []
    opener = open
    if ref_fasta.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(ref_fasta, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    ap = argparse.ArgumentParser(
        description="Tabulate in silico universal-primer detection across a "
                     "reference mtDNA database using seqkit locate."
    )
    ap.add_argument("--ref", required=True,
                     help="Reference multi-FASTA of complete mitochondrial genomes")
    ap.add_argument("--primers", required=True,
                     help="TSV file: name<TAB>sequence<TAB>designed_strand(optional)")
    ap.add_argument("--outdir", default="evaluate_primers_out", help="Output directory")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seqkit", default="seqkit", help="Path to seqkit binary")
    ap.add_argument("--min-primers", type=int, default=2,
                     help="Minimum number of distinct primer-binding sites for a "
                          "genome to be called in silico detectable (default: 2)")
    ap.add_argument("--require-both-strands", dest="require_both_strands",
                     action="store_true", default=True,
                     help="Require at least one hit on each strand for a genome "
                          "to be called in silico detectable (default: on)")
    ap.add_argument("--no-require-both-strands", dest="require_both_strands",
                     action="store_false")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    raw_dir = os.path.join(args.outdir, "raw_locate")
    os.makedirs(raw_dir, exist_ok=True)

    # --- Load primer list ---
    primers = []
    with open(args.primers) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            name, seq = parts[0].strip(), parts[1].strip()
            designed_strand = parts[2].strip() if len(parts) > 2 else "NA"
            primers.append((name, seq, designed_strand))

    if not primers:
        sys.exit("[evaluate_primers] ERROR: no primers loaded from --primers file.")

    # --- Get full genome ID list (denominator for percentages) ---
    print(f"[evaluate_primers] Reading sequence IDs from {args.ref} ...", file=sys.stderr)
    all_ids = get_all_seq_ids(args.ref)
    total_genomes = len(all_ids)
    print(f"[evaluate_primers] Reference contains {total_genomes} sequences.", file=sys.stderr)

    # genome_id -> {primer_name -> set(strands)}
    genome_matrix = defaultdict(lambda: defaultdict(set))
    primer_occurrence = {}

    for name, seq, designed_strand in primers:
        out_path = os.path.join(raw_dir, f"{name}.tsv")
        print(f"[evaluate_primers] seqkit locate: {name} ({seq}) ...", file=sys.stderr)
        run_seqkit_locate(args.seqkit, args.ref, seq, out_path, args.threads)
        hits = parse_locate_output(out_path)
        primer_occurrence[name] = len(hits)  # distinct genomes with >=1 hit, either strand
        for gid, strands in hits.items():
            genome_matrix[gid][name] |= strands
        print(f"[evaluate_primers]   -> {len(hits)} / {total_genomes} genomes matched", file=sys.stderr)

    # --- Output 1: per-primer occurrence summary ---
    summary_path = os.path.join(args.outdir, "primer_occurrence_summary.tsv")
    with open(summary_path, "w") as out:
        out.write("Primer\tSequence\tLength\tDesigned_Strand\tOccurrence\tOccurrence_pct\n")
        for name, seq, designed_strand in primers:
            occ = primer_occurrence[name]
            pct = 100 * occ / total_genomes if total_genomes else 0
            out.write(f"{name}\t{seq}\t{len(seq)}\t{designed_strand}\t{occ}\t{pct:.2f}\n")
    print(f"[evaluate_primers] Wrote {summary_path}", file=sys.stderr)

    # --- Output 2: per-genome hit matrix (-> Supplementary Table S1/S2) ---
    matrix_path = os.path.join(args.outdir, "genome_primer_matrix.tsv")
    primer_names = [p[0] for p in primers]
    n_detectable = 0
    with open(matrix_path, "w") as out:
        header = ["genome_id"] + primer_names + \
            ["n_primers_hit", "has_plus_strand", "has_minus_strand", "detectable"]
        out.write("\t".join(header) + "\n")
        for gid in all_ids:
            row_hits = genome_matrix.get(gid, {})
            cells = []
            strands_seen = set()
            n_primers_hit = 0
            for pname in primer_names:
                strands = row_hits.get(pname, set())
                if strands:
                    n_primers_hit += 1
                    strands_seen |= strands
                    cells.append("/".join(sorted(strands)))
                else:
                    cells.append("")
            has_plus = "+" in strands_seen
            has_minus = "-" in strands_seen
            detectable = (n_primers_hit >= args.min_primers) and (
                (has_plus and has_minus) if args.require_both_strands else True
            )
            if detectable:
                n_detectable += 1
            out.write(
                gid + "\t" + "\t".join(cells) + "\t"
                + f"{n_primers_hit}\t{int(has_plus)}\t{int(has_minus)}\t{int(detectable)}\n"
            )
    print(f"[evaluate_primers] Wrote {matrix_path}", file=sys.stderr)

    # --- Output 3: headline in silico detection rate ---
    detection_summary_path = os.path.join(args.outdir, "overall_detection_summary.txt")
    pct_detectable = 100 * n_detectable / total_genomes if total_genomes else 0
    strand_clause = ", with at least one primer-binding site on each strand" if args.require_both_strands else ""
    with open(detection_summary_path, "w") as out:
        out.write(f"Total genomes in reference: {total_genomes}\n")
        out.write(
            f"Genomes in silico detectable (>= {args.min_primers} distinct primers matched"
            f"{strand_clause}): {n_detectable} ({pct_detectable:.2f}%)\n"
        )
    print(f"[evaluate_primers] Wrote {detection_summary_path}", file=sys.stderr)
    print(
        f"[evaluate_primers] RESULT: {n_detectable}/{total_genomes} genomes "
        f"({pct_detectable:.2f}%) in silico detectable.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
