#!/usr/bin/env python3
"""
build_taxonomy_summary.py

Given genome_primer_matrix.tsv from evaluate_primers.py, resolves the full
taxonomic lineage of each in silico detectable genome (recovered == 1) and
summarizes counts by phylum and by class-within-phylum for a two-ring
donut chart (inner ring = phylum, outer ring = class).

Steps:
  1. Select genomes with detectable == 1 (in silico detectable genomes;
     use --all-genomes to process every genome in the matrix instead --
     this is what Supplementary Table S1/S2 actually do: every genome
     gets a lineage, detectable or not).
  2. Map each genome's accession to an NCBI taxid via
     nucl_gb.accession2taxid(.gz).
  3. Resolve the lineage for each taxid from the NCBI taxdump
     (nodes.dmp + names.dmp) at the standard ranks used in Supplementary
     Table S1/S2: kingdom, phylum, class, order, family, genus, species.
  4. Aggregate counts/percentages by phylum and by class-within-phylum.

Required input files (download once, reusable across projects):
  nucl_gb.accession2taxid.gz
      https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
  taxdump.tar.gz -> extract nodes.dmp and names.dmp
      https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
      tar xzf taxdump.tar.gz nodes.dmp names.dmp

Usage:
    python3 build_taxonomy_summary.py \
        --matrix genome_primer_matrix.tsv \
        --acc2taxid nucl_gb.accession2taxid.gz \
        --nodes nodes.dmp \
        --names names.dmp \
        --outdir taxonomy_out/ \
        --all-genomes

Outputs (in --outdir):
  genome_lineage.tsv    one row per genome: genome_id, taxid, and the
                         7 rank columns (blank if not assigned) -- when
                         joined onto genome_primer_matrix.tsv (see
                         merge_matrix_lineage.py), this reproduces the
                         column layout of Supplementary Table S1/S2.
  phylum_summary.tsv    phylum, count, percent            (inner ring)
  class_summary.tsv     phylum, class, count, percent     (outer ring)
  unmapped_genomes.txt  genome IDs with no taxid found in acc2taxid file
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict, Counter

# Matches the rank columns used in Supplementary Table S1/S2, in order.
MAJOR_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_target_genome_ids(matrix_path, all_genomes=False):
    ids = []
    with open(matrix_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if all_genomes or row.get("detectable", "0") == "1":
                ids.append(row["genome_id"])
    return ids


def strip_version(acc):
    return acc.split(".")[0]


def load_accession2taxid(path, wanted_ids):
    """
    Stream the (large) accession2taxid file once, keeping only rows whose
    accession or accession.version matches something in wanted_ids.
    Returns dict: genome_id (as given in wanted_ids) -> taxid (str)
    """
    wanted_with_version = set(wanted_ids)
    wanted_no_version = {strip_version(x): x for x in wanted_ids}

    result = {}
    print(f"[taxonomy] Scanning {path} for {len(wanted_ids)} accessions ...",
          file=sys.stderr)
    with open_maybe_gz(path) as fh:
        header = fh.readline()  # accession, accession.version, taxid, gi
        for line_num, line in enumerate(fh, start=2):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            accession, accession_version, taxid = parts[0], parts[1], parts[2]

            if accession_version in wanted_with_version:
                result[accession_version] = taxid
            elif accession in wanted_no_version:
                result[wanted_no_version[accession]] = taxid

            if line_num % 20_000_000 == 0:
                print(f"[taxonomy]   ... {line_num:,} lines scanned, "
                      f"{len(result)}/{len(wanted_ids)} resolved so far",
                      file=sys.stderr)

    print(f"[taxonomy] Resolved {len(result)}/{len(wanted_ids)} accessions to taxid.",
          file=sys.stderr)
    return result


def load_nodes(nodes_path):
    """Return parent_of[taxid] = parent_taxid, rank_of[taxid] = rank (both str keys)."""
    parent_of = {}
    rank_of = {}
    print(f"[taxonomy] Loading {nodes_path} ...", file=sys.stderr)
    with open_maybe_gz(nodes_path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            tax_id, parent_id, rank = fields[0], fields[1], fields[2]
            parent_of[tax_id] = parent_id
            rank_of[tax_id] = rank
    print(f"[taxonomy] Loaded {len(parent_of):,} nodes.", file=sys.stderr)
    return parent_of, rank_of


def load_names(names_path):
    """Return name_of[taxid] = scientific name (str keys)."""
    name_of = {}
    print(f"[taxonomy] Loading {names_path} ...", file=sys.stderr)
    with open_maybe_gz(names_path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            tax_id, name_txt, unique_name, name_class = (
                fields[0], fields[1], fields[2], fields[3]
            )
            if name_class == "scientific name":
                name_of[tax_id] = name_txt
    print(f"[taxonomy] Loaded {len(name_of):,} scientific names.", file=sys.stderr)
    return name_of


def resolve_lineage(taxid, parent_of, rank_of, name_of, max_steps=100):
    """
    Walk from taxid up to the root, recording the name at each of the
    MAJOR_RANKS ranks encountered. Returns dict rank -> name (missing
    ranks are simply absent from the dict).
    """
    lineage = {}
    current = taxid
    steps = 0
    while current in parent_of and steps < max_steps:
        rank = rank_of.get(current, "no rank")
        if rank in MAJOR_RANKS and rank not in lineage:
            lineage[rank] = name_of.get(current, "")
        parent = parent_of[current]
        if parent == current:
            break  # reached root (root's parent is itself)
        current = parent
        steps += 1
    return lineage


def main():
    ap = argparse.ArgumentParser(
        description="Resolve taxonomic lineage for in silico detectable genomes "
                    "and summarize phylum/class counts for Figure 1D."
    )
    ap.add_argument("--matrix", required=True,
                     help="genome_primer_matrix.tsv from evaluate_primers.py")
    ap.add_argument("--acc2taxid", required=True,
                     help="nucl_gb.accession2taxid(.gz) file (or a pre-filtered subset)")
    ap.add_argument("--nodes", required=True, help="NCBI taxdump nodes.dmp")
    ap.add_argument("--names", required=True, help="NCBI taxdump names.dmp")
    ap.add_argument("--outdir", default="taxonomy_out")
    ap.add_argument("--all-genomes", action="store_true",
                     help="Process every genome in the matrix, not just detectable==1 "
                          "(use this to reproduce Supplementary Table S1/S2, which give "
                          "every genome a lineage regardless of detection status)")
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Get target genome IDs
    genome_ids = load_target_genome_ids(args.matrix, all_genomes=args.all_genomes)
    subset_label = "all" if args.all_genomes else "in silico detectable (detectable==1)"
    print(f"[taxonomy] {len(genome_ids)} genomes selected ({subset_label}).",
          file=sys.stderr)

    # 2. Map accession -> taxid
    genome_to_taxid = load_accession2taxid(args.acc2taxid, genome_ids)
    unmapped = [g for g in genome_ids if g not in genome_to_taxid]
    if unmapped:
        with open(f"{args.outdir}/unmapped_genomes.txt", "w") as out:
            out.write("\n".join(unmapped) + "\n")
        print(f"[taxonomy] WARNING: {len(unmapped)} genomes had no taxid match "
              f"(see unmapped_genomes.txt).", file=sys.stderr)

    # 3. Load taxonomy dump
    parent_of, rank_of = load_nodes(args.nodes)
    name_of = load_names(args.names)

    # 4. Resolve lineage per genome
    lineage_path = f"{args.outdir}/genome_lineage.tsv"
    phylum_counter = Counter()
    class_counter = defaultdict(Counter)  # phylum -> Counter(class)

    with open(lineage_path, "w") as out:
        header = ["genome_id", "taxid"] + MAJOR_RANKS
        out.write("\t".join(header) + "\n")
        for gid in genome_ids:
            taxid = genome_to_taxid.get(gid, "")
            if not taxid:
                out.write(f"{gid}\t\t" + "\t".join([""] * len(MAJOR_RANKS)) + "\n")
                continue
            lineage = resolve_lineage(taxid, parent_of, rank_of, name_of)
            row = [gid, taxid] + [lineage.get(r, "") for r in MAJOR_RANKS]
            out.write("\t".join(row) + "\n")

            phylum = lineage.get("phylum", "Unclassified") or "Unclassified"
            cls = lineage.get("class", "Unclassified") or "Unclassified"
            phylum_counter[phylum] += 1
            class_counter[phylum][cls] += 1
    print(f"[taxonomy] Wrote {lineage_path}", file=sys.stderr)

    total = sum(phylum_counter.values())

    # 5. Phylum summary (inner ring)
    phylum_path = f"{args.outdir}/phylum_summary.tsv"
    with open(phylum_path, "w") as out:
        out.write("phylum\tcount\tpercent\n")
        for phylum, count in phylum_counter.most_common():
            pct = 100 * count / total if total else 0
            out.write(f"{phylum}\t{count}\t{pct:.2f}\n")
    print(f"[taxonomy] Wrote {phylum_path}", file=sys.stderr)

    # 6. Class-within-phylum summary (outer ring)
    class_path = f"{args.outdir}/class_summary.tsv"
    with open(class_path, "w") as out:
        out.write("phylum\tclass\tcount\tpercent\n")
        for phylum, count in phylum_counter.most_common():
            for cls, ccount in class_counter[phylum].most_common():
                pct = 100 * ccount / total if total else 0
                out.write(f"{phylum}\t{cls}\t{ccount}\t{pct:.2f}\n")
    print(f"[taxonomy] Wrote {class_path}", file=sys.stderr)

    print(f"[taxonomy] DONE. {total} genomes classified into "
          f"{len(phylum_counter)} phyla.", file=sys.stderr)


if __name__ == "__main__":
    main()
