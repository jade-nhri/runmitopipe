# runmitopipe

A batch pipeline for reconstructing mitochondrial genomes from Nanopore long-read sequencing data: automatically selects the best-matching reference, iteratively polishes it with Medaka, and trims unreliable ends based on read depth.

## Pipeline overview

1. **Reference selection**: maps a sample's reads to a reference database with minimap2 and picks the single reference sequence with the most best-hit reads.
2. **Iterative Medaka polishing**: repeatedly polishes the reference with Medaka. Each round filters reads by alignment quality (mapping quality, aligned fraction, soft-clip fraction, sequence identity) and feeds the filtered set into the next round, until the consensus converges. Once converged, unreliable 5'/3' ends are trimmed based on read depth.
3. **Batch execution**: repeats the above two steps for every FASTQ file in the input directory.

## Installation

### Prerequisites

These tools depend on external programs that **pip cannot install automatically**. Please install the ones you need yourself first:

| Command | External tools needed |
|---|---|
| `runmitopipe` | [minimap2](https://github.com/lh3/minimap2), [samtools](http://www.htslib.org/), [medaka](https://github.com/nanoporetech/medaka) |
| `evaluate-primers` | [seqkit](https://github.com/shenwei356/seqkit) |
| `download-ncbi-fasta` | none (pure Python + `requests`/`tqdm`, installed automatically) |
| `build-taxonomy-summary` | none to run, but needs NCBI's `nucl_gb.accession2taxid(.gz)` and `taxdump.tar.gz` (data files, not tools -- see that script's `--help` for download links) |
| `merge-matrix-lineage` | none (pure Python + `pandas`/`openpyxl`, installed automatically) |
| Reproducing the manuscript analysis (below) | additionally needs [sra-tools](https://github.com/ncbi/sra-tools) (`fasterq-dump`) |

`medaka` has fairly complex dependencies (including tensorflow); installing everything via conda avoids most version conflicts:

```bash
conda create -n runmitopipe -c bioconda -c conda-forge minimap2 samtools medaka seqkit sra-tools python=3.10
conda activate runmitopipe
```

### Installing this package

```bash
git clone https://github.com/jade-nhri/runmitopipe.git
cd runmitopipe
pip install -e ".[dev]"
```

## Usage

```bash
runmitopipe \
    -i raw_fastq_dir \
    -r reference_db.fasta \
    -o output_base \
    --threads 64 --min_count 100 --max_iter 5
```

Output files (per sample, under `output_base/`):

| File | Description |
|---|---|
| `{prefix}_ref.fasta` | the reference selected for this sample |
| `{prefix}.consensus_untrimmed.fasta` | consensus before depth trimming |
| `{prefix}.final.fasta` | **the pipeline's final result** — use this file for downstream analysis / submission |
| `{prefix}.filtered_reads.fastq` | the (filtered) read set that produced the final consensus |
| `{prefix}.log` | the full processing log for this sample |

Common options:

- `--min_count`: minimum number of best-hit reads a candidate reference must have to be selected (default: 100)
- `--max_iter`: maximum number of Medaka polishing rounds (default: 10)
- `--min_trim_depth`: depth threshold for end trimming; ends below this depth are cut off (default: 5)
- `--interactive`: if a round's filtering keeps zero reads (common for samples where mitochondrial reads are a small fraction of the input, e.g. whole-genome sequencing data), whether to prompt for a retry with relaxed thresholds (default: off, suitable for unattended batch runs)
- `-m / --medaka_model`: the Medaka basecaller/flowcell model — make sure this matches the actual sequencing chemistry used (an R9 dataset needs the corresponding R9.4.1 model, otherwise Medaka won't warn you, it will just silently produce a lower-quality result)

Run `runmitopipe --help` for the full list of options.

## Additional tools

Besides the main `runmitopipe` pipeline, this package installs four supporting command-line tools used to build the reference database and reproduce the *in silico* primer-detection analysis from the manuscript. Each has its own `--help`; brief pointers below.

### `download-ncbi-fasta`

Batch-downloads FASTA sequences from NCBI for a list of accessions, automatically retrying anything still missing after the first pass.

```bash
download-ncbi-fasta -i accessions.txt -o sequences.fasta
```

### `evaluate-primers`

Runs `seqkit locate` for a panel of universal primers against a reference mitochondrial genome database and tabulates *in silico* primer detection (per-primer occurrence, per-genome hit matrix, overall detection rate).

```bash
evaluate-primers --ref sequences.fasta --primers primers.tsv --outdir results/
```

### `build-taxonomy-summary`

Resolves the NCBI taxonomic lineage (kingdom → species) for genomes listed in a `genome_primer_matrix.tsv` (from `evaluate-primers`) and summarizes counts by phylum and class.

```bash
build-taxonomy-summary \
    --matrix genome_primer_matrix.tsv \
    --acc2taxid nucl_gb.accession2taxid.gz \
    --nodes nodes.dmp --names names.dmp \
    --outdir taxonomy_out/ \
    --all-genomes
```

### `merge-matrix-lineage`

Left-joins the primer-detection matrix with the resolved lineage table and writes the combined table to Excel, matching the manuscript's Supplementary Table layout.

```bash
merge-matrix-lineage \
    --matrix genome_primer_matrix.tsv \
    --lineage taxonomy_out/genome_lineage.tsv \
    --out SupplementaryTables.xlsx \
    --sheet-name TableS1
```

See the "Reproducing the manuscript analysis" section below for how these four fit together end-to-end.

## Testing

```bash
pip install -e ".[dev]"   # only needed once, installs pytest
pytest tests/
```

This runs a small set of unit tests for a few helper functions used in read filtering (`get_softclip_lengths`, `get_max_indel_length`, etc.). They don't require minimap2/samtools/medaka to be installed and run in under a second, so there's no reason not to run them. They're not a substitute for testing on real sequencing data — they only catch regressions in the helper logic itself.

Run this after making any code change, before committing, to catch accidental regressions early.

## Reproducing the manuscript analysis (end-to-end integration test)

Unlike `pytest tests/` above, this is a slow, full end-to-end run against real data — reference database construction, real Nanopore reads from SRA, and the actual Medaka polishing pipeline. It exercises the whole tool chain, not just a few helper functions, but can take anywhere from minutes to hours depending on NCBI/SRA traffic and your machine. Run it after installing (`pip install -e ".[dev]"`) and the external tools in the Prerequisites table above, including `sra-tools`.

**1. Build a reference mitochondrial genome database.**
`tests/refseq_12680_acc.txt` and `tests/refseq_15658_acc.txt` are the accession lists used for the two reference datasets in the manuscript (12,680 and 15,658 accessions; ~150–190 KB as plain text, not sequence data, so they're checked into the repo directly):

```bash
download-ncbi-fasta -i tests/refseq_12680_acc.txt -o mito_12680.fasta
```

This downloads all 12,680 sequences from NCBI, retrying anything missing after the first pass (see `download-ncbi-fasta --help`). For a quicker smoke test instead of the full run, make a small subset first:

```bash
head -20 tests/refseq_12680_acc.txt > tests/acc_subset.txt
download-ncbi-fasta -i tests/acc_subset.txt -o mito_subset.fasta
```

**2. Download real Nanopore reads from SRA.**

```bash
mkdir SRR && cd SRR
fasterq-dump SRR40736960 --threads 64
fasterq-dump SRR40736961 --threads 64
cd ..
```

**3. Run the main pipeline against the downloaded reads and reference.**

```bash
runmitopipe \
    -i SRR/ \
    -r mito_12680.fasta \
    -o output \
    -m r941_min_hac_g507 \
    --min_count 50
```

`-m r941_min_hac_g507` selects an R9.4.1 Medaka model — match this to whatever flow cell/basecaller actually produced the reads you're testing with (see `runmitopipe --help`).

**4. Check the result.** For each sample under `SRR/`, `output/` should contain `{prefix}.final.fasta` (the reconstructed mitochondrial genome) and `{prefix}.log` (the full per-sample run log). A successful run ends with `[+] Final consensus ready: ...` for each sample in its log.

## Citation

<!-- TODO: add citation info once the paper is published, ideally alongside a Zenodo DOI -->

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See [LICENSE](LICENSE) for details.
