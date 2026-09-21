# runmitopipe

A batch pipeline for reconstructing mitochondrial genomes from Nanopore long-read sequencing data: automatically selects the best-matching reference, iteratively polishes it with Medaka, and trims unreliable ends based on read depth.

## Pipeline overview

1. **Reference selection**: maps a sample's reads to a reference database with minimap2 and picks the single reference sequence with the most best-hit reads.
2. **Iterative Medaka polishing**: repeatedly polishes the reference with Medaka. Each round filters reads by alignment quality (mapping quality, aligned fraction, soft-clip fraction, sequence identity) and feeds the filtered set into the next round, until the consensus converges. Once converged, unreliable 5'/3' ends are trimmed based on read depth.
3. **Batch execution**: repeats the above two steps for every FASTQ file in the input directory.

## Installation

### Prerequisites

This tool depends on external programs that **pip cannot install automatically**. Please install these yourself first:

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](http://www.htslib.org/)
- [medaka](https://github.com/nanoporetech/medaka) (has fairly complex dependencies; installing via conda is recommended to avoid tensorflow version conflicts)

The easiest way is to install everything with conda in one go:

```bash
conda create -n runmitopipe -c bioconda -c conda-forge minimap2 samtools medaka python=3.10
conda activate runmitopipe
```

> Note: this project also has other related scripts (not yet packaged here) that additionally depend on [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [seqkit](https://github.com/shenwei356/seqkit). The `runmitopipe` command itself only needs minimap2/samtools/medaka. If those scripts get added to this package later, install the extra tools with:
>
> ```bash
> conda install -c bioconda blast seqkit
> ```

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

## Testing

```bash
pytest tests/
```

## Citation

<!-- TODO: add citation info once the paper is published, ideally alongside a Zenodo DOI -->

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See [LICENSE](LICENSE) for details.
