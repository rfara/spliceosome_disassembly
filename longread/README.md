# Long-read analysis

This directory contains the publication workflow for the single-replicate
TFIP11 and DHX35 Oxford Nanopore libraries. It is deliberately limited to the
two paper panels:

- branchpoint-centered deletion rate;
- read coverage around the 3′ splice site.

The earlier mapping-composition analyses and example-gene plots are not part of
this workflow.

## Inputs

Place the GEO-deposited, demultiplexed and adapter-trimmed FASTQs at:

```text
longread/reads/tfip11.fastq.gz
longread/reads/dhx35.fastq.gz
```

These large files are ignored by Git. Sample and barcode assignments are in
[`samples.tsv`](samples.tsv). The pre-GEO basecalling, demultiplexing, and
trimming commands are recorded in [`UPSTREAM_PROCESSING.md`](UPSTREAM_PROCESSING.md).

## Running the workflow

From the repository root:

```bash
snakemake --snakefile longread/Snakefile --cores 16 --use-conda
```

The workflow downloads the frozen GENCODE v44 transcriptome and GRCh38 primary
assembly after checking their MD5 sums. Local reference paths and resource
settings can be changed in [`config.yaml`](config.yaml).

The mapping order reproduces the cluster analysis:

```text
trimmed FASTQ
  → remove rRNA alignments
  → remove snRNA alignments
  → remove protein-coding transcriptome alignments
  → align the remaining reads to GRCh38
  → summarize metaprofiles
  → render the paper plots
```

Minimap2 uses `--MD -ax map-ont`. The analysis retains primary,
non-supplementary alignments with mapping quality strictly greater than 20,
matching the exploratory R code.

## Outputs

The workflow writes the version-controlled paper source data to
[`source_data/`](source_data/):

- `three_prime_splice_site_coverage.tsv`;
- `branchpoint_deletion_metaprofile.tsv`;
- `branchpoint_deletion_auc.tsv`.

The final figures are written to `longread/results/plots/` directly from these
tables. Analysis counts are written to
`longread/results/source_data/analysis_counts.tsv`.

Intermediate unmapped FASTQs, indexes, BAMs, logs, and benchmarks are ignored by
Git.

Since each condition has one long-read library, the plots are descriptive and
the workflow performs no replicate-based inferential statistics.
