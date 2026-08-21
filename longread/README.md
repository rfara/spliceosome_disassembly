# Long-read analysis

This directory contains the publication workflow for the single-replicate
TFIP11 and DHX35 Oxford Nanopore libraries. It is deliberately limited to the
two paper panels:

- branchpoint-centered deletion rate;
- read coverage around the 3′ splice site.

The earlier mapping-composition analyses and example-gene plots are not part of
this workflow.

## Inputs

The workflow retrieves the GEO-deposited, demultiplexed and adapter-trimmed
FASTQs automatically. The sample mapping is:

| Sample | GEO sample | SRA run | Local FASTQ |
|---|---|---|---|
| `tfip11` | `GSM9980085` | `SRR40278264` | `longread/reads/tfip11.fastq.gz` |
| `dhx35` | `GSM9980086` | `SRR40278263` | `longread/reads/dhx35.fastq.gz` |

The runs belong to
[GEO accession GSE329374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329374).
Snakemake validates each downloaded SRA run, converts it to FASTQ with
`fasterq-dump`, and gzip-compresses it. These large files are ignored by Git and
reused on later runs. Sample and barcode assignments are in
[`samples.tsv`](samples.tsv); SRA accessions and raw-read paths are in
[`config.yaml`](config.yaml). The pre-GEO basecalling, demultiplexing, and
trimming commands remain recorded in
[`UPSTREAM_PROCESSING.md`](UPSTREAM_PROCESSING.md).

## Running the workflow

From the repository root:

```bash
snakemake --snakefile longread/Snakefile --cores 16 --use-conda
```

The workflow first obtains the raw reads from the GEO-linked SRA runs. It also
downloads the frozen GENCODE v44 transcriptome and GRCh38 primary assembly after
checking their MD5 sums. Raw-read paths, SRA accessions, local reference paths,
and resource settings can be changed in [`config.yaml`](config.yaml). SRA
conversion creates an uncompressed temporary FASTQ, so allow substantially more
free space than the final compressed inputs occupy.

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
