# Short-read Pipeline

This directory contains the Snakemake workflow for processing spliceosome disassembly short-read data and generating the R plots used by the project.

## Layout

- `Snakefile`: workflow entry point.
- `config.yaml`: sample paths, references, and analysis parameters.
- `scripts/`: Python scripts that generate processed tables from mapped, deduplicated reads.
- `r/`: R scripts that generate final plots.
- `envs/`: conda environment definition.
- `results/`: generated workflow outputs.

`results/processed_data` contains the table outputs used by the R plots:

- `premrna_mrna/`
- `branchpoint_metaprofiles/`
- `downstream_exon_branchpoint_metaprofiles/`
- `anchored_percentile_overlap/`

`results/plots` contains the final R-generated figures.
It also contains `plot_analysis_counts.tsv`, which summarizes the read and feature counts used for each plot, and `results/plots/source_data/`, which contains compact supplementary tables underlying each plot.

## Running Locally

With Snakemake and Conda available, run the workflow from the repository root:

```bash
snakemake -s shortread/Snakefile --cores 16 --use-conda
```

No manual data download is required. The first run uses the six SRA runs linked
from [GEO accession GSE329374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329374):

| Sample | SRA run |
|---|---|
| `ILS_1` | `SRR38300414` |
| `ILS_2` | `SRR38300413` |
| `ILS_3` | `SRR38300412` |
| `DIS_1` | `SRR38300411` |
| `DIS_2` | `SRR38300410` |
| `DIS_3` | `SRR38300409` |

For each sample, the workflow downloads and validates the SRA run, converts it
to paired FASTQs with `fasterq-dump`, and gzip-compresses the reads under
`shortread/reads/`. These files are ignored by Git and reused on later runs.
The temporary SRA and uncompressed FASTQs are removed after successful
conversion; nevertheless, the download jobs require substantial temporary disk
space.

If NCBI's accession resolver is temporarily unavailable, the workflow uses the
checksum-pinned NCBI [SRA Lite](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/)
URLs in [`config.yaml`](config.yaml). SRA Lite preserves the read sequences but
simplifies base-quality scores. The analysis does not quality-filter reads,
although raw-read FastQC quality distributions can differ when this fallback is
used.

The main sequencing parameter to check before a first run is `star.sjdb_overhang` in `shortread/config.yaml`; it should be read length minus one.

The workflow uses the rDNA FASTA, snRNA FASTA, GENCODE v44 primary-assembly GTF, and branchpoint table bundled under `annotation/`. The GRCh38 primary-assembly genome FASTA is downloaded by Snakemake from GENCODE release 44 and verified with the recorded MD5 checksum before STAR indexing.

The sample-to-run mapping and canonical raw-read paths are recorded in
[`config.yaml`](config.yaml). If the canonical FASTQs already exist, Snakemake
skips their download.

## Running On The Cluster

Submit the Snakemake controller job:

```bash
bash shortread/run_cluster.sh --use-conda
```

The cluster wrapper runs `shortread/Snakefile` from the repository root and writes controller logs to `shortread/logs/controller`.

## Active Analyses

The workflow keeps the analyses that feed the final R plots:

- deduplicated genome BAMs and BAM/QC summaries
- pre-mRNA/mRNA summary and 3' end metaprofiles
- branchpoint-centred and 3'SS-centred metaprofiles on one intron set included in analysis
- branchpoint metaprofiles for reads extending into the downstream exon
- anchored percentile overlap summaries
- final R plots in `results/plots`

Obsolete exploratory analyses have been removed, including snRNA subtype plots, canonical/non-canonical branchpoint metaprofiles, heterogeneity, relative branching features, query-branch-threshold results, residual logos, sequence-context analysis, control-enriched introns, readthrough-event plots, and Python-generated plots.
