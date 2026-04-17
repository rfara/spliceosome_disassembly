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

## Running Locally

Create the conda environment:

```bash
conda env create -f shortread/envs/rnaseq_pipeline.yml
```

Run the workflow:

```bash
snakemake -s shortread/Snakefile --cores 16 --use-conda
```

The main sequencing parameter to check before a first run is `star.sjdb_overhang` in `shortread/config.yaml`; it should be read length minus one.

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
- branchpoint-centred metaprofiles using the unfiltered shared intron set
- branchpoint metaprofiles for reads extending into the downstream exon
- anchored percentile overlap summaries
- final R plots in `results/plots`

Obsolete exploratory analyses have been removed, including snRNA subtype plots, canonical/non-canonical branchpoint metaprofiles, heterogeneity, relative branching features, query-branch-threshold results, residual logos, sequence-context analysis, control-enriched introns, readthrough-event plots, and Python-generated plots.
