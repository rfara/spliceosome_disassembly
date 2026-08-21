# spliceosome_disassembly

This repository contains analysis scripts for RNA sequencing experiments from spliceosomes undergoing disassembly. The results have been published as part of Boreikaite et al. (2026):

Paper reference goes here.

The raw sequencing data are public under [GEO accession GSE329374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329374). The repository contains independent Snakemake workflows for the short-read and long-read datasets. Each workflow now resolves the GEO-linked SRA run accessions, downloads and validates the runs, creates the required FASTQs, downloads the large genome references, and performs the analysis.

After cloning the repository, run either workflow from the repository root with Snakemake and Conda available:

```bash
snakemake --snakefile shortread/Snakefile --cores 16 --use-conda
snakemake --snakefile longread/Snakefile --cores 16 --use-conda
```

The workflows are independent, so these commands may be run separately. The downloaded FASTQs are stored under `shortread/reads/` and `longread/reads/`, are ignored by Git, and are reused on subsequent runs. SRA conversion temporarily requires substantially more free disk space than the final gzip-compressed FASTQs.

See [`shortread/README.md`](shortread/README.md) and [`longread/README.md`](longread/README.md) for workflow-specific details and outputs.
