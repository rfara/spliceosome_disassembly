# Short-read RNA-seq workflow

This Snakemake workflow processes single-end RNA-seq reads through:

1. adapter and homopolymer trimming with `cutadapt`
2. UMI extraction with `umi_tools extract`
3. permissive local pre-mapping to rRNA with `bowtie2`
4. permissive local pre-mapping to snRNA with `bowtie2`
5. end-to-end genomic alignment with `STAR`
6. UMI-aware deduplication with `umi_tools dedup`
7. BAM indexing with `samtools index`

The workflow files and outputs live under `shortread/snake`.

Typical invocation from the repository root:

```bash
conda env create -f shortread/snake/envs/rnaseq_pipeline.yml
conda activate spliceosome-shortread
snakemake -s shortread/snake/Snakefile --cores 16 --use-conda
```

The main parameter you may want to adjust before the first run is `star.sjdb_overhang` in `shortread/snake/config.yaml`, which should match read length minus one.

For SLURM submission with Snakemake 8, each rule exposes an `sbatch` argument string in `params.cluster`. A typical cluster invocation is:

```bash
snakemake -s shortread/snake/Snakefile \
  --use-conda \
  --jobs 50 \
  --latency-wait 60 \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "sbatch {params.cluster}"
```

The workflow currently maps short jobs to `c_short` or `m_short` and medium jobs to `c_medium` or `m_medium`, with time limits configured in `shortread/snake/config.yaml`.
