# Short-read RNA-seq workflow

This Snakemake workflow processes paired-end RNA-seq reads through:

1. paired-end trimming with `cutadapt`
2. UMI extraction from `R1` while retaining synchronized `R2` reads with `umi_tools extract`
3. removal of the retained fixed `CGTGAT` sequence from `R1`
4. permissive paired pre-mapping to rRNA with `bowtie2`
5. permissive paired pre-mapping to snRNA with `bowtie2`
6. paired-end genomic alignment with `STAR`
7. UMI-aware paired-end deduplication with `umi_tools dedup`
8. `FastQC` on trimmed input FASTQs and `umi_ready` FASTQs
9. `samtools flagstat` and `samtools stats` on produced BAMs
10. a per-reference snRNA pre-map count table with mapped read and fragment counts
11. a per-sample RNA-content table combining `rRNA` pre-map, `snRNA` pre-map, and undeduplicated genomic reads
12. `Qualimap RNA-seq QC` on indexed genome-aligned and deduplicated BAMs
13. BAM indexing with `samtools index`
14. STAR-unmapped read diagnostics, including `FastQC` and a summary table

The workflow files and outputs live under `shortread/snake`.

Typical invocation from the repository root:

```bash
conda env create -f shortread/snake/envs/rnaseq_pipeline.yml
conda activate spliceosome-shortread
snakemake -s shortread/snake/Snakefile --cores 16 --use-conda
```

The workflow derives each `R2` path automatically by replacing `_R1_` with `_R2_` in the configured sample path.

The main parameter you may want to adjust before the first run is `star.sjdb_overhang` in `shortread/snake/config.yaml`, which should match read length minus one.

For SLURM submission with Snakemake 8, each rule exposes an `sbatch` argument string in `params.cluster`. A typical cluster invocation is:

```bash
snakemake -s shortread/snake/Snakefile \
  --jobs 50 \
  --cores 200 \
  --latency-wait 60 \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "sbatch {params.cluster}"
```

Add `--use-conda` if you want Snakemake to manage environments instead of using your already activated `spliceosome-shortread` environment.

The workflow currently maps short jobs to `c_short` and medium jobs to `c_medium` or `m_medium`, with time limits configured in `shortread/snake/config.yaml`.

The snRNA count table is written to `shortread/snake/results/qc/snrna_counts/{sample}.snrna_counts.tsv`. It reports one row per snRNA reference sequence, with both mapped read-segment counts and estimated fragment counts from primary `read1` alignments.

The RNA-content table is written to `shortread/snake/results/qc/rna_content/{sample}.rna_content.tsv`. It counts read pairs with this precedence:

1. `rRNA` from the rRNA pre-map
2. `snRNA` from the snRNA pre-map or genomic annotation
3. `mitochondrial`
4. small ncRNA subclasses (`snoRNA`, `scaRNA`, `miRNA`, `tRNA`, `ribozyme`, `vaultRNA`, `misc_RNA`)
5. `lncRNA`
6. `protein_coding_intron`
7. `protein_coding_exon`
8. `pseudogene`
9. `other_annotated_genic`
10. `intergenic`
11. `ambiguous`

This uses the undeduplicated STAR BAM, counts fragments rather than alignment records, ignores strand, and lets intronic small ncRNAs win over host-gene intron annotation.

To help diagnose STAR `too short` failures, the workflow also keeps STAR `Unmapped.out.mate1/2` files and writes:

- `shortread/snake/results/qc/fastqc/star_unmapped/...`
- `shortread/snake/results/qc/star_unmapped/{sample}.summary.tsv`

The summary table reports per-mate length quantiles plus the most common 12 nt prefixes and suffixes, which is usually enough to spot short inserts, residual adapters, or repetitive terminal sequence.
