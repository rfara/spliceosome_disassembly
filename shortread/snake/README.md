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

For SLURM submission with Snakemake 8, each rule exposes an `sbatch` argument string in `params.cluster`. It is safer to submit the Snakemake controller itself through `sbatch` as well, so a dropped interactive session does not kill the workflow. Use the helper script:

```bash
bash shortread/snake/run_cluster.sh
```

Run that from the repository root. Any extra arguments are passed straight through to `snakemake`, for example:

```bash
bash shortread/snake/run_cluster.sh --use-conda
bash shortread/snake/run_cluster.sh shortread/snake/results/branchpoints/combined/branchpoint_5prime_metaprofile.png
```

The wrapper defaults to a controller job on `c_medium` with `1` CPU, `8G` RAM, `2-00:00:00` walltime, `50` Snakemake jobs, and `200` workflow cores. You can override these with environment variables such as `CONTROLLER_PARTITION`, `CONTROLLER_MEM`, `CONTROLLER_TIME`, `SNAKEMAKE_JOBS`, `SNAKEMAKE_CORES`, and `CONDA_ENV_NAME`.

The workflow currently maps short jobs to `c_short` and medium jobs to `c_medium` or `m_medium`, with time limits configured in `shortread/snake/config.yaml`.

The contaminant pre-maps now discard a pair if either mate aligns to the pre-map reference. Only pairs where both mates remain unmapped are passed forward to the next stage.

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

## Branchpoint 5' end analysis

The workflow now includes a branchpoint-centred analysis built from the deduplicated genome BAMs in `shortread/snake/results/dedup`.

It:

1. derives protein-coding introns from MANE Select transcripts in `annotation/gencode.v44.primary_assembly.basic.annotation.gtf.gz`
2. assigns branchpoints from `annotation/colaseq_branchpoints_scores.txt`
3. selects the highest-scoring branchpoint per intron
4. keeps uniquely mapped, primary, non-supplementary, proper-pair fragments
5. anchors fragments whose RNA 3' end falls within `+/- 5 nt` of the intron 3' splice site
6. profiles the `read1` 5' end around the selected branchpoint
7. restricts the combined comparison to introns with at least `5` anchored fragments in every sample
8. reports the primary signal as the fraction of anchored intron-end fragments, with CPM retained as secondary context

Positive offsets in the metaprofile point towards the intron 3' splice site.
The per-sample summary tables report branchpoint-proximal counts at offset `0`, offset `+1`, and combined `0/+1`, both as fractions of anchored fragments and as CPM.
The threshold for the shared intron set is configurable via `branchpoint_analysis.shared_min_reads_all_samples` in `shortread/snake/config.yaml`.

You can run just this analysis from the repository root with:

```bash
snakemake -s shortread/snake/Snakefile \
  --cores 4 \
  --use-conda \
  shortread/snake/results/branchpoints/combined/branchpoint_5prime_metaprofile.png
```

Main outputs:

- `shortread/snake/results/branchpoints/reference/mane_select_top_branchpoints.tsv`
- `shortread/snake/results/branchpoints/samples/{sample}.summary.tsv`
- `shortread/snake/results/branchpoints/samples/{sample}.metaprofile.tsv`
- `shortread/snake/results/branchpoints/samples/{sample}.site_counts.tsv.gz`
- `shortread/snake/results/branchpoints/samples/{sample}.offset_counts.tsv.gz`
- `shortread/snake/results/branchpoints/combined/shared_introns.tsv`
- `shortread/snake/results/branchpoints/combined/metaprofile.by_sample.tsv`
- `shortread/snake/results/branchpoints/combined/metaprofile.by_condition.tsv`
- `shortread/snake/results/branchpoints/combined/summary.by_sample.tsv`
- `shortread/snake/results/branchpoints/combined/summary.by_condition.tsv`
- `shortread/snake/results/branchpoints/combined/branchpoint_5prime_metaprofile.png`
