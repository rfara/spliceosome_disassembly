# Short-read RNA-seq workflow

This Snakemake workflow processes paired-end RNA-seq reads through:

1. paired-end trimming with a two-stage `cutadapt` pass
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
12. BAM indexing with `samtools index`
13. STAR-unmapped read diagnostics, including `FastQC` and a summary table

The workflow files and outputs live under `shortread/snake`.

Typical invocation from the repository root:

```bash
conda env create -f shortread/snake/envs/rnaseq_pipeline.yml
conda activate spliceosome-shortread
snakemake -s shortread/snake/Snakefile --cores 16 --use-conda
```

The workflow derives each `R2` path automatically by replacing `_R1_` with `_R2_` in the configured sample path.

The main parameter you may want to adjust before the first run is `star.sjdb_overhang` in `shortread/snake/config.yaml`, which should match read length minus one.

The trimming defaults are now biased towards rescuing low-quality `R2` tails. The first `cutadapt` pass trims the anchored `R2` 5' fixed sequence plus a wildcarded `R2` 3' adapter (`N{15}AGATCGGAAGAGCGTCGTGC`). A second pass then performs a single terminal `poly-G` cleanup on `R2` before the final minimum-length filter, which avoids repeatedly re-trimming `G` tails after other adapters have already been removed.

`R1` adapter trimming also uses a minimum overlap threshold (`trim.adapter_min_overlap`, default `6`) so short real terminal motifs such as `AG` are not treated as adapter sequence.

STAR itself does not provide a true mate-specific mismatch/score threshold for paired-end reads. The workflow therefore rescues damaged `R2` reads indirectly by running STAR in local mode and lowering the pair-level normalized score/match filters through `star.align_extra`.

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
The combined outputs include both the original 5' end metaprofile and a coverage-style metaprofile that assumes anchored 3' ends align to the 3' splice site.
Anchored fragments are required to have a fragment 3' end near the intron 3' splice site and a `read1` 5' end that still falls inside the selected intron. The stored intron-offset tables retain all such intronic 5' ends, even when they fall outside the plotted metaprofile window, so upstream starts contribute correctly to the coverage metaprofile.
Combined metaprofile plots show each condition as a thick mean trace with a shaded 95% confidence interval derived from replicate-to-replicate variability.

There is also an intron-level heterogeneity analysis that tests whether branchpoint-proximal reads are more unevenly distributed between introns than expected under a single shared branching probability. It reports this separately for each sample and pooled condition, and compares the configured query condition against the configured control condition using shared introns ranked by control branching.

The workflow also includes a sequence-context analysis for branchpoint-associated RT arrest. It extracts branchpoint-flanking sequence and the first intronic bases after the 5' splice site, then uses the configured control condition to estimate sequence-dependent RT readthrough versus arrest.

There is also a relative branching feature analysis that uses the configured control condition as a baseline predictor of per-intron branchiness, then asks which intron features and branchpoint-flanking sequence contexts explain branch enrichment or depletion in the configured query condition.

There is also a branchpoint readthrough-event analysis for reads that extend through the selected branchpoint. It keeps the same anchored fragment assignment, restricts to `read1` alignments whose aligned path covers the branchpoint and whose `read1` 5' end lies upstream of it, then uses the `read1` `MD` tag plus CIGAR to profile mismatches, deletions, and insertions separately relative to the branchpoint. The combined comparison is restricted to introns with at least `branchpoint_analysis.readthrough_shared_min_reads_all_samples` traversing reads in every sample, and the metaprofiles report event frequency as a fraction of traversing-read coverage at each offset. The readthrough plots can be cropped independently of the underlying quantification with `branchpoint_analysis.readthrough_plot_upstream` and `branchpoint_analysis.readthrough_plot_downstream`.
These readthrough metaprofile plots use the same visualization convention: a thick mean trace with a shaded 95% confidence interval for each condition.

The workflow also produces a filtered readthrough replot that blacklists recurrent high-indel introns, excludes loci with extreme total indel burden in any sample, and removes introns with recurrent single-offset indel spikes. The default blacklist is driven by `branchpoint_analysis.readthrough_blacklist_*`, `branchpoint_analysis.readthrough_blacklist_single_offset_*`, and `branchpoint_analysis.readthrough_max_total_indel_percent_any_sample` in `shortread/snake/config.yaml`.

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
- `shortread/snake/results/branchpoints/combined/branchpoint_coverage_metaprofile.png`
- `shortread/snake/results/branchpoints/readthrough_events/samples/{sample}.summary.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/samples/{sample}.metaprofile.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/samples/{sample}.site_counts.tsv.gz`
- `shortread/snake/results/branchpoints/readthrough_events/samples/{sample}.position_counts.tsv.gz`
- `shortread/snake/results/branchpoints/readthrough_events/combined/shared_introns.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/combined/metaprofile.by_sample.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/combined/metaprofile.by_condition.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/combined/summary.by_sample.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/combined/summary.by_condition.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/combined/branchpoint_readthrough_event_metaprofile.png`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/blacklist.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/filter_summary.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/shared_introns.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/metaprofile.by_sample.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/metaprofile.by_condition.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/summary.by_sample.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/summary.by_condition.tsv`
- `shortread/snake/results/branchpoints/readthrough_events/filtered/branchpoint_readthrough_event_metaprofile.filtered.png`
- `shortread/snake/results/branchpoints/heterogeneity/summary.by_sample.tsv`
- `shortread/snake/results/branchpoints/heterogeneity/summary.by_condition.tsv`
- `shortread/snake/results/branchpoints/heterogeneity/condition_comparison.summary.tsv`
- `shortread/snake/results/branchpoints/heterogeneity/branching_heterogeneity.png`
- `shortread/snake/results/branchpoints/sequence_context/summary.tsv`
- `shortread/snake/results/branchpoints/sequence_context/context_effects.tsv`
- `shortread/snake/results/branchpoints/sequence_context/top_bottom_enrichment.tsv`
- `shortread/snake/results/branchpoints/sequence_context/model_summary.tsv`
- `shortread/snake/results/branchpoints/sequence_context/branchpoint_sequence_context.png`
- `shortread/snake/results/branchpoints/relative_branching_features/summary.tsv`
- `shortread/snake/results/branchpoints/relative_branching_features/model_summary.tsv`
- `shortread/snake/results/branchpoints/relative_branching_features/numeric_feature_summary.tsv`
- `shortread/snake/results/branchpoints/relative_branching_features/relative_branching_features.png`
- `shortread/snake/results/branchpoints/relative_branching_features/feature_group_comparison.tsv`
- `shortread/snake/results/branchpoints/relative_branching_features/feature_distributions.png`
