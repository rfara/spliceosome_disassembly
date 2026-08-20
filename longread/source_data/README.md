# Publication source data

These version-controlled tables are generated from the mapped TFIP11 and DHX35
reads by `../scripts/summarize_metaprofiles.py` and are the numerical source of
the two long-read paper figures:

- `three_prime_splice_site_coverage.tsv`: normalized read coverage at positions
  -100 to +99 relative to the 3′ splice site;
- `branchpoint_deletion_metaprofile.tsv`: deletion rate per read at positions
  -60 to +60 relative to the predicted branchpoint;
- `branchpoint_deletion_auc.tsv`: the median background and background-subtracted
  sum from -4 to +2 shown on the deletion figure.

`protein` contains the sample name (`tfip11` or `dhx35`) and
`strand_corrected_position` is in nucleotides after orienting all features in
the direction of transcription. `score` is the read-coverage count at a
position. Deletion counts are normalized first by the number of qualifying
read–branchpoint overlaps and then by local normalized read coverage.

The AUC table records the median background and background-subtracted peak sums
used for the plot labels: 0.147465 for TFIP11 and 0.029153 for DHX35. Running
Snakemake regenerates these tables before plotting them.
