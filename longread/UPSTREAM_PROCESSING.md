# Long-read processing before the deposited FASTQs

The GEO submission starts from the demultiplexed, adapter-trimmed FASTQ files in
`longread/reads/`; basecalling, demultiplexing, and trimming are therefore
documented here but are intentionally outside the Snakemake workflow.

The following commands are transcribed from the cluster scripts in
`20250129_tfip_dhx35`. Dorado 0.9.1 was used. The Cutadapt version was not
recorded in the original scripts and should not be inferred.

## Basecalling

The run used a MinION flow cell and Dorado's `hac` model on CUDA 12.3:

```bash
module load cuda/12.3.0
/groups/plaschka/rupert.faraway/software/dorado-0.9.1-linux-x64/bin/dorado \
  basecaller hac \
  20250129_1507_MN36858_FAW87778_d46613ea/pod5/ \
  -v -x cuda:all > dorado_out/basecalled.bam
```

The Dorado model cache associated with the run identifies the resolved HAC
model as `dna_r10.4.1_e8.2_400bps_hac@v5.0.0`.

## Demultiplexing

```bash
mkdir -p demux_fastq
cutadapt \
  -g CACGACGCTCTTCCGATCTNNNNNGGG \
  --revcomp --minimum-length 18 \
  dorado_out/basecalled.bam \
| cutadapt \
  -a bc1=AAAAAAGATTCAAGATCGGAAGAGCACACGTCTGA \
  -a bc2=AAAAAATGCTAGAGATCGGAAGAGCACACGTCTGA \
  --revcomp --minimum-length 18 \
  -o 'demux_fastq/{name}.fastq' -
```

Barcode assignments were:

| Barcode | Sample |
|---|---|
| `bc1` | TFIP11 |
| `bc2` | DHX35 |

The unassigned `unknown.fastq` file is not part of the paper analysis or GEO
submission.

## Poly(A)-adapter trimming

Each demultiplexed FASTQ was trimmed independently:

```bash
mkdir -p atrim_fastq
for input in demux_fastq/*.fastq; do
  name=$(basename "$input")
  cutadapt \
    -a AAAAAAAAAAAAAAAAAAAAAAAAAAA \
    --revcomp --minimum-length 18 \
    -o "atrim_fastq/${name}.fastq" \
    "$input"
done
```

The two resulting analysis inputs were renamed on import to the repository as
`longread/reads/tfip11.fastq.gz` and `longread/reads/dhx35.fastq.gz`.
