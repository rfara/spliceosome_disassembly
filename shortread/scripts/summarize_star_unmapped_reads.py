#!/usr/bin/env python3

import collections
import sys


def summarize_fastq(path, mate_label):
    length_counts = collections.Counter()
    prefix_counts = collections.Counter()
    suffix_counts = collections.Counter()
    total_reads = 0
    total_bases = 0
    reads_le_20 = 0
    reads_le_30 = 0
    reads_le_40 = 0

    with open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip("\n")
            handle.readline()
            handle.readline()
            read_len = len(seq)
            total_reads += 1
            total_bases += read_len
            length_counts[read_len] += 1
            reads_le_20 += int(read_len <= 20)
            reads_le_30 += int(read_len <= 30)
            reads_le_40 += int(read_len <= 40)
            prefix_counts[seq[:12]] += 1
            suffix_counts[seq[-12:]] += 1

    if total_reads == 0:
        return {
            "mate": mate_label,
            "total_reads": 0,
            "mean_length": 0,
            "min_length": 0,
            "p10_length": 0,
            "p25_length": 0,
            "median_length": 0,
            "p75_length": 0,
            "p90_length": 0,
            "max_length": 0,
            "reads_le_20": 0,
            "reads_le_30": 0,
            "reads_le_40": 0,
            "top_prefix_12": "",
            "top_prefix_12_count": 0,
            "top_suffix_12": "",
            "top_suffix_12_count": 0,
        }

    ordered_lengths = sorted(length_counts.items())

    def quantile(target):
        cutoff = total_reads * target
        seen = 0
        for read_len, count in ordered_lengths:
            seen += count
            if seen >= cutoff:
                return read_len
        return ordered_lengths[-1][0]

    top_prefix, top_prefix_count = prefix_counts.most_common(1)[0]
    top_suffix, top_suffix_count = suffix_counts.most_common(1)[0]

    return {
        "mate": mate_label,
        "total_reads": total_reads,
        "mean_length": round(total_bases / total_reads, 2),
        "min_length": ordered_lengths[0][0],
        "p10_length": quantile(0.10),
        "p25_length": quantile(0.25),
        "median_length": quantile(0.50),
        "p75_length": quantile(0.75),
        "p90_length": quantile(0.90),
        "max_length": ordered_lengths[-1][0],
        "reads_le_20": reads_le_20,
        "reads_le_30": reads_le_30,
        "reads_le_40": reads_le_40,
        "top_prefix_12": top_prefix,
        "top_prefix_12_count": top_prefix_count,
        "top_suffix_12": top_suffix,
        "top_suffix_12_count": top_suffix_count,
    }


def main():
    r1_path, r2_path, output_path = sys.argv[1:4]
    rows = [
        summarize_fastq(r1_path, "R1"),
        summarize_fastq(r2_path, "R2"),
    ]

    fields = [
        "mate",
        "total_reads",
        "mean_length",
        "min_length",
        "p10_length",
        "p25_length",
        "median_length",
        "p75_length",
        "p90_length",
        "max_length",
        "reads_le_20",
        "reads_le_30",
        "reads_le_40",
        "top_prefix_12",
        "top_prefix_12_count",
        "top_suffix_12",
        "top_suffix_12_count",
    ]

    with open(output_path, "w") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[field]) for field in fields) + "\n")


if __name__ == "__main__":
    main()
