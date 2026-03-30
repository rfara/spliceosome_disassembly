#!/usr/bin/env python3

import argparse
import subprocess

import pysam


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        order = list(bam.references)
        counts = {
            ref: {"length": length, "read_count": 0, "fragment_count": 0}
            for ref, length in zip(bam.references, bam.lengths)
        }

    command = ["samtools", "collate", "-f", "-u", "-O", "-@", str(args.threads), args.bam]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    try:
        with pysam.AlignmentFile(process.stdout, "rb") as bam:
            current_name = None
            current_refs = set()
            for alignment in bam.fetch(until_eof=True):
                if alignment.is_unmapped:
                    continue

                ref = alignment.reference_name
                counts[ref]["read_count"] += 1

                if current_name is None:
                    current_name = alignment.query_name
                if alignment.query_name != current_name:
                    for previous_ref in current_refs:
                        counts[previous_ref]["fragment_count"] += 1
                    current_name = alignment.query_name
                    current_refs = set()
                current_refs.add(ref)

            if current_name is not None:
                for previous_ref in current_refs:
                    counts[previous_ref]["fragment_count"] += 1
    finally:
        if process.stdout is not None:
            process.stdout.close()
        return_code = process.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)

    with open(args.output, "w") as handle:
        handle.write("reference\tlength\tread_count\tfragment_count\n")
        for ref in order:
            entry = counts[ref]
            handle.write(
                f"{ref}\t{entry['length']}\t{entry['read_count']}\t{entry['fragment_count']}\n"
            )


if __name__ == "__main__":
    main()
