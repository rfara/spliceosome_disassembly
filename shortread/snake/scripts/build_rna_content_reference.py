#!/usr/bin/env python3

import argparse
import pickle

from classify_rna_content import parse_gtf


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    reference = parse_gtf(args.gtf)
    with open(args.output, "wb") as handle:
        pickle.dump(reference, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    main()
