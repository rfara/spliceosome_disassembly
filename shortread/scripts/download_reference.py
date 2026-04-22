#!/usr/bin/env python3

import argparse
import gzip
import hashlib
import os
import shutil
import tempfile
from pathlib import Path
from urllib.request import urlopen


def parse_args():
    parser = argparse.ArgumentParser(description="Download and optionally decompress a reference file.")
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--compressed-md5")
    parser.add_argument("--decompress-gzip", action="store_true")
    return parser.parse_args()


def download(url, destination):
    md5 = hashlib.md5()
    with urlopen(url) as response, open(destination, "wb") as handle:
        while True:
            chunk = response.read(8 * 1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
            md5.update(chunk)
    return md5.hexdigest()


def main():
    args = parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=output.parent) as tmpdir:
        tmpdir_path = Path(tmpdir)
        downloaded = tmpdir_path / "download"
        observed_md5 = download(args.url, downloaded)

        if args.compressed_md5 and observed_md5 != args.compressed_md5:
            raise SystemExit(
                f"MD5 mismatch for {args.url}: expected {args.compressed_md5}, observed {observed_md5}"
            )

        prepared = tmpdir_path / "prepared"
        if args.decompress_gzip:
            with gzip.open(downloaded, "rb") as source, open(prepared, "wb") as target:
                shutil.copyfileobj(source, target)
        else:
            os.replace(downloaded, prepared)

        os.replace(prepared, output)


if __name__ == "__main__":
    main()
