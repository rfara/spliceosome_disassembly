#!/usr/bin/env python3

import argparse
import gzip
import hashlib
from http.client import IncompleteRead
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen


def parse_args():
    parser = argparse.ArgumentParser(description="Download and optionally decompress a file.")
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--compressed-md5")
    parser.add_argument("--decompress-gzip", action="store_true")
    return parser.parse_args()


def md5sum(path):
    md5 = hashlib.md5()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            md5.update(chunk)
    return md5.hexdigest()


def download_once(url, destination):
    offset = destination.stat().st_size if destination.exists() else 0
    request = Request(url)
    if offset:
        request.add_header("Range", f"bytes={offset}-")

    with urlopen(request, timeout=60) as response:
        append = offset > 0 and getattr(response, "status", None) == 206
        mode = "ab" if append else "wb"
        expected = response.headers.get("Content-Length")
        received = 0
        with open(destination, mode) as handle:
            while True:
                chunk = response.read(8 * 1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                received += len(chunk)

        if expected is not None and received != int(expected):
            raise OSError(
                f"Incomplete download from {url}: expected {expected} bytes, received {received}"
            )


def download(url, destination, attempts=5):
    for attempt in range(1, attempts + 1):
        try:
            download_once(url, destination)
            return md5sum(destination)
        except (IncompleteRead, OSError, URLError) as error:
            if attempt == attempts:
                raise
            delay = min(2**attempt, 30)
            print(
                f"Download attempt {attempt}/{attempts} failed: {error}; retrying in {delay}s",
                file=sys.stderr,
            )
            time.sleep(delay)


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
