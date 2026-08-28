#!/usr/bin/env python3
"""Losslessly retain only columns written by a sampling-only CmdStan schema."""

import argparse
import csv
from pathlib import Path


def read_header(path):
    with path.open(newline="") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                return next(csv.reader([line]))
    raise ValueError(f"No CmdStan header in {path}")


def strip_one(source, destination, wanted, provenance):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open(newline="") as input_handle, destination.open(
        "w", newline=""
    ) as output_handle:
        writer = csv.writer(output_handle, lineterminator="\n")
        header = None
        indices = None
        retained_rows = 0
        for line in input_handle:
            if header is None and line.startswith("#"):
                output_handle.write(line)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            row = next(csv.reader([line]))
            if header is None:
                header = row
                missing = [column for column in wanted if column not in header]
                if missing:
                    raise ValueError(
                        f"{source} lacks slim-schema columns: {', '.join(missing)}"
                    )
                indices = [header.index(column) for column in wanted]
                output_handle.write(
                    f"# streamlined_source = {source.resolve()}\n"
                    f"# streamlined_schema = {provenance}\n"
                    "# streamlined_operation = lossless_column_subset\n"
                )
                writer.writerow(wanted)
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Malformed CmdStan row in {source}: {len(row)} != {len(header)}"
                )
            writer.writerow([row[index] for index in indices])
            retained_rows += 1
    if retained_rows == 0:
        raise ValueError(f"No posterior rows retained from {source}")
    return retained_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--schema", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("sources", nargs="+", type=Path)
    args = parser.parse_args()

    if not args.schema.is_file():
        parser.error(f"Schema CSV does not exist: {args.schema}")
    if any(not source.is_file() for source in args.sources):
        parser.error("Every source CmdStan CSV must exist")
    wanted = read_header(args.schema)
    summaries = []
    for chain, source in enumerate(args.sources, 1):
        destination = args.output_dir / f"{args.output_prefix}-chain{chain}-1.csv"
        rows = strip_one(source, destination, wanted, str(args.schema.resolve()))
        summaries.append((chain, source, destination, rows, destination.stat().st_size))

    manifest = args.output_dir / f"{args.output_prefix}-strip-manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chain", "source", "destination", "rows", "bytes"])
        writer.writerows(summaries)


if __name__ == "__main__":
    main()
