#!/usr/bin/env python3
"""Stream the small policy parameter subset out of large CmdStan CSVs."""

import argparse
import csv
from pathlib import Path


BASE_COLUMNS = [
    "beta_intercept", "beta_ink_effect", "beta_bracelet_effect",
    "wtp_value_utility", "hyper_wtp_mu", "base_mu_rep", "raw_u_sd.1",
    "dist_beta_v.1",
]
BASE_COLUMNS += [f"hyper_beta_1ord.{i}" for i in range(1, 5)]
BASE_COLUMNS += [f"hyper_dist_beta_1ord.{i}" for i in range(1, 5)]


def header_and_rows(path):
    handle = path.open(newline="")
    for line in handle:
        if not line.startswith("#"):
            header = next(csv.reader([line]))
            rows = csv.reader(
                line for line in handle if line.strip() and not line.startswith("#")
            )
            return handle, header, rows
    handle.close()
    raise ValueError(f"No CmdStan header in {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--include-lambda", choices=["none", "grouped", "arm"], default="none")
    parser.add_argument("--include-cluster-shock", type=int, default=0,
                        help="Number of cluster shock elements to retain")
    parser.add_argument("--include-asymmetric", action="store_true",
                        help="Retain the 48 asymmetric-observability parameters")
    parser.add_argument(
        "--asymmetric-structure",
        choices=["none", "full", "f1", "f2", "f3", "u3"],
        default="none",
        help="Retain the parameters present in an observability-ladder fit",
    )
    parser.add_argument("fits", nargs="+")
    args = parser.parse_args()

    wanted = list(BASE_COLUMNS)
    if args.include_lambda == "grouped":
        wanted += ["core_lambda_group_log_ratio_raw.1"]
    elif args.include_lambda == "arm":
        wanted += [f"core_lambda_arm_log_ratio_raw.{i}" for i in range(1, 4)]
    if args.include_cluster_shock:
        wanted += ["core_cluster_shock_sd.1"]
        wanted += [f"core_cluster_shock_raw.{i}" for i in range(1, args.include_cluster_shock + 1)]
    if args.include_asymmetric and args.asymmetric_structure != "none":
        parser.error("Use either --include-asymmetric or --asymmetric-structure")
    asymmetric_structure = (
        "full" if args.include_asymmetric else args.asymmetric_structure
    )
    if asymmetric_structure in ("full", "f1", "u3"):
        wanted += [f"core_recognition_intercept.{i}" for i in range(1, 3)]
    if asymmetric_structure in ("full", "f1"):
        wanted += [f"core_recognition_dist_slope.{i}" for i in range(1, 3)]
        wanted += [f"core_recognition_arm_intercept_raw.{i}.{j}"
                   for j in range(1, 4) for i in range(1, 3)]
        wanted += [f"core_recognition_arm_dist_raw.{i}.{j}"
                   for j in range(1, 4) for i in range(1, 3)]
    if asymmetric_structure in ("full", "f1", "f2"):
        wanted += [f"core_report_intercept.{i}.{j}"
                   for j in range(1, 3) for i in range(1, 3)]
        wanted += [f"core_report_dist_slope.{i}.{j}"
                   for j in range(1, 3) for i in range(1, 3)]
        wanted += [f"core_report_arm_intercept_raw.{i}.{j}"
                   for j in range(1, 7) for i in range(1, 3)]
    if asymmetric_structure == "full":
        wanted += [f"core_report_arm_dist_raw.{i}.{j}"
                   for j in range(1, 7) for i in range(1, 3)]
    if asymmetric_structure in ("f3", "u3"):
        wanted += [f"core_definite_intercept.{i}" for i in range(1, 3)]
        wanted += [f"core_definite_dist_slope.{i}" for i in range(1, 3)]
        wanted += [f"core_definite_arm_intercept_raw.{i}.{j}"
                   for j in range(1, 4) for i in range(1, 3)]
        wanted += ["core_definite_public_signal_dist_slope.1"]
        wanted += [f"core_accuracy_intercept.{i}" for i in range(1, 3)]
        wanted += [f"core_accuracy_arm_intercept_raw.{i}.{j}"
                   for j in range(1, 4) for i in range(1, 3)]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as output_handle:
        writer = csv.writer(output_handle)
        writer.writerow(["chain", "iteration", "source_csv", *wanted])
        for chain, name in enumerate(args.fits, 1):
            path = Path(name)
            handle, header, rows = header_and_rows(path)
            try:
                missing = sorted(set(wanted) - set(header))
                if missing:
                    raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")
                indices = [header.index(column) for column in wanted]
                for iteration, row in enumerate(rows, 1):
                    writer.writerow([chain, iteration, str(path), *[row[i] for i in indices]])
            finally:
                handle.close()


if __name__ == "__main__":
    main()
