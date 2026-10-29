#!/usr/bin/env python3
"""Compare independent Radiant searches with a shared-library batch (requires pyarrow)."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import time

import pyarrow.parquet as parquet


def run(binary, common, inputs, directory):
    directory.mkdir(parents=True, exist_ok=False)
    command = [str(binary), *common, *map(str, inputs), "--output-folder", str(directory)]
    started = time.monotonic()
    with (directory / "run.log").open("w") as log:
        result = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={**os.environ, "QT_HASH_SEED": "0", "QT_QPA_PLATFORM": "offscreen",
                 "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
        )
    text = (directory / "run.log").read_text()
    return {
        "command": command,
        "exit_code": result.returncode,
        "seconds": time.monotonic() - started,
        "library_loads": text.count("Finished reading library"),
        "completed": text.count("COMPLETED:"),
        "failed": text.count("FAILED:"),
        "log": str(directory / "run.log"),
    }


def same_value(left, right, absolute, relative):
    if isinstance(left, float) and isinstance(right, float):
        return (math.isnan(left) and math.isnan(right)) or math.isclose(
            left, right, abs_tol=absolute, rel_tol=relative
        )
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        return len(left) == len(right) and all(
            same_value(a, b, absolute, relative) for a, b in zip(left, right)
        )
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(
            same_value(left[key], right[key], absolute, relative) for key in left
        )
    return left == right


def compare(single, batch, absolute, relative):
    left = parquet.read_table(single)
    right = parquet.read_table(batch)
    if not left.schema.equals(right.schema):
        return {"passed": False, "reason": "Output schemas differ"}
    # Normalize row order by candidate identity, preserving duplicate rows.
    # Rounding is deliberately avoided: all numeric values are checked with the
    # documented tolerances below, including nested arrays.
    keys = [key for key in ("PeptideStringWithMods", "PeptideStringWithModsDecoyOrigin",
                           "Charge", "TargetKey", "ScanNumber", "IsDecoy")
            if key in left.column_names]
    if not keys:
        keys = [key for key in left.column_names if "peptide" in key.lower()
                or "charge" in key.lower() or "scantime" in key.lower()]
    left_rows, right_rows = left.to_pylist(), right.to_pylist()
    duplicate_identities = 0
    if keys:
        order = lambda row: tuple(str(row[key]) for key in keys)
        duplicate_identities = max(
            len(rows) - len({order(row) for row in rows}) for rows in (left_rows, right_rows)
        )
        left_rows.sort(key=order)
        right_rows.sort(key=order)
    different = {column: 0 for column in left.column_names}
    exactly_different = {column: 0 for column in left.column_names}
    maximum_absolute_difference = {}
    for a, b in zip(left_rows, right_rows):
        for column in left.column_names:
            different[column] += not same_value(a[column], b[column], absolute, relative)
            exactly_different[column] += not same_value(a[column], b[column], 0, 0)
            if isinstance(a[column], float) and isinstance(b[column], float):
                if math.isfinite(a[column]) and math.isfinite(b[column]):
                    maximum_absolute_difference[column] = max(
                        maximum_absolute_difference.get(column, 0), abs(a[column] - b[column])
                    )
    mismatches = {column: count for column, count in different.items() if count}
    return {
        "passed": len(left_rows) > 0 and len(left_rows) == len(right_rows) and not mismatches,
        "single_rows": len(left_rows),
        "batch_rows": len(right_rows),
        "columns": len(left.column_names),
        "different_columns": mismatches,
        "exactly_equal": len(left_rows) == len(right_rows) and not any(exactly_different.values()),
        "maximum_absolute_numeric_difference": max(maximum_absolute_difference.values(), default=0),
        "sort_keys": keys,
        "duplicate_identities": duplicate_identities,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "library", "fasta", "config", "output"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-6)
    parser.add_argument("--relative-tolerance", type=float, default=1e-6)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()
    if len(args.inputs) < 2:
        parser.error("Supply at least two input files")
    if len({path.name for path in args.inputs}) != len(args.inputs):
        parser.error("Inputs must have distinct filenames for the shared output directory")
    if any(not math.isfinite(value) or value < 0 for value in
           (args.absolute_tolerance, args.relative_tolerance)):
        parser.error("Tolerances must be finite and nonnegative")
    common = list(map(str, (args.library.resolve(), args.fasta.resolve(), args.config.resolve())))
    binary = args.binary.resolve()
    inputs = [path.resolve() for path in args.inputs]
    args.output.mkdir(parents=True, exist_ok=False)
    runs = []
    for index, path in enumerate(inputs):
        runs.append(run(binary, common, [path], args.output / f"single-{index + 1}"))
        print(f"Single {index + 1}: exit={runs[-1]['exit_code']}, "
              f"seconds={runs[-1]['seconds']:.2f}", flush=True)
    batch = run(binary, common, inputs, args.output / "batch")
    print(f"Batch: exit={batch['exit_code']}, seconds={batch['seconds']:.2f}", flush=True)
    comparisons = []
    for index, path in enumerate(inputs):
        name = path.name + ".radiantDIA"
        single_output = args.output / f"single-{index + 1}" / name
        batch_output = args.output / "batch" / name
        if not single_output.exists() or not batch_output.exists():
            comparisons.append({"input": str(path), "passed": False,
                                "reason": "Expected output file is missing"})
        else:
            comparisons.append({"input": str(path), **compare(
                single_output, batch_output, args.absolute_tolerance, args.relative_tolerance)})
    passed = all(result["exit_code"] == 0 and result["library_loads"] == 1 for result in runs)
    passed &= batch["exit_code"] == 0 and batch["library_loads"] == 1
    passed &= batch["completed"] == len(inputs) and batch["failed"] == 0
    passed &= all(result["passed"] for result in comparisons)
    summary = {
        "passed": passed,
        "absolute_tolerance": args.absolute_tolerance,
        "relative_tolerance": args.relative_tolerance,
        "single_runs": runs,
        "batch_run": batch,
        "comparisons": comparisons,
        "single_total_seconds": sum(result["seconds"] for result in runs),
        "batch_seconds": batch["seconds"],
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
