#!/usr/bin/env python3
"""Compare native Radiant result values exactly, independently of row order.

Requires pyarrow and numpy. No numeric tolerance is used. Matching null values
and matching NaNs compare equal; Parquet encoding and metadata are not results.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


CANDIDATE_COLUMNS = (
    "PeptideStringWithMods",
    "Charge",
    "IsDecoy",
    "TargetKey",
    "ScanNumber",
)
INPUT_IDENTITY_COLUMNS = {
    "FileName",
    "FilePath",
    "InputFile",
    "InputFileName",
    "InputFilePath",
    "RawFile",
    "RawFileName",
    "RawFilePath",
    "SourceFile",
    "SourceFileName",
    "SourceFilePath",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_value(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")
    return value


def identification_counts(table):
    counts = {"rows": table.num_rows, "columns": table.num_columns}
    if "IsDecoy" not in table.column_names:
        return counts
    target = pc.fill_null(pc.equal(table["IsDecoy"], 0), False)
    counts["target_rows"] = pc.sum(pc.cast(target, pa.int64())).as_py() or 0
    counts["decoy_rows"] = table.num_rows - counts["target_rows"]
    for column in ("QValue", "PrecursorQValue", "PeptideQValue", "ProteinQValue"):
        if column not in table.column_names:
            continue
        accepted = pc.and_(target, pc.fill_null(pc.less(table[column], 0.01), False))
        selected = table.filter(accepted)
        result = {"target_rows_below_1_percent": selected.num_rows}
        for field, label in (
            ("PeptideStringWithMods", "unique_peptides_below_1_percent"),
            ("ProteinGroup", "unique_protein_groups_below_1_percent"),
        ):
            if field in selected.column_names:
                result[label] = len(pc.unique(selected[field]).drop_null())
        if all(field in selected.column_names for field in ("PeptideStringWithMods", "Charge")):
            result["unique_precursors_below_1_percent"] = len(
                set(zip(selected["PeptideStringWithMods"].to_pylist(), selected["Charge"].to_pylist()))
            )
        counts[column] = result
    return counts


def canonicalize(table, columns):
    table = table.select(columns).combine_chunks()
    if not columns or table.num_rows == 0:
        return table
    # Candidate identity comes first; every remaining field breaks ties so
    # repeated candidates are compared as a multiset, preserving multiplicity.
    keys = [field for field in CANDIDATE_COLUMNS if field in columns]
    keys.extend(field for field in columns if field not in keys)
    indices = pc.sort_indices(table, sort_keys=[(field, "ascending") for field in keys])
    return table.take(indices)


def matching_values(left, right):
    same = pc.or_(
        pc.fill_null(pc.equal(left, right), False),
        pc.and_(pc.is_null(left), pc.is_null(right)),
    )
    if pa.types.is_floating(left.type):
        same = pc.or_(
            same,
            pc.fill_null(pc.and_(pc.is_nan(left), pc.is_nan(right)), False),
        )
    return same


def compare_results(reference_path, raw_path, ignored_columns=(), example_count=3):
    reference = pq.read_table(reference_path)
    raw = pq.read_table(raw_path)
    ignored_columns = sorted(set(ignored_columns))
    forbidden = set(ignored_columns) - INPUT_IDENTITY_COLUMNS
    if forbidden:
        raise ValueError("Only input identity columns can be ignored: " + ", ".join(sorted(forbidden)))

    reference_names = set(reference.column_names) - set(ignored_columns)
    raw_names = set(raw.column_names) - set(ignored_columns)
    common = sorted(reference_names & raw_names)
    type_mismatches = {
        field: {"mzml": str(reference[field].type), "raw": str(raw[field].type)}
        for field in common
        if reference[field].type != raw[field].type
    }
    columns = [field for field in common if field not in type_mismatches]
    left = canonicalize(reference, columns)
    right = canonicalize(raw, columns)
    overlap = min(left.num_rows, right.num_rows)
    mismatches = {}
    different_rows = np.zeros(overlap, dtype=bool)

    for field in columns:
        left_values = left[field].slice(0, overlap)
        right_values = right[field].slice(0, overlap)
        same = matching_values(left_values, right_values)
        different = np.logical_not(same.to_numpy(zero_copy_only=False))
        count = int(different.sum())
        if not count:
            continue
        different_rows |= different
        entry = {"different_values": count, "examples": []}
        for index in np.flatnonzero(different)[:example_count]:
            index = int(index)
            entry["examples"].append(
                {
                    "sorted_row": index,
                    "mzml_candidate": {
                        key: json_value(left[key][index].as_py())
                        for key in CANDIDATE_COLUMNS
                        if key in columns
                    },
                    "raw_candidate": {
                        key: json_value(right[key][index].as_py())
                        for key in CANDIDATE_COLUMNS
                        if key in columns
                    },
                    "mzml_value": json_value(left_values[index].as_py()),
                    "raw_value": json_value(right_values[index].as_py()),
                }
            )
        if pa.types.is_floating(left_values.type) or pa.types.is_integer(left_values.type):
            left_numeric = left_values.to_numpy(zero_copy_only=False).astype(np.float64)
            right_numeric = right_values.to_numpy(zero_copy_only=False).astype(np.float64)
            with np.errstate(invalid="ignore", over="ignore"):
                deltas = np.abs(left_numeric[different] - right_numeric[different])
            known = deltas[~np.isnan(deltas)]
            entry["max_absolute_delta"] = json_value(float(known.max())) if known.size else None
        mismatches[field] = entry

    identical = (
        reference_names == raw_names
        and not type_mismatches
        and left.num_rows == right.num_rows
        and not mismatches
    )
    return {
        "identical": identical,
        "comparison": "Exact values in every compared field; row order ignored; matching nulls and NaNs equal.",
        "numeric_tolerance": 0,
        "mzml_result": {
            "path": str(reference_path.resolve()),
            "sha256": sha256(reference_path),
            "identifications": identification_counts(reference),
        },
        "raw_result": {
            "path": str(raw_path.resolve()),
            "sha256": sha256(raw_path),
            "identifications": identification_counts(raw),
        },
        "ignored_input_identity_columns": ignored_columns,
        "compared_columns": columns,
        "columns_only_in_mzml": sorted(reference_names - raw_names),
        "columns_only_in_raw": sorted(raw_names - reference_names),
        "column_type_mismatches": type_mismatches,
        "row_count_difference_raw_minus_mzml": right.num_rows - left.num_rows,
        "aligned_rows_compared": overlap,
        "aligned_rows_with_differences": int(different_rows.sum()),
        "different_column_count": len(mismatches),
        "column_mismatches": mismatches,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mzml_result", type=Path, help="Reference mzML search .radiantDIA Parquet file")
    parser.add_argument("raw_result", type=Path, help="Direct RAW search .radiantDIA Parquet file")
    parser.add_argument("--output", type=Path, help="Write the JSON comparison report to this path")
    parser.add_argument(
        "--ignore-input-column",
        action="append",
        default=[],
        help="Ignore an explicitly named input file identity field (no scientific fields allowed)",
    )
    parser.add_argument("--examples", type=int, default=3, help="Maximum example differences per column")
    args = parser.parse_args()
    if args.examples < 0:
        parser.error("--examples must be nonnegative")
    try:
        report = compare_results(
            args.mzml_result, args.raw_result, args.ignore_input_column, args.examples
        )
    except (OSError, ValueError, pa.ArrowException) as error:
        parser.error(str(error))
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")
    return 0 if report["identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
