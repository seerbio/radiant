#!/usr/bin/env python3
"""Plot direct RAW versus mzML identification overlaps from Radiant results.

Requires pyarrow, matplotlib, and matplotlib-venn. Reads the completed Parquet
search results, computes target-only sets using each level's own q-value, and
writes PNG/SVG/PDF figures plus JSON/CSV provenance and overlap counts.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib_venn import venn2, venn2_circles
import pyarrow as pa
import pyarrow.parquet as pq


RAW_COLOR = "#1976b8"
MZML_COLOR = "#d96912"
LEVELS = (
    ("precursor", "Precursors", "PrecursorQValue", ("PeptideStringWithMods", "Charge")),
    ("modified_peptide", "Modified peptides", "PeptideQValue", ("PeptideStringWithMods",)),
    ("protein_group", "Protein groups", "ProteinQValue", ("ProteinGroup",)),
)
COLUMNS = (
    "IsDecoy",
    "PeptideStringWithMods",
    "Charge",
    "ProteinGroup",
    "PrecursorQValue",
    "PeptideQValue",
    "ProteinQValue",
)


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def set_sha256(values):
    """Hash sorted JSON keys with one UTF-8 key per line."""
    digest = hashlib.sha256()
    for value in sorted(values):
        digest.update(json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def identification_key(row, level):
    peptide = row["PeptideStringWithMods"]
    if level in ("precursor", "modified_peptide"):
        if not isinstance(peptide, str) or not peptide.strip():
            raise ValueError("An accepted target has a missing peptide identifier")
        if level == "modified_peptide":
            return peptide
        charge = row["Charge"]
        if (
            charge is None
            or not math.isfinite(charge)
            or charge <= 0
            or int(charge) != charge
        ):
            raise ValueError("An accepted precursor has a missing or invalid charge")
        return (peptide, int(charge))
    group = row["ProteinGroup"]
    if not isinstance(group, str) or not group.strip():
        raise ValueError("An accepted target has a missing protein-group identifier")
    # Match IdLevelQValueAnnotator::cleanedProteinGroup in Radiant.
    return group.strip().replace(",", " ")


def load_sets(path, threshold):
    parquet = pq.ParquetFile(path)
    missing = set(COLUMNS) - set(parquet.schema_arrow.names)
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    table = parquet.read(columns=COLUMNS)
    rows = table.to_pylist()
    targets = [row for row in rows if row["IsDecoy"] == 0]
    sets = {}
    details = {}
    for level, _, q_column, _ in LEVELS:
        accepted = []
        missing_q = 0
        for row in targets:
            q_value = row[q_column]
            if q_value is None or not math.isfinite(q_value):
                missing_q += 1
                continue
            if not 0 <= q_value <= 1:
                raise ValueError(f"{path}: {q_column} contains a value outside [0, 1]")
            if q_value < threshold:
                accepted.append(row)
        values = {identification_key(row, level) for row in accepted}
        sets[level] = values
        details[level] = {
            "q_column": q_column,
            "passing_target_rows": len(accepted),
            "unique_identifications": len(values),
            "excluded_target_rows_with_missing_or_nonfinite_q": missing_q,
            "set_sha256": set_sha256(values),
        }
    provenance = {
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "result_rows": parquet.metadata.num_rows,
        "result_columns": len(parquet.schema_arrow),
        "target_rows": len(targets),
        "required_column_null_counts": {
            column: table[column].null_count for column in COLUMNS
        },
        "identifications": details,
    }
    return sets, provenance


def compare_sets(raw_sets, mzml_sets):
    comparisons = []
    for level, title, q_column, key_columns in LEVELS:
        raw = raw_sets[level]
        mzml = mzml_sets[level]
        shared = raw & mzml
        union = raw | mzml
        comparisons.append(
            {
                "level": level,
                "label": title,
                "q_column": q_column,
                "key_columns": list(key_columns),
                "raw_total": len(raw),
                "mzml_total": len(mzml),
                "raw_only": len(raw - mzml),
                "shared": len(shared),
                "mzml_only": len(mzml - raw),
                "union": len(union),
                "jaccard": len(shared) / len(union) if union else None,
                "identical_sets": raw == mzml,
                "raw_only_keys": sorted(raw - mzml),
                "mzml_only_keys": sorted(mzml - raw),
            }
        )
    return comparisons


def draw_panel(ax, comparison):
    title = comparison["label"]
    shared = comparison["shared"]
    if comparison["identical_sets"] and shared:
        # Identical nonempty sets occupy exactly the same circle.
        ax.add_patch(Circle((0, 0), 0.86, fill=False, color=RAW_COLOR, linewidth=3))
        ax.add_patch(
            Circle(
                (0, 0), 0.86, fill=False, color=MZML_COLOR,
                linewidth=2.2, linestyle=(0, (5, 3)),
            )
        )
        ax.text(0, 0.20, f"{shared:,}", ha="center", va="center", fontsize=22, weight="bold")
        ax.text(0, -0.10, "shared", ha="center", va="center", fontsize=12)
        ax.text(0, -0.39, "100% overlap", ha="center", va="center", fontsize=11)
        ax.text(-0.63, -1.14, "RAW only\n0", ha="center", va="top", fontsize=11, color=RAW_COLOR)
        ax.text(0.63, -1.14, "mzML only\n0", ha="center", va="top", fontsize=11, color=MZML_COLOR)
    elif not comparison["union"]:
        ax.text(0, 0, "No passing identifications", ha="center", va="center", fontsize=12)
        ax.text(0, -0.4, "Both sets empty", ha="center", fontsize=11)
    else:
        subsets = (comparison["raw_only"], comparison["mzml_only"], shared)
        diagram = venn2(subsets=subsets, set_labels=("", ""), ax=ax)
        for subset in ("10", "01", "11"):
            patch = diagram.get_patch_by_id(subset)
            if patch is not None:
                patch.set_facecolor("#edf2f6")
                patch.set_alpha(0.7)
            label = diagram.get_label_by_id(subset)
            if label is not None:
                label.set_fontsize(18 if subset != "11" else 22)
                label.set_fontweight("bold")
        circles = venn2_circles(subsets=subsets, ax=ax, linewidth=2.5)
        circles[0].set_color(RAW_COLOR)
        circles[1].set_color(MZML_COLOR)
        circles[1].set_linestyle((0, (5, 3)))
        overlap = comparison["jaccard"] * 100
        ax.text(0.5, 0, f"{overlap:.1f}% overlap (intersection / union)", transform=ax.transAxes,
                ha="center", fontsize=10)
    ax.set_title(title, fontsize=13, weight="bold", pad=8)
    ax.set_xlim(-1.20, 1.20)
    ax.set_ylim(-1.54, 1.08)
    ax.set_aspect("equal")
    ax.set_axis_off()


def plot(comparisons, output_dir, dataset_label, threshold):
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42}
    )
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 5.1))
    fig.subplots_adjust(left=0.03, right=0.97, top=0.76, bottom=0.25, wspace=0.11)
    fig.suptitle("Radiant search: direct RAW vs mzML", y=0.968, fontsize=16, weight="bold")
    fig.text(0.5, 0.906, dataset_label, ha="center", fontsize=11)
    fig.legend(
        handles=(
            Line2D([0], [0], color=RAW_COLOR, linewidth=3, label="Direct Thermo RAW"),
            Line2D([0], [0], color=MZML_COLOR, linewidth=2.2,
                   linestyle=(0, (5, 3)), label="Reference mzML"),
        ),
        loc="upper center", bbox_to_anchor=(0.5, 0.867), ncol=2,
        frameon=False, fontsize=11, handlelength=2.7, columnspacing=2.5,
    )
    for ax, comparison in zip(axes, comparisons):
        draw_panel(ax, comparison)
    captions = (
        f"Targets only (IsDecoy = 0), q < {threshold:g} using precursor, peptide and protein-group q-values, respectively.",
        "Precursors: modified sequence + charge; peptides: modified sequence; proteins: reported groups.",
        "Identical sets have coincident circles. Panels are scaled independently.",
    )
    for y, caption in zip((0.176, 0.128, 0.080), captions):
        fig.text(0.5, y, caption, ha="center", fontsize=9.5, color="#40464d")
    paths = []
    for extension in ("png", "svg", "pdf"):
        path = output_dir / f"raw-vs-mzml-venn.{extension}"
        fig.savefig(path, dpi=200, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mzml_result", type=Path, help="Reference mzML .radiantDIA Parquet result")
    parser.add_argument("raw_result", type=Path, help="Direct RAW .radiantDIA Parquet result")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for figures and reports")
    parser.add_argument("--q-threshold", type=float, default=0.01, help="Strict q-value cutoff (default: 0.01)")
    parser.add_argument("--dataset-label", default="HeLa 250 pg · replicate 4", help="Figure dataset caption")
    args = parser.parse_args()
    if not math.isfinite(args.q_threshold) or not 0 < args.q_threshold <= 1:
        parser.error("--q-threshold must be finite and in (0, 1]")
    try:
        mzml_sets, mzml_provenance = load_sets(args.mzml_result, args.q_threshold)
        raw_sets, raw_provenance = load_sets(args.raw_result, args.q_threshold)
        comparisons = compare_sets(raw_sets, mzml_sets)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        figures = plot(comparisons, args.output_dir, args.dataset_label, args.q_threshold)
        report = {
            "dataset": args.dataset_label,
            "filter": {"IsDecoy": 0, "q_operator": "<", "q_threshold": args.q_threshold},
            "protein_group_normalization": "Trim surrounding whitespace; replace commas with spaces, matching Radiant.",
            "set_sha256_format": "Sorted identifier keys serialized as compact UTF-8 JSON, each followed by LF.",
            "byte_identical_results": mzml_provenance["sha256"] == raw_provenance["sha256"],
            "mzml_result": mzml_provenance,
            "raw_result": raw_provenance,
            "comparisons": comparisons,
            "figures": [path.name for path in figures],
        }
        json_path = args.output_dir / "raw-vs-mzml-overlap.json"
        json_path.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
        csv_path = args.output_dir / "raw-vs-mzml-overlap.csv"
        fields = (
            "level", "q_column", "raw_total", "mzml_total", "raw_only", "shared",
            "mzml_only", "union", "jaccard", "identical_sets",
        )
        with csv_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(comparisons)
    except (OSError, ValueError, pa.ArrowException) as error:
        parser.error(str(error))
    for comparison in comparisons:
        print(
            f"{comparison['label']}: RAW only {comparison['raw_only']:,}; "
            f"shared {comparison['shared']:,}; mzML only {comparison['mzml_only']:,}"
        )
    print(f"Figures and reports: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
