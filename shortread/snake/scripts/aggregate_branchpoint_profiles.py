#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONDITION_COLORS = {
    "ILS": "#1f77b4",
    "DIS": "#d95f02",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metaprofile", action="append", dest="metaprofiles", required=True)
    parser.add_argument("--summary", action="append", dest="summaries", required=True)
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--intron-offsets", action="append", dest="intron_offsets", required=True)
    parser.add_argument("--shared-min-reads", type=int, default=0)
    parser.add_argument("--output-metaprofile-by-sample", required=True)
    parser.add_argument("--output-metaprofile-by-condition", required=True)
    parser.add_argument("--output-summary-by-sample", required=True)
    parser.add_argument("--output-summary-by-condition", required=True)
    parser.add_argument("--output-shared-introns", required=True)
    parser.add_argument("--output-plot-png", required=True)
    parser.add_argument("--output-plot-pdf", required=True)
    return parser.parse_args()


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_tsv_rows(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, rows, fieldnames):
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def float_mean(values):
    return 0.0 if not values else sum(values) / len(values)


def float_sd(values):
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def float_sem(values):
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values) / math.sqrt(len(values))


def infer_sample_name(rows, path):
    samples = {row["sample"] for row in rows if "sample" in row and row["sample"]}
    if not samples:
        return Path(path).name.split(".")[0]
    if len(samples) != 1:
        raise ValueError(f"Expected one sample in {path}, found {sorted(samples)}")
    return next(iter(samples))


def count_value(row, field):
    raw_value = row.get(field, 0)
    if raw_value in {"", None}:
        return 0
    return int(float(raw_value))


def summarise_condition_rows(summary_rows, condition_order):
    grouped = defaultdict(list)
    for row in summary_rows:
        grouped[row["condition"]].append(row)

    numeric_fields = [field for field in summary_rows[0] if field not in {"sample", "condition"}]
    rows = []
    for condition in condition_order:
        entries = grouped[condition]
        condition_row = {
            "condition": condition,
            "replicate_count": len(entries),
        }
        for field in numeric_fields:
            values = [float(entry[field]) for entry in entries]
            condition_row[f"mean_{field}"] = float_mean(values)
            condition_row[f"sd_{field}"] = float_sd(values)
        rows.append(condition_row)
    return rows


def summarise_condition_profiles(metaprofile_rows, condition_order):
    grouped = defaultdict(lambda: defaultdict(list))
    for row in metaprofile_rows:
        grouped[row["condition"]][int(row["offset_nt"])].append(row)

    condition_rows = []
    for condition in condition_order:
        offsets = sorted(grouped[condition])
        for offset in offsets:
            entries = grouped[condition][offset]
            cpm_values = [float(entry["cpm"]) for entry in entries]
            anchored_fraction_values = [float(entry["anchored_fraction"]) for entry in entries]
            condition_rows.append(
                {
                    "condition": condition,
                    "offset_nt": offset,
                    "replicate_count": len(entries),
                    "mean_cpm": float_mean(cpm_values),
                    "sd_cpm": float_sd(cpm_values),
                    "sem_cpm": float_sem(cpm_values),
                    "mean_anchored_fraction": float_mean(anchored_fraction_values),
                    "sd_anchored_fraction": float_sd(anchored_fraction_values),
                    "mean_anchored_percent": float_mean(anchored_fraction_values) * 100.0,
                    "sd_anchored_percent": float_sd(anchored_fraction_values) * 100.0,
                }
            )
    return condition_rows


def build_shared_intron_set(site_counts_by_sample, sample_order, shared_min_reads):
    qualifying_sets = []
    for sample in sample_order:
        sample_rows = site_counts_by_sample.get(sample, {})
        qualifying_sets.append(
            {
                intron_id
                for intron_id, row in sample_rows.items()
                if count_value(row, "anchored_fragments") >= shared_min_reads
            }
        )
    if not qualifying_sets:
        return set()
    return set.intersection(*qualifying_sets)


def build_shared_introns_rows(shared_introns, site_counts_by_sample, site_metadata, sample_order):
    fieldnames = [
        "intron_id",
        "gene_id",
        "gene_name",
        "transcript_id",
        "intron_number",
        "chrom",
        "strand",
        "intron_start",
        "intron_end",
        "three_prime_ss",
        "branchpoint_position",
        "branchpoint_score",
        "branchpoint_to_3ss_nt",
        "branchpoint_candidates",
        "min_anchored_fragments_all_samples",
    ] + [f"{sample}_anchored_fragments" for sample in sample_order]

    rows = []
    for intron_id in sorted(shared_introns, key=lambda key: (site_metadata[key]["gene_name"], key)):
        metadata = site_metadata[intron_id]
        per_sample_counts = [
            count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "anchored_fragments")
            for sample in sample_order
        ]
        row = {
            "intron_id": intron_id,
            "gene_id": metadata["gene_id"],
            "gene_name": metadata["gene_name"],
            "transcript_id": metadata["transcript_id"],
            "intron_number": metadata["intron_number"],
            "chrom": metadata["chrom"],
            "strand": metadata["strand"],
            "intron_start": metadata["intron_start"],
            "intron_end": metadata["intron_end"],
            "three_prime_ss": metadata["three_prime_ss"],
            "branchpoint_position": metadata["branchpoint_position"],
            "branchpoint_score": metadata["branchpoint_score"],
            "branchpoint_to_3ss_nt": metadata["branchpoint_to_3ss_nt"],
            "branchpoint_candidates": metadata["branchpoint_candidates"],
            "min_anchored_fragments_all_samples": min(per_sample_counts) if per_sample_counts else 0,
        }
        for sample, count in zip(sample_order, per_sample_counts):
            row[f"{sample}_anchored_fragments"] = count
        rows.append(row)
    return rows, fieldnames


def aggregate_offset_counts(intron_offset_counts, shared_introns):
    total_counts = Counter()
    for intron_id in shared_introns:
        for offset, read_count in intron_offset_counts.get(intron_id, {}).items():
            total_counts[offset] += read_count
    return total_counts


def build_sample_metaprofile_rows(
    sample,
    condition,
    library_fragments,
    anchored_fragments,
    offset_range,
    total_offset_counts,
):
    rows = []
    for offset in offset_range:
        read_count = total_offset_counts[offset]
        anchored_fraction = 0.0 if anchored_fragments == 0 else read_count / anchored_fragments
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "offset_nt": offset,
                "read_count": read_count,
                "cpm": 0.0 if library_fragments == 0 else (read_count * 1_000_000.0 / library_fragments),
                "anchored_fraction": anchored_fraction,
                "anchored_percent": anchored_fraction * 100.0,
            }
        )
    return rows


def build_sample_summary_row(
    raw_summary_row,
    site_counts,
    total_offset_counts,
    shared_introns,
    shared_min_reads,
):
    summary_row = dict(raw_summary_row)
    library_fragments = count_value(raw_summary_row, "library_fragments")
    raw_anchored_fragments = count_value(raw_summary_row, "anchored_fragments")
    raw_anchored_introns = count_value(raw_summary_row, "anchored_introns_with_reads")

    eligible_rows = [site_counts[intron_id] for intron_id in shared_introns if intron_id in site_counts]
    anchored_fragments = sum(count_value(row, "anchored_fragments") for row in eligible_rows)
    exact_branchpoint_fragments = sum(count_value(row, "exact_branchpoint_fragments") for row in eligible_rows)
    plus_one_branchpoint_fragments = sum(count_value(row, "plus_one_branchpoint_fragments") for row in eligible_rows)
    zero_or_plus_one_branchpoint_fragments = exact_branchpoint_fragments + plus_one_branchpoint_fragments
    profile_window_fragments = sum(total_offset_counts.values())

    summary_row["shared_min_reads_all_samples"] = shared_min_reads
    summary_row["shared_introns"] = len(shared_introns)
    summary_row["raw_anchored_fragments"] = raw_anchored_fragments
    summary_row["raw_anchored_introns_with_reads"] = raw_anchored_introns
    summary_row["anchored_fragments"] = anchored_fragments
    summary_row["anchored_fragments_cpm"] = (
        0.0 if library_fragments == 0 else anchored_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["anchored_introns_with_reads"] = len(eligible_rows)
    summary_row["exact_branchpoint_fragments"] = exact_branchpoint_fragments
    summary_row["exact_branchpoint_cpm"] = (
        0.0 if library_fragments == 0 else exact_branchpoint_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["exact_branchpoint_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else exact_branchpoint_fragments / anchored_fragments
    )
    summary_row["exact_branchpoint_percent_anchored"] = summary_row["exact_branchpoint_fraction_anchored"] * 100.0
    summary_row["plus_one_branchpoint_fragments"] = plus_one_branchpoint_fragments
    summary_row["plus_one_branchpoint_cpm"] = (
        0.0 if library_fragments == 0 else plus_one_branchpoint_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["plus_one_branchpoint_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else plus_one_branchpoint_fragments / anchored_fragments
    )
    summary_row["plus_one_branchpoint_percent_anchored"] = (
        summary_row["plus_one_branchpoint_fraction_anchored"] * 100.0
    )
    summary_row["zero_or_plus_one_branchpoint_fragments"] = zero_or_plus_one_branchpoint_fragments
    summary_row["zero_or_plus_one_branchpoint_cpm"] = (
        0.0 if library_fragments == 0 else zero_or_plus_one_branchpoint_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["zero_or_plus_one_branchpoint_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else zero_or_plus_one_branchpoint_fragments / anchored_fragments
    )
    summary_row["zero_or_plus_one_branchpoint_percent_anchored"] = (
        summary_row["zero_or_plus_one_branchpoint_fraction_anchored"] * 100.0
    )
    summary_row["profile_window_fragments"] = profile_window_fragments
    summary_row["profile_window_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else profile_window_fragments / anchored_fragments
    )
    summary_row["profile_window_percent_anchored"] = summary_row["profile_window_fraction_anchored"] * 100.0
    return summary_row


def plot_results(
    sample_rows,
    condition_rows,
    summary_rows,
    condition_order,
    shared_min_reads,
    output_png,
    output_pdf,
):
    profile_by_sample = defaultdict(list)
    for row in sample_rows:
        profile_by_sample[row["sample"]].append(row)

    profile_by_condition = defaultdict(list)
    for row in condition_rows:
        profile_by_condition[row["condition"]].append(row)

    figure, (ax_profile, ax_exact) = plt.subplots(
        1,
        2,
        figsize=(12, 4.5),
        gridspec_kw={"width_ratios": [3.4, 1.2]},
        constrained_layout=True,
    )

    for sample in sorted(profile_by_sample):
        ordered_rows = sorted(profile_by_sample[sample], key=lambda row: int(row["offset_nt"]))
        condition = ordered_rows[0]["condition"]
        color = CONDITION_COLORS.get(condition, "#4c4c4c")
        ax_profile.plot(
            [int(row["offset_nt"]) for row in ordered_rows],
            [float(row["anchored_percent"]) for row in ordered_rows],
            color=color,
            alpha=0.25,
            linewidth=1.0,
        )

    for condition in condition_order:
        ordered_rows = sorted(profile_by_condition[condition], key=lambda row: int(row["offset_nt"]))
        color = CONDITION_COLORS.get(condition, "#4c4c4c")
        ax_profile.plot(
            [int(row["offset_nt"]) for row in ordered_rows],
            [float(row["mean_anchored_percent"]) for row in ordered_rows],
            color=color,
            linewidth=2.5,
            label=condition,
        )

    ax_profile.axvline(0, color="#4c4c4c", linestyle="--", linewidth=1)
    ax_profile.set_xlim(-60, 10)
    ax_profile.set_xlabel("Read1 5' end offset from selected branchpoint (nt; + toward intron 3' end)")
    ax_profile.set_ylabel("Anchored shared-intron fragments (%)")
    title = "Branchpoint-centred 5' end metaprofile"
    if shared_min_reads > 0:
        title += f"\nShared introns with >= {shared_min_reads} anchored reads in every sample"
    ax_profile.set_title(title)
    ax_profile.legend(frameon=False)

    summary_by_condition = defaultdict(list)
    for row in summary_rows:
        summary_by_condition[row["condition"]].append(row)

    for idx, condition in enumerate(condition_order):
        entries = summary_by_condition[condition]
        values = [float(entry["zero_or_plus_one_branchpoint_percent_anchored"]) for entry in entries]
        color = CONDITION_COLORS.get(condition, "#4c4c4c")
        if values:
            if len(values) == 1:
                jitter = [idx]
            else:
                step = 0.24 / (len(values) - 1)
                jitter = [idx - 0.12 + step * i for i in range(len(values))]
            ax_exact.scatter(jitter, values, color=color, s=36, zorder=3)
            mean_value = float_mean(values)
            sd_value = float_sd(values)
            ax_exact.hlines(mean_value, idx - 0.18, idx + 0.18, color=color, linewidth=2.5)
            if sd_value > 0:
                ax_exact.vlines(idx, mean_value - sd_value, mean_value + sd_value, color=color, linewidth=1.5)

    ax_exact.set_xticks(range(len(condition_order)))
    ax_exact.set_xticklabels(condition_order)
    ax_exact.set_ylabel("Branchpoint-proximal reads (%)")
    ax_exact.set_title("Offsets 0 / +1 fragments")

    for axis in (ax_profile, ax_exact):
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def main():
    args = parse_args()

    raw_metaprofile_rows = []
    for path in args.metaprofiles:
        raw_metaprofile_rows.extend(read_tsv_rows(path))
    offset_range = sorted({int(row["offset_nt"]) for row in raw_metaprofile_rows})

    raw_summary_rows = []
    condition_order = []
    for path in args.summaries:
        rows = read_tsv_rows(path)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one summary row in {path}")
        row = rows[0]
        raw_summary_rows.append(row)
        if row["condition"] not in condition_order:
            condition_order.append(row["condition"])
    raw_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))
    sample_order = [row["sample"] for row in raw_summary_rows]
    raw_summary_by_sample = {row["sample"]: row for row in raw_summary_rows}

    site_counts_by_sample = {}
    site_metadata = {}
    for path in args.site_counts:
        rows = read_tsv_rows(path)
        sample = infer_sample_name(rows, path)
        sample_rows = {}
        for row in rows:
            sample_rows[row["intron_id"]] = row
            site_metadata.setdefault(row["intron_id"], row)
        site_counts_by_sample[sample] = sample_rows

    intron_offsets_by_sample = {}
    for path in args.intron_offsets:
        rows = read_tsv_rows(path)
        sample = infer_sample_name(rows, path)
        sample_offsets = defaultdict(Counter)
        for row in rows:
            sample_offsets[row["intron_id"]][int(row["offset_nt"])] = int(row["read_count"])
        intron_offsets_by_sample[sample] = sample_offsets

    shared_introns = build_shared_intron_set(site_counts_by_sample, sample_order, args.shared_min_reads)
    shared_intron_rows, shared_intron_fieldnames = build_shared_introns_rows(
        shared_introns,
        site_counts_by_sample,
        site_metadata,
        sample_order,
    )

    sample_profile_rows = []
    sample_summary_rows = []
    for sample in sample_order:
        raw_summary_row = raw_summary_by_sample[sample]
        total_offset_counts = aggregate_offset_counts(intron_offsets_by_sample.get(sample, {}), shared_introns)
        sample_profile_rows.extend(
            build_sample_metaprofile_rows(
                sample,
                raw_summary_row["condition"],
                count_value(raw_summary_row, "library_fragments"),
                sum(count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "anchored_fragments") for intron_id in shared_introns),
                offset_range,
                total_offset_counts,
            )
        )
        sample_summary_rows.append(
            build_sample_summary_row(
                raw_summary_row,
                site_counts_by_sample.get(sample, {}),
                total_offset_counts,
                shared_introns,
                args.shared_min_reads,
            )
        )

    sample_profile_rows.sort(
        key=lambda row: (condition_order.index(row["condition"]), row["sample"], int(row["offset_nt"]))
    )
    sample_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))

    condition_rows = summarise_condition_profiles(sample_profile_rows, condition_order)
    condition_summary_rows = summarise_condition_rows(sample_summary_rows, condition_order)

    write_rows(args.output_metaprofile_by_sample, sample_profile_rows, list(sample_profile_rows[0].keys()))
    write_rows(args.output_metaprofile_by_condition, condition_rows, list(condition_rows[0].keys()))
    write_rows(args.output_summary_by_sample, sample_summary_rows, list(sample_summary_rows[0].keys()))
    write_rows(args.output_summary_by_condition, condition_summary_rows, list(condition_summary_rows[0].keys()))
    write_rows(args.output_shared_introns, shared_intron_rows, shared_intron_fieldnames)

    plot_results(
        sample_profile_rows,
        condition_rows,
        sample_summary_rows,
        condition_order,
        args.shared_min_reads,
        args.output_plot_png,
        args.output_plot_pdf,
    )

    print(f"Shared introns retained: {len(shared_introns)}")
    print(f"Minimum anchored reads in all samples: {args.shared_min_reads}")
    print(f"Metaprofile rows aggregated: {len(sample_profile_rows)}")
    print(f"Summary rows aggregated: {len(sample_summary_rows)}")


if __name__ == "__main__":
    main()
