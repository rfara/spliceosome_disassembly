#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
import numpy as np
import pysam

from splice_site_maxent import score3, score5

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DNA_BASES = ["A", "C", "G", "T"]
RNA_BASE_LABELS = ["A", "C", "G", "U"]
QUERY_COLOR = "#d95f02"
CONTROL_COLOR = "#1f77b4"
TOP_COLOR = "#b2182b"
BOTTOM_COLOR = "#2166ac"
MODEL_ORDER = {
    "null_intercept_only": 0,
    "baseline_control_only": 1,
    "baseline_plus_numeric": 2,
    "baseline_plus_bp_flank": 3,
    "baseline_plus_numeric_plus_bp_flank": 4,
}
FEATURE_LABELS = {
    "log_intron_length": "Intron length",
    "branchpoint_to_3ss_nt": "BP to 3'SS distance",
    "branchpoint_score": "BP score",
    "donor_maxent": "5'SS MaxEnt",
    "acceptor_maxent": "3'SS MaxEnt",
    "anchored_enrichment_log2_fold_change": "Anchored abundance log2(DIS/ILS)",
    "intron_gc": "Intron GC",
    "intron_gc_5p_window": "5' intron GC (100 nt)",
    "intron_gc_3p_window": "3' intron GC (100 nt)",
    "bp_to_3ss_pyrimidine_fraction": "BP to 3'SS pyrimidine fraction",
    "three_prime_window_pyrimidine_fraction": "3' window pyrimidine fraction",
    "branchpoint_candidates": "BP candidate count",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--signal-field", default="zero_or_plus_one_branchpoint_fragments")
    parser.add_argument("--control-condition", default="ILS")
    parser.add_argument("--query-condition", default="DIS")
    parser.add_argument("--control-min-reads", type=int, default=20)
    parser.add_argument("--query-min-reads", type=int, default=20)
    parser.add_argument("--require-bp-center-base", default="A")
    parser.add_argument("--bp-upstream", type=int, default=5)
    parser.add_argument("--bp-downstream", type=int, default=5)
    parser.add_argument("--intron-gc-window", type=int, default=100)
    parser.add_argument("--three-prime-window", type=int, default=20)
    parser.add_argument("--residual-quantile-fraction", type=float, default=0.1)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-shared-introns", required=True)
    parser.add_argument("--output-model-summary", required=True)
    parser.add_argument("--output-feature-summary", required=True)
    parser.add_argument("--output-bp-sequence-effects", required=True)
    parser.add_argument("--output-bp-sequence-enrichment", required=True)
    parser.add_argument("--output-plot-png", required=True)
    parser.add_argument("--output-plot-pdf", required=True)
    parser.add_argument("--output-group-comparison", required=True)
    parser.add_argument("--output-feature-distribution-plot-png", required=True)
    parser.add_argument("--output-feature-distribution-plot-pdf", required=True)
    return parser.parse_args()


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_tsv_rows(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def count_value(row, field):
    raw_value = row.get(field, 0)
    if raw_value in {"", None}:
        return 0
    return int(float(raw_value))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def signal_label(signal_field):
    labels = {
        "exact_branchpoint_fragments": "exact branchpoint-terminated fragments",
        "plus_one_branchpoint_fragments": "+1 branchpoint-proximal fragments",
        "zero_or_plus_one_branchpoint_fragments": "0/+1 branchpoint-proximal fragments",
    }
    return labels.get(signal_field, signal_field.replace("_", " "))


def load_condition_rows(paths, signal_field):
    condition_rows = defaultdict(dict)
    for path in paths:
        rows = read_tsv_rows(path)
        if not rows:
            continue
        condition = rows[0]["condition"]
        for row in rows:
            intron_id = row["intron_id"]
            pooled = condition_rows[condition].setdefault(
                intron_id,
                {
                    "condition": condition,
                    "intron_id": intron_id,
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                    "transcript_id": row["transcript_id"],
                    "intron_number": row["intron_number"],
                    "chrom": row["chrom"],
                    "strand": row["strand"],
                    "intron_start": int(row["intron_start"]),
                    "intron_end": int(row["intron_end"]),
                    "intron_length": int(row["intron_length"]),
                    "three_prime_ss": int(row["three_prime_ss"]),
                    "branchpoint_position": int(row["branchpoint_position"]),
                    "branchpoint_score": float(row["branchpoint_score"]),
                    "branchpoint_to_3ss_nt": int(row["branchpoint_to_3ss_nt"]),
                    "branchpoint_type": row.get("branchpoint_type", ""),
                    "branchpoint_candidates": int(row["branchpoint_candidates"]),
                    "anchored_fragments": 0,
                    "signal_fragments": 0,
                },
            )
            pooled["anchored_fragments"] += count_value(row, "anchored_fragments")
            pooled["signal_fragments"] += count_value(row, signal_field)
    return condition_rows


def fetch_raw_sequence(fasta, chrom, start_0based, end_0based):
    if start_0based < 0 or end_0based <= start_0based:
        return None
    try:
        ref_length = fasta.get_reference_length(chrom)
    except KeyError:
        return None
    if end_0based > ref_length:
        return None
    return fasta.fetch(chrom, start_0based, end_0based).upper()


def fetch_oriented_sequence(fasta, chrom, start_1based, end_1based, strand):
    sequence = fetch_raw_sequence(fasta, chrom, start_1based - 1, end_1based)
    if sequence is None:
        return None
    if strand == "-":
        return reverse_complement(sequence)
    return sequence


def donor_window_sequence(record, fasta):
    if record["strand"] == "+":
        sequence = fetch_raw_sequence(
            fasta,
            record["chrom"],
            record["intron_start"] - 1 - 3,
            record["intron_start"] - 1 + 6,
        )
    else:
        sequence = fetch_raw_sequence(
            fasta,
            record["chrom"],
            record["intron_end"] - 6,
            record["intron_end"] + 3,
        )
        if sequence is not None:
            sequence = reverse_complement(sequence)
    if sequence is None or len(sequence) != 9:
        return None
    return sequence


def acceptor_window_sequence(record, fasta):
    if record["strand"] == "+":
        sequence = fetch_raw_sequence(
            fasta,
            record["chrom"],
            record["intron_end"] - 20,
            record["intron_end"] + 3,
        )
    else:
        sequence = fetch_raw_sequence(
            fasta,
            record["chrom"],
            record["intron_start"] - 1 - 3,
            record["intron_start"] - 1 + 20,
        )
        if sequence is not None:
            sequence = reverse_complement(sequence)
    if sequence is None or len(sequence) != 23:
        return None
    return sequence


def gc_fraction(sequence):
    if not sequence:
        return math.nan
    valid = sum(1 for base in sequence if base in {"A", "C", "G", "T"})
    if valid == 0:
        return math.nan
    gc = sum(1 for base in sequence if base in {"C", "G"})
    return gc / valid


def pyrimidine_fraction(sequence):
    if not sequence:
        return math.nan
    valid = sum(1 for base in sequence if base in {"A", "C", "G", "T"})
    if valid == 0:
        return math.nan
    pyrimidines = sum(1 for base in sequence if base in {"C", "T"})
    return pyrimidines / valid


def empirical_probability(successes, totals):
    return (successes + 0.5) / (totals + 1.0)


def logit(probability):
    probability = min(max(probability, 1e-9), 1.0 - 1e-9)
    return math.log(probability / (1.0 - probability))


def sigmoid(values):
    positive = values >= 0
    result = np.empty_like(values, dtype=float)
    result[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exp_values = np.exp(values[~positive])
    result[~positive] = exp_values / (1.0 + exp_values)
    return result


def binomial_deviance(successes, totals, probabilities):
    probabilities = np.clip(probabilities, 1e-9, 1.0 - 1e-9)
    fitted_successes = totals * probabilities
    fitted_failures = totals - fitted_successes
    failures = totals - successes
    deviance = np.zeros_like(probabilities, dtype=float)
    nonzero_successes = successes > 0
    nonzero_failures = failures > 0
    deviance[nonzero_successes] += successes[nonzero_successes] * np.log(
        successes[nonzero_successes] / fitted_successes[nonzero_successes]
    )
    deviance[nonzero_failures] += failures[nonzero_failures] * np.log(
        failures[nonzero_failures] / fitted_failures[nonzero_failures]
    )
    return float(2.0 * np.sum(deviance))


def fit_binomial_glm(design_matrix, successes, totals, ridge=1e-6, max_iter=100, tolerance=1e-8):
    beta = np.zeros(design_matrix.shape[1], dtype=float)
    initial_probability = empirical_probability(float(successes.sum()), float(totals.sum()))
    beta[0] = logit(initial_probability)

    penalty = np.zeros(design_matrix.shape[1], dtype=float)
    penalty[1:] = ridge
    covariance = None

    for _ in range(max_iter):
        eta = design_matrix @ beta
        probabilities = np.clip(sigmoid(eta), 1e-9, 1.0 - 1e-9)
        weights = totals * probabilities * (1.0 - probabilities)
        weights = np.maximum(weights, 1e-9)
        working_response = eta + (successes - totals * probabilities) / weights

        weighted_design = design_matrix * np.sqrt(weights)[:, None]
        weighted_response = working_response * np.sqrt(weights)
        lhs = weighted_design.T @ weighted_design
        lhs.flat[:: lhs.shape[0] + 1] += penalty
        rhs = weighted_design.T @ weighted_response

        try:
            new_beta = np.linalg.solve(lhs, rhs)
            covariance = np.linalg.inv(lhs)
        except np.linalg.LinAlgError:
            new_beta = np.linalg.lstsq(lhs, rhs, rcond=None)[0]
            covariance = np.linalg.pinv(lhs)

        if np.max(np.abs(new_beta - beta)) < tolerance:
            beta = new_beta
            break
        beta = new_beta

    probabilities = np.clip(sigmoid(design_matrix @ beta), 1e-9, 1.0 - 1e-9)
    if covariance is None:
        weights = totals * probabilities * (1.0 - probabilities)
        weighted_design = design_matrix * np.sqrt(np.maximum(weights, 1e-9))[:, None]
        lhs = weighted_design.T @ weighted_design
        lhs.flat[:: lhs.shape[0] + 1] += penalty
        covariance = np.linalg.pinv(lhs)
    return beta, probabilities, covariance


def weighted_r_squared(observed_fraction, predicted_fraction, weights):
    weighted_mean = float(np.average(observed_fraction, weights=weights))
    residual_sum = float(np.sum(weights * (observed_fraction - predicted_fraction) ** 2))
    total_sum = float(np.sum(weights * (observed_fraction - weighted_mean) ** 2))
    if total_sum == 0.0:
        return 0.0
    return 1.0 - (residual_sum / total_sum)


def build_shared_rows(
    condition_records,
    fasta_path,
    control_condition,
    query_condition,
    control_min_reads,
    query_min_reads,
    require_bp_center_base,
    bp_upstream,
    bp_downstream,
    intron_gc_window,
    three_prime_window,
):
    control_records = condition_records[control_condition]
    query_records = condition_records[query_condition]
    total_control_anchored = sum(record["anchored_fragments"] for record in control_records.values())
    total_query_anchored = sum(record["anchored_fragments"] for record in query_records.values())
    if total_control_anchored <= 0 or total_query_anchored <= 0:
        raise ValueError("Condition-level anchored fragment totals must be positive")
    shared_introns = sorted(set(control_records) & set(query_records))
    rows = []

    with pysam.FastaFile(fasta_path) as fasta:
        for intron_id in shared_introns:
            control = control_records[intron_id]
            query = query_records[intron_id]
            if control["anchored_fragments"] < control_min_reads or query["anchored_fragments"] < query_min_reads:
                continue

            intron_sequence = fetch_oriented_sequence(
                fasta,
                control["chrom"],
                control["intron_start"],
                control["intron_end"],
                control["strand"],
            )
            branchpoint_context = fetch_oriented_sequence(
                fasta,
                control["chrom"],
                control["branchpoint_position"] - bp_upstream,
                control["branchpoint_position"] + bp_downstream,
                control["strand"],
            )
            donor_sequence = donor_window_sequence(control, fasta)
            acceptor_sequence = acceptor_window_sequence(control, fasta)

            if intron_sequence is None or branchpoint_context is None:
                continue
            if len(branchpoint_context) != bp_upstream + bp_downstream + 1:
                continue
            if any(base not in DNA_BASES for base in branchpoint_context):
                continue

            branchpoint_center_base = branchpoint_context[bp_upstream]
            if require_bp_center_base and branchpoint_center_base != require_bp_center_base:
                continue

            if control["strand"] == "+":
                branchpoint_index = control["branchpoint_position"] - control["intron_start"]
            else:
                branchpoint_index = control["intron_end"] - control["branchpoint_position"]
            if branchpoint_index < 0 or branchpoint_index >= len(intron_sequence):
                continue

            bp_to_3ss_sequence = intron_sequence[branchpoint_index + 1 : -2]
            three_prime_window_sequence = intron_sequence[-min(three_prime_window, len(intron_sequence)) :]
            five_prime_gc_sequence = intron_sequence[: min(intron_gc_window, len(intron_sequence))]
            three_prime_gc_sequence = intron_sequence[-min(intron_gc_window, len(intron_sequence)) :]

            control_fraction = 0.0 if control["anchored_fragments"] == 0 else control["signal_fragments"] / control["anchored_fragments"]
            query_fraction = 0.0 if query["anchored_fragments"] == 0 else query["signal_fragments"] / query["anchored_fragments"]
            control_smoothed = empirical_probability(control["signal_fragments"], control["anchored_fragments"])
            query_smoothed = empirical_probability(query["signal_fragments"], query["anchored_fragments"])
            control_global_anchored_fraction = control["anchored_fragments"] / total_control_anchored
            query_global_anchored_fraction = query["anchored_fragments"] / total_query_anchored

            row = {
                "intron_id": intron_id,
                "gene_id": control["gene_id"],
                "gene_name": control["gene_name"],
                "transcript_id": control["transcript_id"],
                "intron_number": control["intron_number"],
                "chrom": control["chrom"],
                "strand": control["strand"],
                "intron_start": control["intron_start"],
                "intron_end": control["intron_end"],
                "intron_length": control["intron_length"],
                "branchpoint_position": control["branchpoint_position"],
                "branchpoint_to_3ss_nt": control["branchpoint_to_3ss_nt"],
                "branchpoint_score": control["branchpoint_score"],
                "branchpoint_type": control["branchpoint_type"],
                "branchpoint_candidates": control["branchpoint_candidates"],
                "branchpoint_context": branchpoint_context,
                "branchpoint_center_base": branchpoint_center_base,
                "donor_sequence": donor_sequence or "",
                "acceptor_sequence": acceptor_sequence or "",
                "control_anchored_fragments": control["anchored_fragments"],
                "control_branched_fragments": control["signal_fragments"],
                "control_branch_fraction": control_fraction,
                "control_smoothed_branch_fraction": control_smoothed,
                "control_empirical_logit": logit(control_smoothed),
                "query_anchored_fragments": query["anchored_fragments"],
                "query_branched_fragments": query["signal_fragments"],
                "query_branch_fraction": query_fraction,
                "query_smoothed_branch_fraction": query_smoothed,
                "query_empirical_logit": logit(query_smoothed),
                "delta_branch_fraction": query_fraction - control_fraction,
                "control_global_anchored_fraction": control_global_anchored_fraction,
                "query_global_anchored_fraction": query_global_anchored_fraction,
                "anchored_enrichment_log2_fold_change": math.log2(
                    query_global_anchored_fraction / control_global_anchored_fraction
                ),
                "log_intron_length": math.log10(max(control["intron_length"], 1)),
                "donor_maxent": score5(donor_sequence),
                "acceptor_maxent": score3(acceptor_sequence),
                "intron_gc": gc_fraction(intron_sequence),
                "intron_gc_5p_window": gc_fraction(five_prime_gc_sequence),
                "intron_gc_3p_window": gc_fraction(three_prime_gc_sequence),
                "bp_to_3ss_pyrimidine_fraction": pyrimidine_fraction(bp_to_3ss_sequence),
                "three_prime_window_pyrimidine_fraction": pyrimidine_fraction(three_prime_window_sequence),
            }
            rows.append(row)
    return rows


def finite_feature_rows(rows, feature_names):
    filtered = []
    for row in rows:
        if any(not math.isfinite(float(row[feature])) for feature in feature_names):
            continue
        filtered.append(row)
    return filtered


def build_design_matrix(rows, numeric_features=(), bp_positions=()):
    weights = np.array([row["query_anchored_fragments"] for row in rows], dtype=float)
    columns = [
        np.ones(len(rows), dtype=float),
        np.array([row["control_empirical_logit"] for row in rows], dtype=float),
    ]
    metadata = [
        {"kind": "intercept"},
        {"kind": "control_logit"},
    ]

    for feature in numeric_features:
        values = np.array([row[feature] for row in rows], dtype=float)
        mean = float(np.average(values, weights=weights))
        variance = float(np.average((values - mean) ** 2, weights=weights))
        std = math.sqrt(max(variance, 1e-12))
        standardized = (values - mean) / std
        columns.append(standardized)
        metadata.append(
            {
                "kind": "numeric",
                "feature": feature,
                "mean": mean,
                "std": std,
            }
        )

    for index, position in enumerate(bp_positions):
        bases = [row["branchpoint_context"][index] for row in rows]
        counts = Counter(bases)
        if len(counts) < 2:
            continue
        baseline_base = sorted(counts, key=lambda base: (-counts[base], base))[0]
        for base in sorted(counts):
            if base == baseline_base:
                continue
            column = np.array([1.0 if observed == base else 0.0 for observed in bases], dtype=float)
            if column.sum() == 0:
                continue
            columns.append(column)
            metadata.append(
                {
                    "kind": "bp_flank",
                    "position": position,
                    "base": base,
                    "baseline_base": baseline_base,
                }
            )

    return np.column_stack(columns), metadata


def fit_model(rows, numeric_features=(), bp_positions=(), label=""):
    design_matrix, metadata = build_design_matrix(rows, numeric_features=numeric_features, bp_positions=bp_positions)
    successes = np.array([row["query_branched_fragments"] for row in rows], dtype=float)
    totals = np.array([row["query_anchored_fragments"] for row in rows], dtype=float)
    observed_fraction = successes / totals
    beta, probabilities, covariance = fit_binomial_glm(design_matrix, successes, totals)
    return {
        "label": label,
        "design_matrix": design_matrix,
        "metadata": metadata,
        "beta": beta,
        "covariance": covariance,
        "predicted_fraction": probabilities,
        "residual_deviance": binomial_deviance(successes, totals, probabilities),
        "weighted_r_squared": weighted_r_squared(observed_fraction, probabilities, totals),
    }


def normal_pvalue_from_z(z_score):
    return math.erfc(abs(z_score) / math.sqrt(2.0))


def overall_fraction(rows, count_field, total_field):
    numerator = sum(row[count_field] for row in rows)
    denominator = sum(row[total_field] for row in rows)
    return 0.0 if denominator == 0 else numerator / denominator


def build_model_summary(
    rows,
    control_condition,
    query_condition,
    model_lookup,
    null_deviance,
    baseline_deviance,
):
    total_query_anchored = sum(row["query_anchored_fragments"] for row in rows)
    total_query_branched = sum(row["query_branched_fragments"] for row in rows)
    summary_rows = []
    for name, model in model_lookup.items():
        residual_deviance = model["residual_deviance"]
        total_explained = 0.0 if null_deviance == 0.0 else 1.0 - (residual_deviance / null_deviance)
        added_over_baseline = max(0.0, baseline_deviance - residual_deviance)
        added_fraction_of_baseline = 0.0 if baseline_deviance == 0.0 else added_over_baseline / baseline_deviance
        summary_rows.append(
            {
                "model": name,
                "control_condition": control_condition,
                "query_condition": query_condition,
                "introns": len(rows),
                "query_anchored_fragments": total_query_anchored,
                "query_branched_fragments": total_query_branched,
                "feature_columns": max(0, model["design_matrix"].shape[1] - 2),
                "residual_deviance": residual_deviance,
                "total_deviance_explained": total_explained,
                "weighted_r_squared": model["weighted_r_squared"],
                "added_deviance_over_baseline": added_over_baseline,
                "added_fraction_of_baseline_residual": added_fraction_of_baseline,
            }
        )
    return summary_rows


def build_numeric_feature_summary(rows, feature_names, null_deviance, baseline_model):
    baseline_deviance = baseline_model["residual_deviance"]
    baseline_r2 = baseline_model["weighted_r_squared"]
    summary_rows = []
    for feature in feature_names:
        model = fit_model(rows, numeric_features=(feature,), label=f"baseline_plus_{feature}")
        coefficient_index = 2
        coefficient = float(model["beta"][coefficient_index])
        standard_error = float(math.sqrt(max(model["covariance"][coefficient_index, coefficient_index], 0.0)))
        z_score = 0.0 if standard_error == 0.0 else coefficient / standard_error
        residual_deviance = model["residual_deviance"]
        summary_rows.append(
            {
                "feature": feature,
                "coefficient_log_odds_per_sd": coefficient,
                "odds_ratio_per_sd": math.exp(coefficient),
                "standard_error": standard_error,
                "z_score": z_score,
                "pvalue_two_sided": normal_pvalue_from_z(z_score),
                "residual_deviance": residual_deviance,
                "total_deviance_explained": 0.0 if null_deviance == 0.0 else 1.0 - (residual_deviance / null_deviance),
                "added_deviance_over_baseline": max(0.0, baseline_deviance - residual_deviance),
                "added_fraction_of_baseline_residual": 0.0
                if baseline_deviance == 0.0
                else max(0.0, baseline_deviance - residual_deviance) / baseline_deviance,
                "weighted_r_squared": model["weighted_r_squared"],
                "delta_weighted_r_squared_over_baseline": model["weighted_r_squared"] - baseline_r2,
            }
        )
    summary_rows.sort(key=lambda row: (-row["added_fraction_of_baseline_residual"], row["feature"]))
    return summary_rows


def add_baseline_predictions(rows, baseline_model):
    predicted = baseline_model["predicted_fraction"]
    for row, predicted_fraction in zip(rows, predicted, strict=True):
        row["baseline_expected_query_branch_fraction"] = float(predicted_fraction)
        row["baseline_residual_query_branch_fraction"] = row["query_branch_fraction"] - float(predicted_fraction)
        predicted_logit = logit(float(predicted_fraction))
        row["baseline_expected_query_logit"] = predicted_logit
        row["baseline_residual_query_logit"] = row["query_empirical_logit"] - predicted_logit


def build_bp_sequence_effect_rows(rows, bp_positions):
    overall_residual = np.average(
        np.array([row["baseline_residual_query_branch_fraction"] for row in rows], dtype=float),
        weights=np.array([row["query_anchored_fragments"] for row in rows], dtype=float),
    )
    effect_rows = []
    for index, position in enumerate(bp_positions):
        for base in DNA_BASES:
            matching = [row for row in rows if row["branchpoint_context"][index] == base]
            if not matching:
                continue
            query_weights = np.array([row["query_anchored_fragments"] for row in matching], dtype=float)
            control_weights = np.array([row["control_anchored_fragments"] for row in matching], dtype=float)
            residuals = np.array([row["baseline_residual_query_branch_fraction"] for row in matching], dtype=float)
            observed = np.array([row["query_branch_fraction"] for row in matching], dtype=float)
            expected = np.array([row["baseline_expected_query_branch_fraction"] for row in matching], dtype=float)
            effect_rows.append(
                {
                    "position": position,
                    "base": base,
                    "intron_count": len(matching),
                    "query_anchored_fragments": int(query_weights.sum()),
                    "control_anchored_fragments": int(control_weights.sum()),
                    "mean_query_branch_fraction": float(np.average(observed, weights=query_weights)),
                    "mean_expected_query_branch_fraction": float(np.average(expected, weights=query_weights)),
                    "mean_residual_query_branch_fraction": float(np.average(residuals, weights=query_weights)),
                    "delta_residual_from_dataset_mean": float(np.average(residuals, weights=query_weights) - overall_residual),
                }
            )
    return effect_rows


def build_bp_sequence_enrichment_rows(rows, bp_positions, quantile_fraction):
    ordered = sorted(
        rows,
        key=lambda row: (
            -row["baseline_residual_query_branch_fraction"],
            -row["query_anchored_fragments"],
            row["gene_name"],
            row["intron_id"],
        ),
    )
    group_size = max(1, int(round(len(ordered) * quantile_fraction)))
    top_rows = ordered[:group_size]
    bottom_rows = ordered[-group_size:]

    enrichment_rows = []
    for index, position in enumerate(bp_positions):
        for base in DNA_BASES:
            top_count = sum(1 for row in top_rows if row["branchpoint_context"][index] == base)
            bottom_count = sum(1 for row in bottom_rows if row["branchpoint_context"][index] == base)
            top_other = len(top_rows) - top_count
            bottom_other = len(bottom_rows) - bottom_count
            odds_ratio = ((top_count + 0.5) / (top_other + 0.5)) / ((bottom_count + 0.5) / (bottom_other + 0.5))
            enrichment_rows.append(
                {
                    "position": position,
                    "base": base,
                    "top_group_introns": len(top_rows),
                    "bottom_group_introns": len(bottom_rows),
                    "top_group_base_count": top_count,
                    "bottom_group_base_count": bottom_count,
                    "top_group_base_fraction": top_count / len(top_rows),
                    "bottom_group_base_fraction": bottom_count / len(bottom_rows),
                    "log2_odds_ratio_top_vs_bottom": math.log2(odds_ratio),
                }
            )
    return enrichment_rows, top_rows, bottom_rows


def distribution_specs(feature_summary_rows):
    specs = []
    for row in feature_summary_rows:
        model_feature = row["feature"]
        if model_feature == "branchpoint_candidates":
            continue
        plot_feature = "intron_length" if model_feature == "log_intron_length" else model_feature
        specs.append(
            {
                "model_feature": model_feature,
                "plot_feature": plot_feature,
                "label": FEATURE_LABELS.get(model_feature, model_feature.replace("_", " ")),
                "log_scale": plot_feature == "intron_length",
            }
        )
    return specs


def build_group_comparison_rows(rows, top_rows, specs):
    top_ids = {row["intron_id"] for row in top_rows}
    other_rows = [row for row in rows if row["intron_id"] not in top_ids]
    comparison_rows = []
    for spec in specs:
        plot_feature = spec["plot_feature"]
        top_values = np.array([float(row[plot_feature]) for row in top_rows], dtype=float)
        other_values = np.array([float(row[plot_feature]) for row in other_rows], dtype=float)
        top_mean = float(np.mean(top_values))
        other_mean = float(np.mean(other_values))
        top_median = float(np.median(top_values))
        other_median = float(np.median(other_values))
        top_var = float(np.var(top_values, ddof=1)) if len(top_values) > 1 else 0.0
        other_var = float(np.var(other_values, ddof=1)) if len(other_values) > 1 else 0.0
        pooled_sd = math.sqrt(max((top_var + other_var) / 2.0, 1e-12))
        comparison_rows.append(
            {
                "model_feature": spec["model_feature"],
                "plot_feature": plot_feature,
                "label": spec["label"],
                "top_group_introns": len(top_rows),
                "other_introns": len(other_rows),
                "top_mean": top_mean,
                "other_mean": other_mean,
                "delta_mean": top_mean - other_mean,
                "top_median": top_median,
                "other_median": other_median,
                "delta_median": top_median - other_median,
                "top_q25": float(np.quantile(top_values, 0.25)),
                "top_q75": float(np.quantile(top_values, 0.75)),
                "other_q25": float(np.quantile(other_values, 0.25)),
                "other_q75": float(np.quantile(other_values, 0.75)),
                "standardized_mean_difference": (top_mean - other_mean) / pooled_sd,
            }
        )
    return comparison_rows


def rows_to_matrix(rows, value_field, positions):
    position_to_index = {position: index for index, position in enumerate(positions)}
    base_to_index = {base: index for index, base in enumerate(DNA_BASES)}
    matrix = np.full((len(DNA_BASES), len(positions)), np.nan, dtype=float)
    for row in rows:
        matrix[base_to_index[row["base"]], position_to_index[row["position"]]] = row[value_field]
    return matrix


def draw_heatmap(axis, matrix, positions, title, colorbar_label, cmap, value_limit):
    image = axis.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=-value_limit,
        vmax=value_limit,
    )
    axis.set_xticks(range(len(positions)))
    axis.set_xticklabels(positions)
    axis.set_yticks(range(len(RNA_BASE_LABELS)))
    axis.set_yticklabels(RNA_BASE_LABELS)
    axis.set_xlabel("BP flank position")
    axis.set_title(title)
    return image, colorbar_label


def plot_results(
    rows,
    control_condition,
    query_condition,
    signal_field,
    bp_positions,
    model_summary_rows,
    feature_summary_rows,
    bp_effect_rows,
    bp_enrichment_rows,
    top_rows,
    bottom_rows,
    output_png,
    output_pdf,
):
    figure, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    baseline_model_row = next(row for row in model_summary_rows if row["model"] == "baseline_control_only")
    full_model_row = next(row for row in model_summary_rows if row["model"] == "baseline_plus_numeric_plus_bp_flank")

    scatter_axis = axes[0, 0]
    control_fraction = np.array([row["control_branch_fraction"] for row in rows], dtype=float)
    query_fraction = np.array([row["query_branch_fraction"] for row in rows], dtype=float)
    scatter_axis.scatter(control_fraction, query_fraction, s=8, color="#bbbbbb", alpha=0.35, linewidths=0)
    scatter_axis.scatter(
        [row["control_branch_fraction"] for row in bottom_rows],
        [row["query_branch_fraction"] for row in bottom_rows],
        s=14,
        color=BOTTOM_COLOR,
        alpha=0.6,
        linewidths=0,
        label="Bottom residual decile",
    )
    scatter_axis.scatter(
        [row["control_branch_fraction"] for row in top_rows],
        [row["query_branch_fraction"] for row in top_rows],
        s=14,
        color=TOP_COLOR,
        alpha=0.6,
        linewidths=0,
        label="Top residual decile",
    )

    x_grid = np.linspace(0.001, 0.999, 400)
    baseline_beta = np.array([baseline_model_row["intercept_beta"], baseline_model_row["control_logit_beta"]], dtype=float)
    predicted_curve = sigmoid(baseline_beta[0] + baseline_beta[1] * np.array([logit(value) for value in x_grid], dtype=float))
    scatter_axis.plot(x_grid, predicted_curve, color="black", linewidth=2.0, label=f"Baseline {query_condition}~{control_condition}")
    scatter_axis.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=1.0, alpha=0.4)
    scatter_axis.set_xlabel(f"{control_condition} branch fraction")
    scatter_axis.set_ylabel(f"{query_condition} branch fraction")
    scatter_axis.set_title(f"{query_condition} vs {control_condition} by intron")
    scatter_axis.legend(frameon=False, loc="upper left")

    bar_axis = axes[0, 1]
    top_features = feature_summary_rows[:8]
    feature_labels = [row["feature"].replace("_", "\n") for row in reversed(top_features)]
    feature_values = [row["added_fraction_of_baseline_residual"] * 100.0 for row in reversed(top_features)]
    bar_axis.barh(feature_labels, feature_values, color=QUERY_COLOR, alpha=0.85)
    bar_axis.set_xlabel("Residual deviance removed beyond baseline (%)")
    bar_axis.set_title("Top numeric intron features")
    bar_axis.axvline(0.0, color="black", linewidth=0.8)
    bar_axis.text(
        0.98,
        0.02,
        "\n".join(
            [
                f"Signal: {signal_label(signal_field)}",
                f"Baseline explained: {baseline_model_row['total_deviance_explained'] * 100.0:.1f}%",
                f"Full model explained: {full_model_row['total_deviance_explained'] * 100.0:.1f}%",
                f"Full added over baseline: {full_model_row['added_fraction_of_baseline_residual'] * 100.0:.1f}%",
            ]
        ),
        transform=bar_axis.transAxes,
        ha="right",
        va="bottom",
    )

    bp_effect_matrix = rows_to_matrix(bp_effect_rows, "mean_residual_query_branch_fraction", bp_positions) * 100.0
    effect_limit = max(
        float(np.nanmax(np.abs(bp_effect_matrix))) if not np.isnan(bp_effect_matrix).all() else 0.0,
        0.5,
    )
    image, label = draw_heatmap(
        axes[1, 0],
        bp_effect_matrix,
        bp_positions,
        f"BP context effect on {query_condition} residual\nObserved minus baseline-expected branch fraction (pp)",
        "Residual branch fraction (pp)",
        "coolwarm",
        effect_limit,
    )
    figure.colorbar(image, ax=axes[1, 0], label=label, shrink=0.86)

    bp_enrichment_matrix = rows_to_matrix(bp_enrichment_rows, "log2_odds_ratio_top_vs_bottom", bp_positions)
    enrichment_limit = max(
        float(np.nanmax(np.abs(bp_enrichment_matrix))) if not np.isnan(bp_enrichment_matrix).all() else 0.0,
        0.5,
    )
    image, label = draw_heatmap(
        axes[1, 1],
        bp_enrichment_matrix,
        bp_positions,
        "Top vs bottom residual introns\nBP context log2 odds ratio",
        "Log2 odds ratio",
        "PiYG",
        enrichment_limit,
    )
    figure.colorbar(image, ax=axes[1, 1], label=label, shrink=0.86)

    for axis in axes.flat:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        f"{query_condition} residual branchiness relative to {control_condition}\n"
        f"{len(rows)} shared introns after filters, baseline explained {baseline_model_row['total_deviance_explained'] * 100.0:.1f}% of query deviance, "
        f"full model {full_model_row['total_deviance_explained'] * 100.0:.1f}%",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def plot_group_feature_distributions(
    rows,
    top_rows,
    specs,
    query_condition,
    control_condition,
    output_png,
    output_pdf,
):
    top_ids = {row["intron_id"] for row in top_rows}
    other_rows = [row for row in rows if row["intron_id"] not in top_ids]
    n_features = len(specs)
    ncols = 2
    nrows = int(math.ceil(n_features / ncols))
    figure, axes = plt.subplots(nrows, ncols, figsize=(12.5, 2.8 * nrows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for axis, spec in zip(axes, specs, strict=False):
        plot_feature = spec["plot_feature"]
        top_values = np.array([float(row[plot_feature]) for row in top_rows], dtype=float)
        other_values = np.array([float(row[plot_feature]) for row in other_rows], dtype=float)
        boxplot = axis.boxplot(
            [other_values, top_values],
            vert=False,
            patch_artist=True,
            showfliers=False,
            widths=0.6,
            tick_labels=["Other", "Top residual decile"],
        )
        for patch, color in zip(boxplot["boxes"], ["#c7c7c7", QUERY_COLOR], strict=True):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        for median in boxplot["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)
        if spec["log_scale"]:
            axis.set_xscale("log")
            axis.set_xlabel("nt")
        elif "fraction" in plot_feature or plot_feature.startswith("intron_gc"):
            axis.set_xlim(0.0, 1.0)
            axis.set_xlabel("fraction")
        elif plot_feature.endswith("_maxent") or plot_feature == "branchpoint_score":
            axis.set_xlabel("score")
        else:
            axis.set_xlabel("value")
        axis.set_title(spec["label"])
        delta_median = float(np.median(top_values) - np.median(other_values))
        if plot_feature == "intron_length":
            median_text = f"{np.median(top_values):.0f} vs {np.median(other_values):.0f} nt"
            delta_text = f"log10 delta {np.log10(np.median(top_values)) - np.log10(np.median(other_values)):.2f}"
        else:
            median_text = f"{np.median(top_values):.3g} vs {np.median(other_values):.3g}"
            delta_text = f"delta median {delta_median:.3g}"
        axis.text(
            0.98,
            0.06,
            f"Top vs other median:\n{median_text}\n{delta_text}",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    for axis in axes[n_features:]:
        axis.axis("off")

    figure.suptitle(
        f"{query_condition} residual-high introns vs other shared introns\n"
        f"Top residual decile defined relative to baseline from {control_condition}",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def main():
    args = parse_args()
    if not (0.0 < args.residual_quantile_fraction <= 0.5):
        raise ValueError("--residual-quantile-fraction must be in (0, 0.5]")

    condition_records = load_condition_rows(args.site_counts, args.signal_field)
    if args.control_condition not in condition_records:
        raise ValueError(f"Missing control condition {args.control_condition}")
    if args.query_condition not in condition_records:
        raise ValueError(f"Missing query condition {args.query_condition}")

    shared_rows = build_shared_rows(
        condition_records,
        args.genome_fasta,
        args.control_condition,
        args.query_condition,
        args.control_min_reads,
        args.query_min_reads,
        args.require_bp_center_base.upper() if args.require_bp_center_base else "",
        args.bp_upstream,
        args.bp_downstream,
        args.intron_gc_window,
        args.three_prime_window,
    )

    numeric_features = [
        "anchored_enrichment_log2_fold_change",
        "log_intron_length",
        "branchpoint_to_3ss_nt",
        "branchpoint_score",
        "branchpoint_candidates",
        "donor_maxent",
        "acceptor_maxent",
        "intron_gc",
        "intron_gc_5p_window",
        "intron_gc_3p_window",
        "bp_to_3ss_pyrimidine_fraction",
        "three_prime_window_pyrimidine_fraction",
    ]
    shared_rows = finite_feature_rows(shared_rows, numeric_features)
    if not shared_rows:
        raise ValueError("No shared introns remained after feature filtering")

    bp_positions = list(range(-args.bp_upstream, args.bp_downstream + 1))
    successes = np.array([row["query_branched_fragments"] for row in shared_rows], dtype=float)
    totals = np.array([row["query_anchored_fragments"] for row in shared_rows], dtype=float)
    overall_query_probability = float(successes.sum() / totals.sum())
    overall_query_smoothed = empirical_probability(float(successes.sum()), float(totals.sum()))
    null_probabilities = np.full(len(shared_rows), overall_query_probability, dtype=float)
    null_deviance = binomial_deviance(successes, totals, null_probabilities)

    null_model = {
        "label": "null_intercept_only",
        "design_matrix": np.ones((len(shared_rows), 1), dtype=float),
        "metadata": [{"kind": "intercept"}],
        "beta": np.array([logit(overall_query_smoothed)], dtype=float),
        "covariance": np.array([[0.0]], dtype=float),
        "predicted_fraction": null_probabilities,
        "residual_deviance": null_deviance,
        "weighted_r_squared": 0.0,
    }
    baseline_model = fit_model(shared_rows, label="baseline_control_only")
    numeric_model = fit_model(shared_rows, numeric_features=numeric_features, label="baseline_plus_numeric")
    bp_model = fit_model(shared_rows, bp_positions=bp_positions, label="baseline_plus_bp_flank")
    full_model = fit_model(
        shared_rows,
        numeric_features=numeric_features,
        bp_positions=bp_positions,
        label="baseline_plus_numeric_plus_bp_flank",
    )

    model_lookup = {
        "null_intercept_only": null_model,
        "baseline_control_only": baseline_model,
        "baseline_plus_numeric": numeric_model,
        "baseline_plus_bp_flank": bp_model,
        "baseline_plus_numeric_plus_bp_flank": full_model,
    }
    model_summary_rows = build_model_summary(
        shared_rows,
        args.control_condition,
        args.query_condition,
        model_lookup,
        null_deviance,
        baseline_model["residual_deviance"],
    )

    baseline_summary_row = next(row for row in model_summary_rows if row["model"] == "baseline_control_only")
    baseline_summary_row["intercept_beta"] = float(baseline_model["beta"][0])
    baseline_summary_row["control_logit_beta"] = float(baseline_model["beta"][1])
    for row in model_summary_rows:
        if row["model"] != "baseline_control_only":
            row["intercept_beta"] = ""
            row["control_logit_beta"] = ""

    feature_summary_rows = build_numeric_feature_summary(shared_rows, numeric_features, null_deviance, baseline_model)
    add_baseline_predictions(shared_rows, baseline_model)
    bp_effect_rows = build_bp_sequence_effect_rows(shared_rows, bp_positions)
    bp_enrichment_rows, top_rows, bottom_rows = build_bp_sequence_enrichment_rows(
        shared_rows,
        bp_positions,
        args.residual_quantile_fraction,
    )
    feature_specs = distribution_specs(feature_summary_rows)
    group_comparison_rows = build_group_comparison_rows(shared_rows, top_rows, feature_specs)

    shared_rows.sort(
        key=lambda row: (
            -row["baseline_residual_query_branch_fraction"],
            -row["query_anchored_fragments"],
            row["gene_name"],
            row["intron_id"],
        )
    )
    model_summary_rows.sort(key=lambda row: MODEL_ORDER.get(row["model"], 999))
    bp_effect_rows.sort(key=lambda row: (row["position"], row["base"]))
    bp_enrichment_rows.sort(key=lambda row: (row["position"], row["base"]))
    group_comparison_rows.sort(key=lambda row: [spec["model_feature"] for spec in feature_specs].index(row["model_feature"]))

    best_numeric_feature = feature_summary_rows[0]["feature"] if feature_summary_rows else ""
    best_numeric_delta = feature_summary_rows[0]["added_fraction_of_baseline_residual"] if feature_summary_rows else 0.0
    full_model_row = next(row for row in model_summary_rows if row["model"] == "baseline_plus_numeric_plus_bp_flank")
    summary_rows = [
        {
            "control_condition": args.control_condition,
            "query_condition": args.query_condition,
            "signal_field": args.signal_field,
            "require_bp_center_base": args.require_bp_center_base or "",
            "introns_analysed": len(shared_rows),
            "control_anchored_fragments": sum(row["control_anchored_fragments"] for row in shared_rows),
            "control_branched_fragments": sum(row["control_branched_fragments"] for row in shared_rows),
            "control_branch_fraction": overall_fraction(shared_rows, "control_branched_fragments", "control_anchored_fragments"),
            "query_anchored_fragments": sum(row["query_anchored_fragments"] for row in shared_rows),
            "query_branched_fragments": sum(row["query_branched_fragments"] for row in shared_rows),
            "query_branch_fraction": overall_fraction(shared_rows, "query_branched_fragments", "query_anchored_fragments"),
            "baseline_total_deviance_explained": baseline_summary_row["total_deviance_explained"],
            "baseline_weighted_r_squared": baseline_summary_row["weighted_r_squared"],
            "full_total_deviance_explained": full_model_row["total_deviance_explained"],
            "full_weighted_r_squared": full_model_row["weighted_r_squared"],
            "full_added_fraction_of_baseline_residual": full_model_row["added_fraction_of_baseline_residual"],
            "best_numeric_feature": best_numeric_feature,
            "best_numeric_feature_added_fraction_of_baseline_residual": best_numeric_delta,
        }
    ]

    write_rows(args.output_summary, summary_rows, list(summary_rows[0].keys()))
    write_rows(args.output_shared_introns, shared_rows, list(shared_rows[0].keys()))
    write_rows(args.output_model_summary, model_summary_rows, list(model_summary_rows[0].keys()))
    write_rows(args.output_feature_summary, feature_summary_rows, list(feature_summary_rows[0].keys()))
    write_rows(args.output_bp_sequence_effects, bp_effect_rows, list(bp_effect_rows[0].keys()))
    write_rows(args.output_bp_sequence_enrichment, bp_enrichment_rows, list(bp_enrichment_rows[0].keys()))
    write_rows(args.output_group_comparison, group_comparison_rows, list(group_comparison_rows[0].keys()))

    plot_results(
        shared_rows,
        args.control_condition,
        args.query_condition,
        args.signal_field,
        bp_positions,
        model_summary_rows,
        feature_summary_rows,
        bp_effect_rows,
        bp_enrichment_rows,
        top_rows,
        bottom_rows,
        args.output_plot_png,
        args.output_plot_pdf,
    )
    plot_group_feature_distributions(
        shared_rows,
        top_rows,
        feature_specs,
        args.query_condition,
        args.control_condition,
        args.output_feature_distribution_plot_png,
        args.output_feature_distribution_plot_pdf,
    )

    print(f"Shared introns analysed: {len(shared_rows)}")
    print(f"Best numeric feature: {best_numeric_feature}")
    print(
        f"Full model added residual deviance over baseline: "
        f"{full_model_row['added_fraction_of_baseline_residual'] * 100.0:.2f}%"
    )


if __name__ == "__main__":
    main()
