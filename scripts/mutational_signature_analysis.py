#!/usr/bin/env python3
"""Exploratory COSMIC SBS signature analysis for the local targeted-panel MAF.

This script deliberately does not install packages or reference data. It uses a
preconfigured SigProfiler environment and keeps all sample-level results under
the private output directory. Public outputs contain aggregate counts only.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import fisher_exact, mannwhitneyu

from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as mat_gen


PROCESSES = {
    "Platinum-associated SBS": ["SBS31", "SBS35"],
    "HRD-associated SBS": ["SBS3"],
    "MMR-associated SBS": ["SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44"],
}
SIGNATURE_TO_PROCESS = {
    signature: process for process, signatures in PROCESSES.items() for signature in signatures
}
TARGET_SIGNATURES = list(SIGNATURE_TO_PROCESS)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--maf", required=True)
    p.add_argument("--cohort", required=True)
    p.add_argument("--clinical", required=True)
    p.add_argument("--clinical-sheet", default="Main")
    p.add_argument("--output", required=True, help="Private sample-level output directory")
    p.add_argument("--public-output", required=True, help="Manuscript output root")
    p.add_argument("--reference-volume", required=True)
    p.add_argument("--min-snv", type=int, default=10)
    p.add_argument("--bootstraps", type=int, default=100)
    p.add_argument("--min-count", type=int, default=5)
    p.add_argument("--min-proportion", type=float, default=0.20)
    p.add_argument("--min-cosine", type=float, default=0.90)
    p.add_argument("--min-stability", type=float, default=0.80)
    p.add_argument("--seed", type=int, default=20260928)
    return p.parse_args()


def holm(values: pd.Series) -> pd.Series:
    p = np.asarray(values, dtype=float)
    result = np.full(len(p), np.nan)
    keep = np.isfinite(p)
    idx = np.where(keep)[0]
    if not len(idx):
        return pd.Series(result, index=values.index)
    order = idx[np.argsort(p[idx])]
    adjusted = np.maximum.accumulate(np.minimum(1.0, p[order] * (len(order) - np.arange(len(order)))))
    result[order] = adjusted
    return pd.Series(result, index=values.index)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denom) if denom > 0 else np.nan


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def maf_audit(maf: pd.DataFrame, cohort_ids: set[str], args: argparse.Namespace) -> pd.DataFrame:
    is_snv = maf["Variant_Type"].astype(str).eq("SNV")
    silent = maf["Variant_Classification"].astype(str).str.contains("Silent|Synonymous", case=False, regex=True)
    cohort = maf[maf["Tumor_Sample_Barcode"].astype(str).isin(cohort_ids)]
    cohort_snv = cohort[cohort["Variant_Type"].astype(str).eq("SNV")]
    per_sample = cohort_snv.groupby("Tumor_Sample_Barcode").size().reindex(sorted(cohort_ids), fill_value=0)
    rows = [
        ("Clinical workbook", str(Path(args.clinical).name), "Main analysis input"),
        ("Clinical sheet", args.clinical_sheet, "Configured sheet"),
        ("MAF file", str(Path(args.maf).name), "Genomic input"),
        ("Genome build", ", ".join(sorted(maf["NCBI_Build"].astype(str).unique())), "MAF field"),
        ("MAF rows", len(maf), "All retained variants"),
        ("MAF samples", maf["Tumor_Sample_Barcode"].nunique(), "All MAF samples"),
        ("Analysis cohort samples", len(cohort_ids), "After duplicate-patient exclusion"),
        ("Cohort samples represented in MAF", cohort["Tumor_Sample_Barcode"].nunique(), "Matched by Tumor_Sample_Barcode"),
        ("SNV rows in full MAF", int(is_snv.sum()), "Sequence-context eligible before reference checks"),
        ("Synonymous or silent rows", int(silent.sum()), "Zero indicates a retained-coding-variant input"),
        ("Median cohort SNVs per sample", float(per_sample.median()), "Retained MAF SNVs, not mutations/Mb"),
        (f"Cohort samples with >= {args.min_snv} raw MAF SNVs", int((per_sample >= args.min_snv).sum()), "Before exome-context generation"),
        ("Panel BED/callable territory", "Not available", "No mutations/Mb signature burden is calculated"),
        ("Interpretation", "Exploratory targeted-panel signature analysis", "Not a definitive signature diagnosis"),
    ]
    return pd.DataFrame(rows, columns=["item", "value", "definition"])


def assignment_activities(assignment_dir: Path) -> pd.DataFrame:
    hits = list(assignment_dir.glob("**/Assignment_Solution_Activities.txt"))
    if not hits:
        raise FileNotFoundError(f"No Assignment_Solution_Activities.txt found below {assignment_dir}")
    act = pd.read_csv(hits[0], sep="\t")
    sample_col = act.columns[0]
    return act.rename(columns={sample_col: "sample_id"}).set_index("sample_id")


def odds_ratio_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    cells = np.array([a, b, c, d], dtype=float)
    if np.any(cells == 0):
        cells += 0.5
    aa, bb, cc, dd = cells
    estimate = aa * dd / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return estimate, math.exp(math.log(estimate) - 1.96 * se), math.exp(math.log(estimate) + 1.96 * se)


def rank_biserial(x_yes: np.ndarray, x_no: np.ndarray, u_yes: float) -> float:
    return 2 * u_yes / (len(x_yes) * len(x_no)) - 1 if len(x_yes) and len(x_no) else np.nan


def pooled_sensitivity_analysis(
    sbs96: pd.DataFrame, cohort: pd.DataFrame, reference: pd.DataFrame,
    private: Path, args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    groups = {}
    for group in ["No", "Yes"]:
        ids = cohort.loc[cohort["Preop_platinum_exposure"].eq(group), "sample_id"].astype(str)
        groups[group] = [x for x in ids if x in sbs96.columns]
    pooled = pd.DataFrame({g: sbs96[ids].sum(axis=1) for g, ids in groups.items()})
    pooled_file = private / "SBS96_pooled_by_platinum.tsv"
    pooled.to_csv(pooled_file, sep="\t", index_label="MutationType")
    pooled_assignment = private / "pooled_sigprofiler_assignment"
    if pooled_assignment.exists():
        shutil.rmtree(pooled_assignment)
    Analyze.cosmic_fit(
        samples=str(pooled_file), output=str(pooled_assignment), genome_build="GRCh37",
        cosmic_version=3.6, make_plots=False, verbose=False, exome=True,
        input_type="matrix", context_type="96", export_probabilities=False,
        export_probabilities_per_mutation=False, sample_reconstruction_plots=False,
        cpu=2, add_background_signatures=True,
    )
    baseline = assignment_activities(pooled_assignment).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    active = [s for s in baseline.columns if baseline[s].sum() > 0 and s in reference.columns]
    targets = sorted(set(TARGET_SIGNATURES) | {"SBS1", "SBS5"})
    candidates = list(dict.fromkeys(active + [s for s in targets if s in reference.columns]))
    sigmat = reference[candidates].to_numpy(float)
    entities = {**PROCESSES, **{s: [s] for s in TARGET_SIGNATURES}}

    def fit_spectrum(x: np.ndarray, n_patients: int) -> dict[str, dict[str, float]]:
        coef, _ = nnls(sigmat, x)
        rec = sigmat @ coef
        total = float(x.sum())
        result = {}
        for entity, signatures in entities.items():
            count = float(sum(coef[candidates.index(s)] for s in signatures if s in candidates))
            result[entity] = {
                "burden_per_patient": count / n_patients,
                "proportion": count / total if total else np.nan,
                "attributed_count": count,
                "cosine": cosine(x, rec),
            }
        return result

    observed = {g: fit_spectrum(pooled[g].to_numpy(float), len(groups[g])) for g in groups}
    reps = max(1000, args.bootstraps)
    rng = np.random.default_rng(args.seed + 101)
    boot = {e: {m: {g: [] for g in groups} for m in ["burden_per_patient", "proportion"]} for e in entities}
    for group, ids in groups.items():
        matrix = sbs96[ids].to_numpy(float)
        n = len(ids)
        for _ in range(reps):
            weights = rng.multinomial(n, np.repeat(1 / n, n))
            fitted = fit_spectrum(matrix @ weights, n)
            for entity in entities:
                for metric in ["burden_per_patient", "proportion"]:
                    boot[entity][metric][group].append(fitted[entity][metric])

    all_ids = groups["No"] + groups["Yes"]
    all_matrix = sbs96[all_ids].to_numpy(float)
    n_yes = len(groups["Yes"]); n_no = len(groups["No"])
    perm_diffs = {e: {m: [] for m in ["burden_per_patient", "proportion"]} for e in entities}
    for _ in range(reps):
        order = rng.permutation(len(all_ids))
        yes_x = all_matrix[:, order[:n_yes]].sum(axis=1)
        no_x = all_matrix[:, order[n_yes:]].sum(axis=1)
        fy = fit_spectrum(yes_x, n_yes); fn = fit_spectrum(no_x, n_no)
        for entity in entities:
            for metric in ["burden_per_patient", "proportion"]:
                perm_diffs[entity][metric].append(fy[entity][metric] - fn[entity][metric])

    def result_rows(selected: dict[str, list[str]], entity_column: str) -> pd.DataFrame:
        rows = []
        label = {"burden_per_patient": "Pooled attributed retained SNVs per patient",
                 "proportion": "Pooled signature proportion"}
        for entity, signatures in selected.items():
            for metric in ["burden_per_patient", "proportion"]:
                no_value = observed["No"][entity][metric]
                yes_value = observed["Yes"][entity][metric]
                diff = yes_value - no_value
                boot_diff = np.asarray(boot[entity][metric]["Yes"]) - np.asarray(boot[entity][metric]["No"])
                null = np.asarray(perm_diffs[entity][metric])
                p = (1 + np.sum(np.abs(null) >= abs(diff))) / (reps + 1)
                row = {
                    entity_column: entity, "metric": label[metric],
                    "method": f"Patient bootstrap CI and label permutation test ({reps} replicates)",
                    "no_n": n_no, "yes_n": n_yes, "no_value": no_value, "yes_value": yes_value,
                    "difference_yes_minus_no": diff,
                    "difference_lower95": np.quantile(boot_diff, .025),
                    "difference_upper95": np.quantile(boot_diff, .975),
                    "no_lower95": np.quantile(boot[entity][metric]["No"], .025),
                    "no_upper95": np.quantile(boot[entity][metric]["No"], .975),
                    "yes_lower95": np.quantile(boot[entity][metric]["Yes"], .025),
                    "yes_upper95": np.quantile(boot[entity][metric]["Yes"], .975),
                    "permutation_p": p,
                    "no_constrained_cosine": observed["No"][entity]["cosine"],
                    "yes_constrained_cosine": observed["Yes"][entity]["cosine"],
                    "no_standard_assignment_count": sum(baseline.loc["No"].get(s, 0.0) for s in signatures),
                    "yes_standard_assignment_count": sum(baseline.loc["Yes"].get(s, 0.0) for s in signatures),
                    "interpretation_scope": "Group-pooled hypothesis-directed sensitivity analysis; not a patient-level positive rate or treatment-induced effect",
                }
                if entity_column == "signature":
                    row["biological_process"] = SIGNATURE_TO_PROCESS[entity]
                rows.append(row)
        result = pd.DataFrame(rows)
        correction_n = len(result)
        correction_col = f"p_holm_{correction_n}_tests"
        result[correction_col] = holm(result["permutation_p"])
        result["p_display"] = result["permutation_p"].map(
            lambda x: "<0.001" if x < .001 else f"{x:.3f}")
        result["adjusted_p_display"] = result[correction_col].map(
            lambda x: "<0.001" if x < .001 else f"{x:.3f}")
        return result

    process_result = result_rows(PROCESSES, "process")
    signature_result = result_rows({s: [s] for s in TARGET_SIGNATURES}, "signature")
    audit_rows = []
    for signature in TARGET_SIGNATURES:
        audit_rows.append({
            "signature": signature,
            "biological_process": SIGNATURE_TO_PROCESS[signature],
            "no_standard_assignment_count": baseline.loc["No"].get(signature, 0.0),
            "no_constrained_refit_count": observed["No"][signature]["attributed_count"],
            "yes_standard_assignment_count": baseline.loc["Yes"].get(signature, 0.0),
            "yes_constrained_refit_count": observed["Yes"][signature]["attributed_count"],
            "interpretation": "Constrained values are hypothesis-directed sensitivity estimates and do not establish patient-level signature positivity",
        })
    assignment_audit = pd.DataFrame(audit_rows)
    pd.DataFrame({"candidate_signature": candidates, "source": ["Standard pooled assignment or prespecified target/background"] * len(candidates)}).to_csv(
        private / "pooled_constrained_refit_signatures.tsv", sep="\t", index=False)
    return process_result, signature_result, assignment_audit


def requested_metric_status(tests: pd.DataFrame) -> pd.DataFrame:
    rows = []
    display = {
        "Signature-positive proportion": "Positive rate",
        "Signature-attributed retained SNV count": "Patient-level burden",
        "Signature proportion": "Patient-level proportion",
    }
    for _, item in tests.iterrows():
        evaluable = int(item["no_n"]) + int(item["yes_n"])
        estimable = int(item["no_n"]) > 0 and int(item["yes_n"]) > 0
        rows.append({
            "process": item["process"], "requested_metric": display[item["metric"]],
            "no_evaluable_n": int(item["no_n"]), "yes_evaluable_n": int(item["yes_n"]),
            "total_evaluable_n": evaluable,
            "status": "Estimated" if estimable else "Not estimable",
            "result": (f"No={item['no_value']:.4g}; Yes={item['yes_value']:.4g}"
                       if estimable else "NE"),
            "reason": ("Both exposure groups contain quality-evaluable patients"
                       if estimable else
                       "No patient passed both the retained-SNV and reconstruction-cosine criteria"),
            "quality_rule": "Retained SNVs >= minimum and reconstruction cosine >= minimum; positive classification additionally requires count, proportion, and bootstrap stability thresholds",
        })
    return pd.DataFrame(rows)


def make_public_plots(samples: pd.DataFrame, tests: pd.DataFrame, pooled: pd.DataFrame,
                      individual: pd.DataFrame, status: pd.DataFrame,
                      assignment_audit: pd.DataFrame, figures: Path,
                      args: argparse.Namespace) -> None:
    plt.rcParams.update({"font.family": "Pretendard", "axes.titlesize": 12, "axes.labelsize": 10})
    colors = {"No": "#A8A8A8", "Yes": "#246A8A"}
    processes = list(PROCESSES)
    fig, axes = plt.subplots(len(processes), 2, figsize=(11, 12))
    for row, process in enumerate(processes):
        for col, metric, ylabel in [(0, "Pooled attributed retained SNVs per patient", "Attributed retained SNVs per patient"),
                                    (1, "Pooled signature proportion", "Proportion of retained SNVs")]:
            ax = axes[row, col]
            q = pooled[(pooled["process"] == process) & (pooled["metric"] == metric)].iloc[0]
            vals = [q["no_value"], q["yes_value"]]
            lows = [q["no_lower95"], q["yes_lower95"]]
            highs = [q["no_upper95"], q["yes_upper95"]]
            yerr = np.array([[v - lo for v, lo in zip(vals, lows)], [hi - v for v, hi in zip(vals, highs)]])
            ax.bar(["No", "Yes"], vals, color=[colors["No"], colors["Yes"]], yerr=yerr, capsize=4, alpha=.85)
            ax.text(.02, .98, f"Permutation p={q['p_display']}\nHolm p={q['adjusted_p_display']}\nDifference Yes-No={q['difference_yes_minus_no']:.3f}",
                    transform=ax.transAxes, va="top", fontsize=9)
            ax.set_ylabel(ylabel)
            ax.set_title(process)
            if col == 1:
                ax.set_ylim(bottom=0)
    fig.suptitle("Group-pooled hypothesis-directed SBS refit", fontsize=16, fontweight="bold")
    pooled_note = (
        "Bars are group-pooled constrained refit estimates with patient-bootstrap 95% intervals. "
        "P values use patient-label permutation and Holm correction across six tests. This sensitivity "
        "analysis is not a patient-level signature-positive analysis. Pre-treatment sampling precludes "
        "a treatment-induced interpretation."
    )
    fig.text(.01, .008, textwrap.fill(pooled_note, width=165), fontsize=9, va="bottom")
    fig.tight_layout(rect=[0, .075, 1, .96], h_pad=2.2, w_pad=1.7)
    fig.savefig(figures / "Supplementary_Figure_Platinum_HRD_MMR_signatures.png", dpi=300, facecolor="white")
    plt.close(fig)

    signature_order = TARGET_SIGNATURES
    process_colors = {
        "Platinum-associated SBS": "#C17438",
        "HRD-associated SBS": "#246A8A",
        "MMR-associated SBS": "#438569",
    }
    fig, axes = plt.subplots(1, 2, figsize=(13, 7.8), sharey=True)
    for ax, metric, title in [
        (axes[0], "Pooled attributed retained SNVs per patient", "Attributed retained SNVs per patient"),
        (axes[1], "Pooled signature proportion", "Proportion of retained SNVs"),
    ]:
        q = individual[individual["metric"].eq(metric)].set_index("signature").loc[signature_order]
        y = np.arange(len(q))
        vals = q["difference_yes_minus_no"].to_numpy(float)
        low = q["difference_lower95"].to_numpy(float)
        high = q["difference_upper95"].to_numpy(float)
        colors = [process_colors[p] for p in q["biological_process"]]
        ax.axvline(0, color="#777777", linestyle="--", linewidth=1)
        ax.hlines(y, low, high, color=colors, linewidth=1.4)
        ax.scatter(vals, y, c=colors, s=42, zorder=3)
        for yi, (_, row) in zip(y, q.iterrows()):
            ax.text(high[yi], yi, f"  Holm p={row['adjusted_p_display']}", va="center", fontsize=7.6)
        ax.set_title(title)
        ax.set_xlabel("Difference: preoperative platinum Yes - No")
        ax.set_yticks(y, signature_order)
        ax.invert_yaxis()
        ax.grid(axis="x", color="#E5E5E5", linewidth=.6)
    fig.suptitle("Individual COSMIC SBS pooled sensitivity analysis", fontsize=15, fontweight="bold")
    fig.text(.01, .012, textwrap.fill(
        "Each point is a group-pooled constrained-refit difference with a patient-bootstrap 95% interval. "
        "Holm correction covers 20 prespecified signature-metric tests. These are sensitivity estimates, "
        "not patient-level signature positivity or evidence that preoperative platinum induced the signatures.",
        width=165), fontsize=9, va="bottom")
    fig.tight_layout(rect=[0, .09, 1, .95], w_pad=2.5)
    fig.savefig(figures / "Supplementary_Figure_Individual_SBS_Pooled_Sensitivity.png", dpi=300, facecolor="white")
    plt.close(fig)

    audit = assignment_audit.set_index("signature").loc[signature_order]
    audit_matrix = audit[["no_standard_assignment_count", "no_constrained_refit_count",
                          "yes_standard_assignment_count", "yes_constrained_refit_count"]].to_numpy(float)
    fig, ax = plt.subplots(figsize=(10.5, 7.6))
    image = ax.imshow(np.log1p(audit_matrix), cmap="Blues", aspect="auto")
    for i in range(audit_matrix.shape[0]):
        for j in range(audit_matrix.shape[1]):
            value = audit_matrix[i, j]
            ax.text(j, i, f"{value:.1f}", ha="center", va="center",
                    color="white" if np.log1p(value) > np.log1p(audit_matrix).max() * .55 else "#222222",
                    fontsize=9)
    ax.set_xticks(range(4), ["No\nstandard", "No\nconstrained", "Yes\nstandard", "Yes\nconstrained"])
    ax.set_yticks(range(len(signature_order)), signature_order)
    ax.set_title("Target SBS assignment depends on hypothesis-directed constrained refit", fontweight="bold")
    cbar = fig.colorbar(image, ax=ax, fraction=.025, pad=.03)
    cbar.set_label("log(1 + pooled attributed SNV count)")
    fig.text(.01, .012, textwrap.fill(
        "Cell labels are pooled attributed SNV counts. Standard SigProfilerAssignment and constrained refit "
        "answer different questions; a signal appearing only after forced inclusion is exploratory and cannot "
        "be used as a patient-level positive call.", width=145), fontsize=9, va="bottom")
    fig.tight_layout(rect=[0, .085, 1, .96])
    fig.savefig(figures / "Supplementary_Figure_Signature_Assignment_Audit.png", dpi=300, facecolor="white")
    plt.close(fig)

    metric_order = ["Positive rate", "Patient-level burden", "Patient-level proportion"]
    process_order = list(PROCESSES)
    status_grid = status.pivot(index="process", columns="requested_metric", values="total_evaluable_n").loc[process_order, metric_order]
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.imshow(status_grid.to_numpy(float) > 0, cmap=plt.matplotlib.colors.ListedColormap(["#E7E7E7", "#7BAAC1"]), vmin=0, vmax=1)
    for i, process in enumerate(process_order):
        for j, metric in enumerate(metric_order):
            n_eval = int(status_grid.loc[process, metric])
            ax.text(j, i, f"{'Estimated' if n_eval else 'NE'}\n{n_eval} evaluable", ha="center", va="center", fontsize=10)
    ax.set_xticks(range(3), ["Positive rate", "Patient-level\nburden", "Patient-level\nproportion"])
    ax.set_yticks(range(3), process_order)
    ax.set_title("Availability of the requested patient-level signature comparisons", fontweight="bold")
    fig.text(.01, .015, textwrap.fill(
        "NE means not estimable under the prespecified quality criteria; it does not mean that the signature "
        "is biologically absent. The thresholds were not relaxed after observing the results.", width=140),
        fontsize=9, va="bottom")
    fig.tight_layout(rect=[0, .12, 1, .94])
    fig.savefig(figures / "Supplementary_Figure_Signature_Requested_Metrics.png", dpi=300, facecolor="white")
    plt.close(fig)

    qc = samples.drop_duplicates("sample_id")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8))
    bins = np.arange(0, max(20, min(65, int(qc["total_snv"].max()) + 2)), 2)
    axes[0].hist(qc["total_snv"], bins=bins, color="#246A8A", alpha=.8)
    axes[0].axvline(args.min_snv, color="#C17438", linestyle="--", label=f"Minimum={args.min_snv}")
    axes[0].set_xlabel("Retained SNVs per patient")
    axes[0].set_ylabel("Patients")
    axes[0].set_title("Available SNV information")
    axes[0].legend(frameon=False)
    eligible = qc[qc["total_snv"] >= args.min_snv]
    axes[1].hist(eligible["reconstruction_cosine"].dropna(), bins=np.linspace(.35, .92, 18), color="#438569", alpha=.85)
    axes[1].axvline(args.min_cosine, color="#C17438", linestyle="--", label=f"Minimum={args.min_cosine:.2f}")
    axes[1].set_xlabel("COSMIC reconstruction cosine")
    axes[1].set_ylabel("Patients with retained SNVs ≥ threshold")
    axes[1].set_title("No patient met both quality criteria")
    axes[1].legend(frameon=False)
    fig.suptitle("Mutational signature feasibility and quality control", fontsize=15, fontweight="bold")
    qc_note = (
        f"{len(eligible)}/{len(qc)} patients had at least {args.min_snv} retained SNVs; maximum reconstruction "
        f"cosine={eligible['reconstruction_cosine'].max():.3f}. Therefore patient-level positive rate, burden, "
        "and proportion comparisons are not estimable under the prespecified criteria."
    )
    fig.text(.01, .008, textwrap.fill(qc_note, width=125), fontsize=9, va="bottom")
    fig.tight_layout(rect=[0, .14, 1, .93], w_pad=2.5)
    fig.savefig(figures / "Supplementary_Figure_Signature_QC.png", dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    private = Path(args.output).resolve()
    public = Path(args.public_output).resolve()
    tables = public / "tables"
    figures = public / "figures"
    private.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    cohort_ids = set(cohort["sample_id"].dropna().astype(str))
    maf = pd.read_csv(args.maf, sep="\t", comment="#", low_memory=False)
    required = {"Tumor_Sample_Barcode", "NCBI_Build", "Chromosome", "Start_Position", "End_Position",
                "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Variant_Type", "Variant_Classification"}
    missing = sorted(required - set(maf.columns))
    if missing:
        raise ValueError(f"Required MAF fields missing: {missing}")
    audit = maf_audit(maf, cohort_ids, args)
    write_tsv(audit, tables / "Input_MAF_signature_audit.tsv")

    input_dir = private / "matrix_input"
    if input_dir.exists():
        shutil.rmtree(input_dir)
    input_dir.mkdir(parents=True)
    # SBS96 uses single-base substitutions only. Explicitly exclude DNP/TNP,
    # indels and other multi-base events rather than relying on converter rules.
    cohort_maf = maf[
        maf["Tumor_Sample_Barcode"].astype(str).isin(cohort_ids)
        & maf["Variant_Type"].astype(str).eq("SNV")
    ].copy()
    cohort_maf.to_csv(input_dir / "cohort.maf", sep="\t", index=False)

    matrices = mat_gen.SigProfilerMatrixGeneratorFunc(
        "PAAD_targeted_panel", "GRCh37", str(input_dir), exome=True, plot=False,
        tsb_stat=False, seqInfo=True, volume=str(Path(args.reference_volume).resolve())
    )
    sbs96 = matrices["96"].copy()
    sbs96.index = sbs96.index.astype(str)
    sbs96.columns = sbs96.columns.astype(str)
    matrix_file = private / "SBS96_all_samples.tsv"
    sbs96.to_csv(matrix_file, sep="\t", index_label="MutationType")

    all_ids = sorted(cohort_ids)
    total_snv = sbs96.sum(axis=0).reindex(all_ids, fill_value=0).astype(int)
    eligible_ids = total_snv[total_snv >= args.min_snv].index.tolist()
    audit = pd.concat([audit, pd.DataFrame([
        ("SBS96 context-counted cohort SNVs", int(total_snv.sum()), "After GRCh37 context validation and exome filtering"),
        (f"Cohort samples with >= {args.min_snv} context-counted SNVs", len(eligible_ids), "Actual SigProfilerAssignment input"),
    ], columns=audit.columns)], ignore_index=True)
    write_tsv(audit, tables / "Input_MAF_signature_audit.tsv")
    if not eligible_ids:
        raise RuntimeError("No cohort samples met the minimum retained-SNV threshold.")
    eligible_matrix = sbs96.reindex(columns=eligible_ids)
    eligible_file = private / "SBS96_eligible.tsv"
    eligible_matrix.to_csv(eligible_file, sep="\t", index_label="MutationType")

    assignment_dir = private / "sigprofiler_assignment"
    if assignment_dir.exists():
        shutil.rmtree(assignment_dir)
    Analyze.cosmic_fit(
        samples=str(eligible_file), output=str(assignment_dir), genome_build="GRCh37",
        cosmic_version=3.6, make_plots=False, verbose=False, exome=True,
        input_type="matrix", context_type="96", export_probabilities=False,
        export_probabilities_per_mutation=False, sample_reconstruction_plots=False,
        cpu=4, add_background_signatures=True,
    )
    activities = assignment_activities(assignment_dir)
    activities.index = activities.index.astype(str)
    activities = activities.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    ref_path = Path(__import__("SigProfilerAssignment").__file__).resolve().parent / \
        "data" / "Reference_Signatures" / "GRCh37" / "COSMIC_v3.6_SBS_GRCh37_exome.txt"
    reference = pd.read_csv(ref_path, sep="\t").set_index("Type")
    common = [x for x in eligible_matrix.index if x in reference.index]
    if len(common) != 96:
        raise RuntimeError(f"Expected 96 aligned SBS contexts, found {len(common)}")
    reference = reference.loc[common].astype(float)
    observed = eligible_matrix.loc[common].astype(float)

    target_signatures = sorted({s for sigs in PROCESSES.values() for s in sigs} | {"SBS1", "SBS5"})
    sample_rows: list[dict] = []
    rng = np.random.default_rng(args.seed)
    for sample_id in all_ids:
        total = int(total_snv.get(sample_id, 0))
        exposure = activities.loc[sample_id] if sample_id in activities.index else pd.Series(dtype=float)
        active = [s for s, value in exposure.items() if value > 0 and s in reference.columns]
        candidates = list(dict.fromkeys(active + [s for s in target_signatures if s in reference.columns]))
        x = observed[sample_id].to_numpy(float) if sample_id in observed.columns else np.zeros(96)
        reconstruction = np.zeros(96)
        for sig, value in exposure.items():
            if sig in reference.columns:
                reconstruction += reference[sig].to_numpy(float) * float(value)
        cos = cosine(x, reconstruction) if total >= args.min_snv else np.nan
        quality = total >= args.min_snv and np.isfinite(cos) and cos >= args.min_cosine
        quality_status = "Evaluable" if quality else ("Low retained-SNV count" if total < args.min_snv else "Low reconstruction cosine")

        bootstrap_hits = {process: 0 for process in PROCESSES}
        if quality and candidates and total > 0:
            sigmat = reference[candidates].to_numpy(float)
            probs = x / x.sum()
            for _ in range(args.bootstraps):
                xb = rng.multinomial(total, probs)
                coef, _ = nnls(sigmat, xb)
                recon_b = sigmat @ coef
                cos_b = cosine(xb, recon_b)
                for process, sigs in PROCESSES.items():
                    count_b = sum(coef[candidates.index(s)] for s in sigs if s in candidates)
                    prop_b = count_b / total if total else 0
                    if count_b >= args.min_count and prop_b >= args.min_proportion and cos_b >= args.min_cosine:
                        bootstrap_hits[process] += 1
        for process, sigs in PROCESSES.items():
            count = float(sum(exposure.get(s, 0.0) for s in sigs))
            proportion = count / total if total else np.nan
            stability = bootstrap_hits[process] / args.bootstraps if quality else np.nan
            if not quality:
                classification = "Unclassifiable"
            elif count >= args.min_count and proportion >= args.min_proportion:
                classification = "Positive" if stability >= args.min_stability else "Indeterminate"
            else:
                classification = "Not positive"
            sample_rows.append({
                "sample_id": sample_id, "process": process, "total_snv": total,
                "reconstruction_cosine": cos, "quality_evaluable": quality,
                "quality_status": quality_status, "attributed_count": count,
                "proportion": proportion, "bootstrap_stability": stability,
                "classification": classification,
            })
    sample_results = pd.DataFrame(sample_rows).merge(
        cohort[["sample_id", "Preop_platinum_exposure"]], on="sample_id", how="left"
    )
    write_tsv(sample_results, private / "signature_sample_results.tsv")

    quality_summary = (sample_results.drop_duplicates("sample_id")
        .groupby(["Preop_platinum_exposure", "quality_status"], dropna=False).size()
        .rename("n").reset_index())
    quality_summary["percent_within_exposure"] = quality_summary["n"] / quality_summary.groupby("Preop_platinum_exposure")["n"].transform("sum") * 100
    write_tsv(quality_summary, tables / "Input_Signature_sample_quality_summary.tsv")

    class_summary = (sample_results.groupby(["process", "Preop_platinum_exposure", "classification"], dropna=False).size()
        .rename("n").reset_index())
    class_summary["percent_within_process_exposure"] = class_summary["n"] / class_summary.groupby(["process", "Preop_platinum_exposure"])["n"].transform("sum") * 100
    write_tsv(class_summary, tables / "Input_Signature_classification_summary.tsv")

    test_rows: list[dict] = []
    descriptive_rows: list[dict] = []
    for process in PROCESSES:
        z = sample_results[(sample_results["process"] == process) & sample_results["Preop_platinum_exposure"].isin(["No", "Yes"])].copy()
        usable = z[z["classification"].isin(["Positive", "Not positive"])]
        pos = pd.crosstab(usable["Preop_platinum_exposure"], usable["classification"] == "Positive").reindex(index=["No", "Yes"], columns=[False, True], fill_value=0)
        no_neg, no_pos = int(pos.loc["No", False]), int(pos.loc["No", True])
        yes_neg, yes_pos = int(pos.loc["Yes", False]), int(pos.loc["Yes", True])
        if (no_pos + no_neg) > 0 and (yes_pos + yes_neg) > 0:
            _, p = fisher_exact([[yes_pos, yes_neg], [no_pos, no_neg]], alternative="two-sided")
            or_est, lo, hi = odds_ratio_ci(yes_pos, yes_neg, no_pos, no_neg)
        else:
            p = or_est = lo = hi = np.nan
        test_rows.append({"process": process, "metric": "Signature-positive proportion", "method": "Fisher exact",
            "no_n": no_pos + no_neg, "yes_n": yes_pos + yes_neg,
            "no_value": 100 * no_pos / (no_pos + no_neg) if (no_pos + no_neg) else np.nan,
            "yes_value": 100 * yes_pos / (yes_pos + yes_neg) if (yes_pos + yes_neg) else np.nan,
            "effect": or_est, "lower95": lo, "upper95": hi,
            "effect_definition": "Odds ratio for positive classification: Yes / No", "p": p})

        evalz = z[z["quality_evaluable"]]
        for metric, label in [("attributed_count", "Signature-attributed retained SNV count"), ("proportion", "Signature proportion")]:
            x_no = evalz.loc[evalz["Preop_platinum_exposure"] == "No", metric].dropna().to_numpy(float)
            x_yes = evalz.loc[evalz["Preop_platinum_exposure"] == "Yes", metric].dropna().to_numpy(float)
            if len(x_no) and len(x_yes):
                result = mannwhitneyu(x_yes, x_no, alternative="two-sided")
                p = float(result.pvalue); effect = rank_biserial(x_yes, x_no, float(result.statistic))
            else:
                p = effect = np.nan
            test_rows.append({"process": process, "metric": label, "method": "Wilcoxon rank-sum / Mann-Whitney U",
                "no_n": len(x_no), "yes_n": len(x_yes), "no_value": np.median(x_no) if len(x_no) else np.nan,
                "yes_value": np.median(x_yes) if len(x_yes) else np.nan, "effect": effect, "lower95": np.nan, "upper95": np.nan,
                "effect_definition": "Rank-biserial correlation; positive means higher in Yes", "p": p})
            for group, values in [("No", x_no), ("Yes", x_yes)]:
                descriptive_rows.append({"process": process, "metric": label, "group": group, "n": len(values),
                    "median": np.median(values) if len(values) else np.nan,
                    "q1": np.quantile(values, .25) if len(values) else np.nan,
                    "q3": np.quantile(values, .75) if len(values) else np.nan})
    tests = pd.DataFrame(test_rows)
    tests["p_holm_9_tests"] = holm(tests["p"])
    tests["p_display"] = tests["p"].map(lambda x: "NE" if not np.isfinite(x) else ("<0.001" if x < .001 else f"{x:.3f}"))
    tests["adjusted_p_display"] = tests["p_holm_9_tests"].map(lambda x: "NE" if not np.isfinite(x) else ("<0.001" if x < .001 else f"{x:.3f}"))
    tests["effect_display"] = tests["effect"].map(lambda x: "NE" if not np.isfinite(x) else f"{x:.2f}")
    tests["interpretation_scope"] = "Exploratory association; pre-treatment tissue prevents treatment-induced signature interpretation"
    write_tsv(tests, tables / "Input_Preop_platinum_signature_tests.tsv")
    write_tsv(pd.DataFrame(descriptive_rows), tables / "Input_Preop_platinum_signature_descriptive.tsv")

    metric_status = requested_metric_status(tests)
    write_tsv(metric_status, tables / "Input_Signature_requested_metrics_status.tsv")

    pooled, individual, assignment_audit = pooled_sensitivity_analysis(
        sbs96.loc[common], cohort, reference, private, args)
    write_tsv(pooled, tables / "Input_Preop_platinum_pooled_signature_sensitivity.tsv")
    write_tsv(individual, tables / "Input_Preop_platinum_pooled_signature_by_SBS.tsv")
    write_tsv(assignment_audit, tables / "Input_Signature_target_assignment_audit.tsv")

    thresholds = pd.DataFrame([
        ("Minimum retained SNVs", args.min_snv), ("Minimum reconstruction cosine", args.min_cosine),
        ("Minimum process-attributed retained SNVs", args.min_count), ("Minimum process proportion", args.min_proportion),
        ("Minimum bootstrap stability", args.min_stability), ("Bootstrap replicates", args.bootstraps),
        ("Multiple-testing correction", "Holm across 9 prespecified process-metric tests"),
        ("Pooled process sensitivity correction", "Holm across 6 process-metric permutation tests"),
        ("Pooled individual-SBS sensitivity correction", "Holm across 20 signature-metric permutation tests"),
        ("Platinum process", "SBS31 + SBS35"), ("HRD process", "SBS3"),
        ("MMR process", "SBS6 + SBS14 + SBS15 + SBS20 + SBS21 + SBS26 + SBS44"),
        ("Burden unit", "Signature-attributed retained SNV count; not mutations/Mb"),
        ("Specimen timing", "Preoperative treatment-naive per investigator; no treatment-induced interpretation"),
        ("Patient-level result", "Not estimable when no patient passes both retained-SNV and reconstruction thresholds"),
        ("Pooled sensitivity result", "Hypothesis-directed constrained refit; patient-bootstrap CI and label-permutation p value"),
    ], columns=["parameter", "value"])
    write_tsv(thresholds, tables / "Input_Signature_analysis_definitions.tsv")
    make_public_plots(sample_results, tests, pooled, individual, metric_status,
                      assignment_audit, figures, args)

    run_manifest = {
        "clinical_file": str(Path(args.clinical).resolve()), "clinical_sheet": args.clinical_sheet,
        "maf_file": str(Path(args.maf).resolve()), "reference_build": "GRCh37",
        "cosmic_version": "3.6", "signature_assignment": "SigProfilerAssignment 1.1.5",
        "matrix_generation": "SigProfilerMatrixGenerator 1.4.0", "targeted_panel_bed_available": False,
        "interpretation": "exploratory", "parameters": vars(args),
    }
    (private / "run_manifest.json").write_text(json.dumps(run_manifest, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
