#!/usr/bin/env python3
"""Read-only contract checks for the current five-figure manuscript package."""
import argparse
import csv
import json
import math
from pathlib import Path
from zipfile import ZipFile


ap = argparse.ArgumentParser()
ap.add_argument("--results", type=Path, required=True)
args = ap.parse_args()
root = args.results.resolve()


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table(name):
    return read_tsv(root / "tables" / name)


prov = json.loads((root / "private" / "provenance.json").read_text(encoding="utf-8"))
assert prov["kras_source"] == "workbook_maf"
assert prov["kras_source_column"] == "KRAS_subtype_MAF"
assert prov["source_comparisons_included"] is False
assert prov["main_figures"] == 5 and prov["main_tables"] == 2
assert any(Path(path).name == "260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx" for path in prov["input_md5"])

flow = {row["item"]: int(float(row["n"])) for row in table("Input_Cohort_flow.tsv")}
n = flow["Analysis patients"]
assert flow["Source records"] - flow["Duplicate-ID records excluded"] == n == 1011
assert flow["M0"] + flow["M1"] == n
assert flow["TMB observed"] + flow["TMB missing"] == n

dataset = read_tsv(root / "private" / "analysis_dataset.tsv")
assert len(dataset) == n
assert len({row["patient_id"] for row in dataset}) == n
assert all(row["patient_duplicate"] == "FALSE" for row in dataset)
workbook_to_analysis = {
    "Not detected": "Not detected", "G12D": "G12D", "G12V": "G12V",
    "G12R": "G12R", "Other": "Other",
}
assert all(workbook_to_analysis[row["KRAS_subtype_MAF"]] == row["KRAS_group"] for row in dataset)

audit = read_tsv(root / "private" / "variant_count_audit.tsv")[0]
patient_counts = read_tsv(root / "private" / "Figure1B_Oncoplot_patient_counts.tsv")
gene_counts = table("Input_Figure1B_Oncoplot_gene_counts.tsv")
assert any(row["gene"] == "KRAS" for row in gene_counts)
assert len(patient_counts) == n
assert sum(int(row["variant_count"]) for row in patient_counts) == int(audit["unique_patient_events"])

track_contracts = {
    "Figure1B_Oncoplot": ["KRAS_group", "Sex", "Age_group", "Differentiation", "Neoadjuvant", "T_stage", "N_stage", "LN_positive", "M_stage", "BMI", "CA19_9", "CEA", "AJCC_stage", "Stage_Group", "Tumor_size", "TMB"],
    "Supplementary_Figure1_HRD": ["KRAS_group", "Age", "Sex", "Neoadjuvant", "Preop_platinum_exposure", "T_stage", "N_stage", "M_stage", "Differentiation", "R_status", "MSI", "TMB"],
    "Supplementary_Figure2_MMR": ["KRAS_group", "Age", "Sex", "Neoadjuvant", "Preop_platinum_exposure", "T_stage", "N_stage", "M_stage", "Differentiation", "R_status", "MSI", "TMB"],
}
for figure, expected in track_contracts.items():
    tracks = table(f"Input_{figure}_clinical_tracks.tsv")
    assert [row["column"] for row in tracks] == expected
    assert all(int(row["observed"]) + int(row["missing"]) == n for row in tracks)
    assert all(row["color_mapping"] and row["missing_color"] == "#BDBDBD" for row in tracks)
    assert len({row["color_mapping"] for row in tracks}) == len(tracks)

assert not list((root / "tables").glob("Supplementary_Table*"))
required_tables = [
    "Input_Preop_platinum_population.tsv", "Input_Preop_platinum_global_comparisons.tsv",
    "Input_Preop_platinum_pairwise_comparisons.tsv", "Input_Preop_platinum_Cox_models.tsv",
    "Input_Recurrence_pattern_distribution.tsv", "Input_Recurrence_pattern_global_comparisons.tsv",
    "Input_Recurrence_pattern_pairwise_comparisons.tsv", "Input_Recurrence_pattern_post_OS_pairwise_logrank.tsv",
    "Input_Recurrence_pattern_Cox_models.tsv", "Table1_Cohort_characteristics.tsv",
    "Table2_Main_survival_models.tsv", "Input_Distant_only_pattern_post_OS_logrank.tsv",
    "Input_Distant_only_pattern_post_OS_Cox_models.tsv", "Input_Expanded_recurrence_post_OS_logrank.tsv",
    "Input_Expanded_recurrence_post_OS_Cox_models.tsv", "Input_Preop_platinum_HRD_MMR_set_tests.tsv",
    "Input_Preop_platinum_HRD_MMR_gene_tests.tsv", "Input_MAF_signature_audit.tsv",
    "Input_Signature_sample_quality_summary.tsv", "Input_Preop_platinum_pooled_signature_sensitivity.tsv",
    "Input_Preop_platinum_pooled_signature_by_SBS.tsv", "Input_Signature_requested_metrics_status.tsv",
    "Input_Signature_target_assignment_audit.tsv",
]
assert all((root / "tables" / name).is_file() for name in required_tables)

signature_status = table("Input_Signature_requested_metrics_status.tsv")
assert len(signature_status) == 9
assert all(row["status"] == "Not estimable" and row["result"] == "NE" for row in signature_status)
individual_sbs = table("Input_Preop_platinum_pooled_signature_by_SBS.tsv")
assert len(individual_sbs) == 20
assert {row["signature"] for row in individual_sbs} == {
    "SBS31", "SBS35", "SBS3", "SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44"
}
assert all(float(row["p_holm_20_tests"]) >= .05 for row in individual_sbs)
assignment_audit = table("Input_Signature_target_assignment_audit.tsv")
assert len(assignment_audit) == 10
assert all(float(row["no_standard_assignment_count"]) == 0 and
           float(row["yes_standard_assignment_count"]) == 0 for row in assignment_audit)

platinum = table("Input_Preop_platinum_Cox_models.tsv")
assert len(platinum) == 4
assert all(math.isfinite(float(row["HR"])) and row["direction"] for row in platinum)
assert float(platinum[2]["PH_global_p"]) < 0.05
assert float(platinum[3]["PH_global_p"]) > 0.05

recurrence = table("Input_Recurrence_pattern_distribution.tsv")
assert sum(int(row["n"]) for row in recurrence) == 659
assert {row["Recurrence_group"] for row in recurrence} == {"Distant only", "Local only", "Local + distant"}
main_survival = table("Table2_Main_survival_models.tsv")
assert len(main_survival) == 21
assert all(row["comparison_group"] and row["reference_group"] and row["direction"] for row in main_survival)

required_figures = [
    "Figure1A_Cohort.png", "Figure1B_Oncoplot.png", "Figure1C_KRAS_distribution.png",
    "Figure2A_Pathology_matrix.png", "Figure2_TP53_differentiation.png", "Figure2_CDKN2A_N_stage.png", "Figure2_CDKN2A_N_category.png",
    "Figure3_Platinum_context.png", "Figure3_Platinum_OS.png", "Figure3_Platinum_forest.png",
    "Figure4_Recurrence_distribution.png", "Figure4_Recurrence_context.png", "Figure4_Post_recurrence_OS.png", "Figure4_Post_recurrence_forest.png",
    "Figure3A_Gene_OS.png", "Figure4A_KRAS_OS.png", "Figure4B_KRAS_pairwise.png",
    "Table1_Cohort_characteristics.png", "Table2_Main_survival_models.png",
    "Supplementary_Figure1_HRD.png", "Supplementary_Figure2_MMR.png", "Supplementary_Figure4_Reported_TMB.png",
    "Figure3_Platinum_HRD_MMR.png", "Figure4_Distant_only_pattern_KM.png", "Figure4_Distant_only_pattern_forest.png",
    "Figure4_Expanded_recurrence_KM.png", "Figure4_Expanded_recurrence_forest.png",
    "Supplementary_Figure_Signature_QC.png", "Supplementary_Figure_Platinum_HRD_MMR_signatures.png",
    "Supplementary_Figure_Individual_SBS_Pooled_Sensitivity.png",
    "Supplementary_Figure_Signature_Assignment_Audit.png",
    "Supplementary_Figure_Signature_Requested_Metrics.png",
]
assert all((root / "figures" / name).is_file() for name in required_figures)

docx_candidates = sorted((root / "deliverables").glob("PAAD_Updated_Manuscript_Report_*_KR_Explained_v*.docx"))
pptx_candidates = sorted((root / "deliverables").glob("PAAD_Updated_Manuscript_Presentation_*_KR_Explained_v*.pptx"))
assert docx_candidates and pptx_candidates
docx = docx_candidates[-1]
pptx = pptx_candidates[-1]
with ZipFile(docx) as archive:
    styles = archive.read("word/styles.xml").decode()
    body = archive.read("word/document.xml").decode()
    assert "Pretendard" in styles
    assert "260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx" in body
    assert "KRAS_subtype_MAF" in body
    assert "Supplementary Table 1" not in body
with ZipFile(pptx) as archive:
    slides = [name for name in archive.namelist() if name.startswith("ppt/slides/slide") and name.endswith(".xml")]
    joined = "".join(archive.read(name).decode(errors="ignore") for name in slides)
    assert len(slides) == 44
    assert "260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx" in joined
    assert "KRAS_subtype_MAF" in joined
    assert "Pretendard" in joined

print(
    f"PASS: specified workbook; {n} unique patients; 5 main figures / 2 main tables; "
    "platinum and recurrence extensions; no Supplementary Table series; DOCX/PPTX package checks."
)
