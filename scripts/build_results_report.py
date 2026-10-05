#!/usr/bin/env python3
"""Build an editable manuscript-planning DOCX from completed, aggregate R outputs.
Run AFTER run_all.R. Requires python-docx, pandas, Pillow; never reads raw data.
The Word report is a draft for scientific review, not a submission-ready claim.
"""
import argparse
from pathlib import Path
from io import BytesIO
import csv
import math
import re
import hashlib
import pandas as pd
from report_interpretations import expanded_interpretation
from PIL import Image
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ap = argparse.ArgumentParser()
ap.add_argument("--results", type=Path, required=True)
ap.add_argument("--output", type=Path)
ap.add_argument("--legacy-internal-qc", action="store_true", help="Explicitly allow the historical internal-QC report; not a manuscript report.")
args = ap.parse_args()
if not args.legacy_internal_qc:
    raise SystemExit("This historical report includes internal source comparisons. For the current manuscript use scripts/build_manuscript_documents.py --results outputs/.../manuscript. Use --legacy-internal-qc only to reproduce the old internal report.")
out = args.results.resolve()
target = args.output or out / "PAAD_Manuscript_Results_Report_20260928.docx"
def read(name):
    return pd.read_csv(out / name, sep="\t")
def fp(x):
    if pd.isna(x) or not math.isfinite(float(x)): return "NE"
    return "<0.001" if x < .001 else f"{x:.3f}"
def hr(r):
    return f"{r.HR:.2f} ({r.lower95:.2f}–{r.upper95:.2f})"
def yes(x): return str(x).lower() == "true"
def clean(x):
    if pd.isna(x): return "NE"
    return str(x)
def display(x):
    replacements = {"MAF_variant_count":"Variant count", "KRAS_clinical_group":"Clinical KRAS",
      "KRAS_MAF_group":"MAF KRAS", "KRAS_MAF":"MAF KRAS",
      "Age10":"Age (per 10 years)", "Neoadjuvant":"Neoadjuvant", "T_stage":"T category",
      "N_stage":"N category", "M_stage":"M category", "R_status":"Margin",
      "Operation":"Surgery type", "Adjuvant_treatment":"Adjuvant treatment",
      "Preop_platinum_exposure":"Preoperative platinum", "driver_group":"Driver gene count"}
    return replacements.get(x, x)

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Inches(8.5), Inches(11)
sec.top_margin = sec.bottom_margin = Inches(.65)
sec.left_margin = sec.right_margin = Inches(.75)
sec.header_distance = sec.footer_distance = Inches(.28)
for name in ["Normal", "Body Text", "Caption", "Title", "Subtitle", "Heading 1", "Heading 2"]:
    st = doc.styles[name]
    st.font.name = "Pretendard"
    st.font.color.rgb = RGBColor(0, 0, 0)
    st.element.get_or_add_rPr().get_or_add_rFonts().set(qn("w:eastAsia"), "Pretendard")
doc.styles["Normal"].font.size = Pt(11)
# Prevent theme fonts and the default Word Title border from overriding Pretendard.
for st in doc.styles:
    if st.type in (1, 2):
        st.font.name = "Pretendard"
        fonts = st.element.get_or_add_rPr().get_or_add_rFonts()
        for attr in list(fonts.attrib):
            if "Theme" in attr or "theme" in attr: del fonts.attrib[attr]
        for key in ("ascii", "hAnsi", "eastAsia", "cs"): fonts.set(qn("w:"+key), "Pretendard")
        for border in list(st.element.iter(qn("w:pBdr"))): border.getparent().remove(border)
doc.styles["Normal"].paragraph_format.space_after = Pt(6)
doc.styles["Normal"].paragraph_format.line_spacing = 1.12
for name, size in [("Title",24),("Heading 1",16),("Heading 2",12)]:
    doc.styles[name].font.size = Pt(size)
    doc.styles[name].font.bold = True
    doc.styles[name].paragraph_format.space_after = Pt(8)
head = sec.header.paragraphs[0]
head.text = "PAAD  |  UPDATED CLINICAL & GENOMIC ANALYSIS"
head.runs[0].font.size = Pt(8)
foot = sec.footer.paragraphs[0]
foot.alignment = 2
foot.add_run("PAAD 분석 결과  ·  ")
fld = OxmlElement("w:fldSimple"); fld.set(qn("w:instr"), "PAGE")
foot._p.append(fld)
doc.core_properties.title = "PAAD 임상·유전체 분석 결과 및 원고용 Figure·Table"
doc.core_properties.subject = "All-pairs inference, source comparison, manuscript tables"
doc.core_properties.author = "PAAD analysis"

def para(text, size=None, bold=False):
    p = doc.add_paragraph()
    r = p.add_run(text); r.bold = bold
    if size: r.font.size = Pt(size)
    return p
def newpage(title):
    doc.add_page_break()
    doc.add_heading(title, 1)
def table(headers, rows, widths=None, size=8.5):
    t = doc.add_table(rows=1, cols=len(headers)); t.autofit = False
    widths = widths or [7/len(headers)]*len(headers)
    for col,w in zip(t.columns,widths): col.width = Inches(w)
    for cell,w in zip(t.rows[0].cells,widths): cell.width = Inches(w)
    for c,s in zip(t.rows[0].cells,headers): c.text = str(s)
    h = OxmlElement("w:tblHeader"); t.rows[0]._tr.get_or_add_trPr().append(h)
    for row in rows:
        cells = t.add_row().cells
        for c,w,v in zip(cells,widths,row): c.width = Inches(w); c.text = clean(v)
    pr = t._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge in ["top","left","bottom","right","insideH","insideV"]:
        el=OxmlElement("w:"+edge); el.set(qn("w:val"),"single")
        el.set(qn("w:sz"),"4"); el.set(qn("w:color"),"D9D9D9"); borders.append(el)
    pr.append(borders)
    for ri,row in enumerate(t.rows):
        cant=OxmlElement("w:cantSplit"); row._tr.get_or_add_trPr().append(cant)
        for cell in row.cells:
            tcpr=cell._tc.get_or_add_tcPr()
            margins=OxmlElement("w:tcMar")
            for e,v in [("top",65),("bottom",65),("left",70),("right",70)]:
                x=OxmlElement("w:"+e); x.set(qn("w:w"),str(v)); x.set(qn("w:type"),"dxa"); margins.append(x)
            tcpr.append(margins)
            if ri==0:
                sh=OxmlElement("w:shd"); sh.set(qn("w:fill"),"E9EDF0"); tcpr.append(sh)
            for p in cell.paragraphs:
                p.paragraph_format.space_after=Pt(0); p.paragraph_format.line_spacing=1.05
                for r in p.runs: r.font.size=Pt(size); r.bold=ri==0
    doc.add_paragraph().paragraph_format.space_after = Pt(0)
    return t
supp = []
def refs(group, files):
    for name in files:
        existing=next((r for r in supp if r["file"]==name),None)
        if existing: label=existing["label"]
        else:
            number=sum(r["group"]==group for r in supp)+1
            label=f"Supplementary table{group}-{number}."
            supp.append(dict(group=group,label=label,file=name))
        assert (out/name).exists(), name
        assert "_private" not in name, "Patient-level data must not enter publication attachments"
        p=para(f"{label} {name}",8)
        p.paragraph_format.space_after=Pt(3)
def export_table(name, headers, rows):
    pd.DataFrame(rows,columns=headers).to_csv(out/name,sep="\t",index=False)
def picture(name, maxheight=4.1):
    im=Image.open(out/name); im.thumbnail((2400,2400),Image.Resampling.LANCZOS)
    stream=BytesIO(); im.save(stream,format="PNG"); stream.seek(0)
    w=min(7.0,maxheight*im.width/im.height)
    doc.add_picture(stream,width=Inches(w))
    doc.paragraphs[-1].alignment=1
def pair_text(variable,endpoint=None,file="26_KRAS_pairwise_logrank.tsv"):
    p=read(file); p=p[p.variable==variable]
    if endpoint: p=p[p.endpoint==endpoint]
    good=p[p.p_holm<.05]
    if good.empty: return f"총 {len(p)}개 쌍의 Holm 보정 후 유의한 비교는 없었다."
    parts=[f"{r.group1} vs {r.group2} (보정 p={fp(r.p_holm)})" for r in good.itertuples()]
    return f"전체 {len(p)}개 쌍 중 보정 후 유의한 비교: "+"; ".join(parts)+"."
def figure(num,title,file,method,result,files,limit=""):
    newpage(f"Fig{num}. {title}")
    picture(file, maxheight=5.2 if num in (21,22) else 4.1)
    para("분석  "+method,9.5)
    para("해석  "+result,10.5)
    extra = expanded_interpretation(num, read)
    if extra: para(extra,10.5)
    if limit: para("주의  "+limit,9)
    para("Figure 원본: "+file,8)
    refs(num,files)

flow=read("18_Cohort_flow.tsv")
cox=read("20_OS_Cox_coefficients.tsv")
sc=read("29_Source_specific_Cox_coefficients.tsv")
cmh=read("19_Neoadjuvant_adjusted_pathology_CMH.tsv")
agreement=read("25_KRAS_agreement_summary.tsv")
corr=read("27_TMB_MAF_rank_correlation.tsv").iloc[0]
burden=read("27_TMB_distribution_summary.tsv")
kmt=read("26_KRAS_logrank_tests.tsv")
genes=["KRAS","TP53","SMAD4","CDKN2A"]
n=int(flow.loc[flow.step=="Primary M0+M1 unique patients","n"].iloc[0])
rec_n=int(flow.loc[flow.step=="Valid recorded recurrence follow-up","n"].iloc[0])
doc.add_heading("PAAD 임상·유전체 분석 결과",0)
para("논문용 Figure·Table 구성안 및 상세 결과 해석",14,True)
para("2026-09-28  |  updated_v3_20260928",10)
para("연구 질문: 주요 driver 유전자 및 KRAS 아형은 병리 특성과 생존에 어떤 관련이 있으며, 임상 보고와 MAF에서 얻은 지표는 얼마나 일관적인가?")
doc.add_heading("핵심 결과",1)
para(f"주 분석은 중복 ID에 해당하는 4행을 제외한 {n:,}명이다. M0 972명과 M1 39명을 포함하며, OS 사망은 631건이다. 기록된 재발 분석은 {rec_n:,}명, 재발 659건이다. M1을 전체 분석에서 제외하지 않았다.")
a=cox[(cox.family=="Adjusted genomic") & cox.is_exposure.map(yes) & cox.exposure.isin(genes)]
para("주요 유전자별 보정 OS 결과: "+"; ".join(f"{r.exposure} HR {hr(r)}, p={fp(r.p)}, BH q={fp(r.q_BH)}" for r in a.itertuples())+".")
para(f"KRAS 축약 범주는 {int(agreement.iloc[0].agreement_n)}/{int(agreement.iloc[0].n)}명에서 일치했다. 정확한 변이형을 판정할 수 있는 쌍에서는 {int(agreement.iloc[1].agreement_n)}/{int(agreement.iloc[1].n)}쌍이 일치했다. 두 분모를 혼용하면 안 된다.")
para(f"임상 TMB는 836명에서 관측되고 175명에서 결측이다. Variant count는 1,011명 모두 관측되며, 동일 환자 836쌍의 Spearman ρ={corr.spearman_rho:.3f} (95% CI {corr.lower95:.3f}–{corr.upper95:.3f})이다. 높은 순위 상관만으로 서로 대체 가능한 TMB라고 판단할 수 없다.")
doc.add_heading("읽는 순서",1)
para("Fig1–3: 코호트·유전체·병리 → Fig4–5: 주요 유전자와 OS → Fig6–12: KRAS 출처별 분석과 일치도 → Fig13–17: TMB·Variant count 및 임상 연관성 → Fig18–20: 재발 양상과 DNA 복구 유전자 → Fig21–25: HRD 관련·MMR oncoplot과 임상 분포.")
para("Table1은 임상 특성, Table2는 유전자별 병리 특성, Table3은 임상 OS 모형, Table4는 유전체 OS 모형, Table5는 두 수치 지표의 관측·분포를, Table6–7은 HRD 관련·MMR 유전자의 빈도와 임상 특성을 제시한다. 표는 Word에서 편집 가능하다.")
para("논문 제출 전 단위, 표본 추출·시퀀싱 시점, 검체 범위, 치료 시점 및 endpoint 정의를 연구자가 최종 확인해야 한다.",9)

newpage("분석 원칙과 참고논문 적용")
table(["항목","이번 분석의 적용"],[
["대상","새 임상정보 기준, 중복 ID 4행 제외. OS 주 분석에 M0와 M1 포함."],
["전체 검정","KM: log-rank; 연속형: Kruskal–Wallis/2군 Wilcoxon; 범주형: Fisher; 병리: neoadjuvant 층화 CMH."],
["쌍별 검정","모든 그룹 조합을 저장. log-rank, Dunn(동순위 보정), Fisher, 층화 CMH, Cox Wald contrast. 전체 p와 무관하게 계산."],
["다중비교","각 변수×그룹변수 또는 모형×endpoint 내 전체 쌍에 Holm 보정. 기존 omnibus family의 BH q와 다르며 연구 전체 오류율을 통제하는 것은 아님."],
["표시","KM에는 global log-rank p와 Holm p<0.05인 쌍. 전 조합의 n, 사건 수(생존), 통계량, 원 p, 보정 p, 계산 불가 사유는 TSV."],
["Cox","기본: 나이/10, 성별, T, Neoadjuvant, M. 완전사례 및 endpoint 내 변이가 없는 항목 제외. 확장·M0-only 민감도는 별도."],
["신뢰구간","Cox·RMST 95% CI는 개별(pointwise) 구간이며 동시 신뢰구간이 아님. 유의성 판단에는 해당 보정 p/q 확인."],
["재발","기록된 재발 사건과 재발 관찰기간. 사망 포함 DFS/RFS나 경쟁위험 누적발생률 아님."],
["출처 구분","Clinical KRAS와 MAF KRAS 별도 분석. No KRAS call은 검증된 biological WT가 아님."],
["수치 지표","Clinical TMB=mutations/Mb; Variant count=유지된 MAF 행 수. 코드 내부명 MAF_variant_count는 추적성을 위해 유지."],
], [1.2,5.8],9)
para("참고논문 적용: McIntyre의 임상 특성·driver/아형별 OS·재발 및 DNA 복구 유전자 구성을 참고했다. Campbell의 neoadjuvant 보정 병리 비교와 임상/분자 OS 표 구조를 반영했다. 원 논문 수치·표를 복제하지 않고 현재 코호트의 집계값으로 작성했다.",9.5)
para("차이와 제한: 두 논문의 절제·국소성 중심 코호트와 달리 M1을 포함한다. LOH/CNV, tumor purity 및 tumor cellularity 관련 분석은 이번 분석 범위에서 제외했다. 생식세포 변이, 병원성 판정, Clavien–Dindo 및 재발 전 사망의 경쟁위험 정의는 확보되지 않았다. 보조치료 시점·시퀀싱 진입 시점이 없어 치료 인과효과나 지연 진입 보정은 주장하지 않는다.",9.5)

figure(1,"분석 코호트와 endpoint별 자료 가용성","18_Cohort_flow.png",
 "새 임상자료의 표본 수, 중복 제외, M1 포함 및 재발 추적 가능 인원을 집계했다.",
 "1,015행에서 중복 ID 4행을 제외한 1,011명을 사용한다. 재발 분석의 969명은 전 코호트와 구분한다.",
 ["18_Cohort_flow.tsv","00_Cohort_audit_counts.tsv"],"각 막대는 서로 더하는 구성요소가 아니다.")

# Table1: adapt demographic/therapy/pathology sections, with explicit missingness.
base=read("18_Extended_clinical_Table1.tsv")
cohorts=["Primary M0+M1 unique patients","Primary: no neoadjuvant","Primary: neoadjuvant"]
groups=[
("A. 인구학적·치료 특성",["Age","Sex","BMI","Diabetes","ASA","CA19_9","CEA","Operation","Preop_platinum_exposure","Adjuvant_treatment"]),
("B. 병리 특성",["Tumor_size","T_stage","N_stage","M_stage","Differentiation","LVI","PNI","R_status"]),
]
allbaseline=[]
for panel,variables in groups:
    newpage("Table1. 임상·병리 특성 — "+panel)
    headers=["Characteristic","Overall","No neoadjuvant","Neoadjuvant"]
    rows=[]
    for var in variables:
        b=base[(base.variable==var)&(base.cohort==cohorts[0])]
        for r in b.itertuples():
            category=r.category
            vals=[]
            for c in cohorts:
                q=base[(base.variable==var)&(base.category==category)&(base.cohort==c)]
                if q.empty: vals.append("0 (0.0%)"); continue
                rr=q.iloc[0]
                vals.append(f"{rr['median']:.1f} [{rr.q1:.1f}, {rr.q3:.1f}]; m={int(rr.missing_n)}"
                  if category=="Continuous" else f"{int(rr.n)} ({rr.percent:.1f}%)")
            label=display(var)+(f": {category}" if category!="Continuous" else "")
            rows.append([label,*vals])
    table(headers,rows,[2.25,1.58,1.58,1.59],8)
    allbaseline.extend(rows)
    para("연속형: median [IQR]; m=결측 수. 범주형: n(%), 각 열의 전체 코호트 분모(결측 포함). 전체 1,011명, No neoadjuvant 664명, Neoadjuvant 347명. 미측정은 No로 합치지 않았다. 종양 크기·혈청 표지자 단위는 입력 사전 확인 후 원고에 확정 기재할 것.",8.5)
    refs(1,["18_Extended_clinical_Table1.tsv","18_Clinical_Table1_display.tsv"])
export_table("Table1_Manuscript_baseline.tsv",headers,allbaseline)

figure(2,"주요 유전자 변이와 임상 특성의 개관","1_Oncoplot_Main_PDAC.png",
 "유지된 MAF 변이와 임상 주석을 환자별로 정렬한 oncoplot이다.",
 "KRAS·TP53 등 빈번한 driver를 출발점으로 이후 병리 및 생존 분석을 구성한다. 색은 변이 분류이며 빈도 자체가 예후 효과는 아니다.",
 ["28_Categorical_clinical_summary.tsv","28_Continuous_clinical_summary.tsv"],
 "열은 검체/환자이며 그림의 no-call을 검증된 WT로 단정하지 않는다. 상세 원본 PNG를 함께 확인할 것.")
figure(3,"Neoadjuvant 보정 후 driver와 병리 특성의 연관성","19_Adjusted_pathology_associations.png",
 "4개 유전자×6개 병리 특성의 CMH 검정, 24개 전체 검정에 BH 보정. 범주 간 쌍별 CMH는 Holm 보정.",
 "BH q<0.05인 조합은 "+("; ".join(f"{r.gene}–{r.characteristic} (q={fp(r.q_BH)})" for r in cmh[cmh.q_BH<.05].itertuples()) or "없었다")+".",
 ["19_Neoadjuvant_adjusted_pathology_CMH.tsv","19_Pathology_by_driver_counts.tsv","19_Pathology_pairwise_CMH.tsv","13_Driver_pathology_pairwise_Fisher.tsv"],
 "조정은 neoadjuvant 층화에 한정되며 인과관계나 모든 교란 통제를 의미하지 않는다.")
counts=read("19_Pathology_by_driver_counts.tsv")
table2=[]
for panel,gg in [("A. KRAS / TP53",genes[:2]),("B. SMAD4 / CDKN2A",genes[2:])]:
    newpage("Table2. Driver별 병리 특성 — "+panel)
    heads=["Characteristic"]+[f"{g}\nNo call / Variant\nn (%)" for g in gg]+[f"{g}\np / BH q" for g in gg]
    rows=[]
    for var in ["T_stage","N_stage","Differentiation","LVI","PNI","R_status"]:
        cats=sorted(counts[counts.characteristic==var].category.astype(str).unique())
        for j,cat in enumerate(cats):
            vals=[]; ps=[]
            for g in gg:
                cc=counts[(counts.gene==g)&(counts.characteristic==var)&(counts.category.astype(str)==cat)]
                vals.append(" / ".join(f"{int(z.iloc[0].n)} ({z.iloc[0].percent:.1f})" if not z.empty else "0 (0.0)" for status in ["No retained call","Retained variant"] for z in [cc[cc.gene_status==status]]))
                rr=cmh[(cmh.gene==g)&(cmh.characteristic==var)].iloc[0]
                ps.append(f"{fp(rr.p)} / {fp(rr.q_BH)}" if j==0 else "—")
            rows.append([display(var)+": "+cat,*vals,*ps])
    table(heads,rows,[1.75,1.65,1.65,.975,.975],8)
    table2.extend([[panel,*r] for r in rows])
    para("No call / Variant 순서. 비율은 해당 특성에 값이 있는 유전자 상태군 내 분모이다. p/q는 개별 행이 아니라 해당 특성 전체 분포에 대한 CMH 검정이며 neoadjuvant로 층화했다. 모든 범주 쌍의 보정 p는 별도 표에 제공한다.",8.5)
    refs(3,["19_Pathology_by_driver_counts.tsv","19_Neoadjuvant_adjusted_pathology_CMH.tsv","19_Pathology_pairwise_CMH.tsv"])
export_table("Table2_Manuscript_pathology.tsv",["Panel","Characteristic","Gene1 no call / variant","Gene2 no call / variant","Gene1 p / q","Gene2 p / q"],table2)

figure(4,"주요 유전체 지표와 전체생존의 보정 연관성","20_Adjusted_OS_forest.png",
 "각 유전체 지표를 별도 Cox 모형에 넣고 나이, 성별, T, neoadjuvant, M을 보정했다.",
 "점은 HR, 선은 95% CI이다. HR>1은 비교군에서 높은 사망 hazard와의 연관성을 뜻한다. 다중 범주의 모든 쌍은 별도 Wald contrast 표에 있다.",
 ["20_OS_Cox_coefficients.tsv","20_OS_Cox_PH_tests.tsv","20_OS_Cox_model_diagnostics.tsv","20_Cox_pairwise_contrasts.tsv","21_RMST_pairwise_contrasts.tsv"],
 "별표가 붙은 항목은 비례위험 가정 위반 신호가 있어 단일 HR을 시간 전체에 일정한 효과로 해석하면 안 된다.")
# Table3: align clinical univariable terms with joint model.
newpage("Table3. 임상 특성과 OS — 단변량 및 다변량 Cox")
uni=cox[(cox.family=="Univariable clinical")&cox.is_exposure.map(yes)]
joint=cox[cox.family=="Joint clinical model"]
rows=[]
for r in uni.itertuples():
    match=joint[joint.term==r.term]
    rr=match.iloc[0] if not match.empty else None
    rows.append([r.term,f"{int(r.n)}/{int(r.events)}",hr(r),fp(r.p),
      hr(rr) if rr is not None else "Not in model",fp(rr.p) if rr is not None else "—"])
heads=["Term (vs reference)","Univ.\nN/deaths","Univariable\nHR (95% CI)","p","Adjusted\nHR (95% CI)","p"]
table(heads,rows,[1.65,.6,1.65,.55,1.85,.7],8)
para("기준: Sex F, T 1, N 0, Neoadjuvant n, LVI/PNI Negative, Margin R0, Differentiation WD, OP CompTP, Diabetes n. Age는 10세 증가당. 실제 term 표기를 보존했으며 범주별 표본 수는 Table1/Supplementary에 제공한다.",8.5)
para("N/deaths 열은 단변량 분석의 표본 수이다. 다변량은 1,000명/사망 623건으로 나이, 성별, T, N, neoadjuvant, M, LVI, PNI, margin, differentiation의 공동 모형이다. 이 표에 없는 공변량을 포함한 전체 계수와 PH 진단은 원자료 참조. Not in model은 검정 결과가 유의하지 않다는 뜻이 아니다.",8.5)
refs(4,["20_OS_Cox_coefficients.tsv","20_OS_Cox_model_diagnostics.tsv","20_Cox_pairwise_contrasts.tsv"])
export_table("Table3_Manuscript_clinical_OS.tsv",heads,rows)

newpage("Table4. 유전체 지표와 OS — 단변량 및 보정 Cox")
rows=[]
for r in cox[(cox.family=="Adjusted genomic")&cox.is_exposure.map(yes)].itertuples():
    u=cox[(cox.family=="Univariable genomic")&(cox.exposure==r.exposure)&(cox.term==r.term)]
    label=r.term.replace("KRAS_MAF","KRAS ").replace("driver_group","Driver genes ").replace("KRAS_TP53","KRAS/TP53 ")
    rows.append([label,f"{int(r.n)}/{int(r.events)}",hr(u.iloc[0]) if not u.empty else "NE",hr(r),fp(r.p),fp(r.q_BH),fp(r.PH_exposure_p)])
heads=["Comparison term","N/deaths","Univariable\nHR (95% CI)","Adjusted\nHR (95% CI)","p","BH q","PH p"]
table(heads,rows,[1.55,.6,1.55,1.55,.55,.55,.65],7.8)
para("기준: 각 driver/Any_driver는 no retained call(0); KRAS subtype은 No KRAS call; driver_group은 1개; TP53 truncation은 Nontruncating; KRAS/TP53는 Neither; G12D는 Other KRAS. HR은 표의 수준/기준 수준이다. 나이/10, 성별, T, neoadjuvant, M 보정. PH p는 노출변수의 비례위험 검정이다.",8.5)
refs(4,["20_OS_Cox_coefficients.tsv","20_OS_Cox_PH_tests.tsv","20_Cox_pairwise_contrasts.tsv"])
export_table("Table4_Manuscript_genomic_OS.tsv",heads,rows)

figure(5,"Driver 유전자 수와 전체생존","21_OS_driver_count_0_1_2_3plus.png",
 "4개 driver 중 retained call이 있는 유전자 수를 0/1/2/3+로 분류한 KM 곡선이다.",
 pair_text("driver_group",file="21_OS_pairwise_logrank.tsv"),
 ["21_OS_logrank_tests.tsv","21_OS_pairwise_logrank.tsv","21_OS_medians_and_landmark_estimates.tsv","21_Unadjusted_RMST_36months.tsv"],
 "유전자 수와 Variant count는 서로 다른 변수이다. 드문 군의 표본 수와 위험집단 감소를 함께 확인한다.")
figure(6,"임상 보고와 MAF의 KRAS 아형 분포","26_KRAS_source_distributions.png",
 "동일한 1,011명의 KRAS를 출처별로 No call/G12D/G12V/G12R/Other로 요약했다.",
 "전체 분포가 비슷해도 개별 환자의 분류가 일치한다는 뜻은 아니다. 세부 아형 빈도 및 교차표를 함께 확인한다.",
 ["26_KRAS_group_frequencies.tsv","26_KRAS_subtype_detailed_frequencies.tsv","26_KRAS_subtype_raw_detailed_frequencies.tsv"],
 "Other에는 희귀·다중 변이가 포함되며 임상 원문과 MAF raw 문자열을 별도로 보존했다.")
for num,var,endpoint,title in [
(7,"KRAS_clinical_group","OS","임상 KRAS 아형과 전체생존"),
(8,"KRAS_MAF_group","OS","MAF KRAS 아형과 전체생존"),
(9,"KRAS_clinical_group","Recorded_recurrence","임상 KRAS 아형과 기록된 재발"),
(10,"KRAS_MAF_group","Recorded_recurrence","MAF KRAS 아형과 기록된 재발")]:
    r=kmt[(kmt.variable==var)&(kmt.endpoint==endpoint)].iloc[0]
    figure(num,title,f"26_{var}_{endpoint}.png",
      f"Kaplan–Meier, N={int(r.n)}, 사건 {int(r.events)}건; global log-rank p={fp(r.p)}, endpoint 내 BH q={fp(r.q_BH)}.",
      pair_text(var,endpoint),
      ["26_KRAS_pairwise_logrank.tsv","26_KRAS_KM_medians.tsv","26_KRAS_KM_number_at_risk.tsv","26_KRAS_logrank_tests.tsv"],
      "곡선 차이는 미보정 연관성이다." if endpoint=="OS" else "사망을 사건으로 추가하지 않은 기록된 재발 endpoint이다. 사망 포함 DFS/RFS 또는 CIF로 명명하지 않는다.")
figure(11,"KRAS 출처별 보정 예후 모형","29_KRAS_source_adjusted_outcomes.png",
 "두 출처를 각각 사용한 OS/기록된 재발 Cox 모형. 아형 기준은 No call, G12D 비교 기준은 Other KRAS이다.",
 "출처별 모형의 HR 방향·정밀도를 비교한다. 한 출처의 p<0.05와 다른 출처의 p≥0.05만으로 출처 간 효과 차이를 입증하지 않는다.",
 ["29_Source_specific_Cox_coefficients.tsv","29_Source_Cox_pairwise_contrasts.tsv","29_Source_specific_Cox_PH_tests.tsv","29_Source_specific_Cox_diagnostics.tsv"],
 "같은 환자의 상관된 분석이며 독립 검증이 아니다. 전체 보정항목과 PH 진단을 함께 확인한다.")
figure(12,"KRAS 두 출처의 일치와 불일치","25_KRAS_group_heatmap.png",
 "동일 환자의 축약 KRAS 범주를 교차표로 비교하고, 별도 판정 가능한 변이 양성 쌍에서 정확한 아형을 비교했다.",
 f"축약 범주 일치 {int(agreement.iloc[0].agreement_n)}/{n}명 ({agreement.iloc[0].agreement_fraction*100:.1f}%), κ={agreement.iloc[0].kappa:.3f}. 정확 비교는 866/870쌍 일치이며 4쌍 불일치이다.",
 ["25_KRAS_group_cross_table.tsv","25_KRAS_detailed_cross_table.tsv","25_KRAS_agreement_summary.tsv","25_KRAS_comparison_status.tsv"],
 "분류 일치도와 예후 연관성은 서로 다른 질문이다. 상세 판정 불가·양쪽 no call은 정확한 변이 일치 분모에서 분리한다.")

figure(13,"Clinical TMB와 Variant count의 개별 분포","27_TMB_source_distributions.png",
 "각 지표를 원래 단위로 그린 히스토그램이다. 임상 TMB는 mutations/Mb, Variant count는 retained MAF 행 수이다.",
 "임상 TMB 836명, Variant count 1,011명을 사용한다. 표본 수가 달라 이후 동일 환자 836쌍의 비교를 추가했다.",
 ["27_TMB_distribution_summary.tsv","27_TMB_analysis_limits.tsv"],
 "가로축 단위와 분모가 달라 분포의 절대 위치를 직접 동등 비교할 수 없다.")
newpage("Table5. 두 수치 지표의 관측률과 분포")
heads=["Cohort / measure","Observed / total","Missing","Median [IQR]","Range"]
rows=[]
for r in burden.itertuples():
    rows.append([r.cohort+"\n"+display(r.metric),f"{int(r.observed_n)}/{int(r.total_n)}",str(int(r.missing_n)),
      f"{r.median:.2f} [{r.q1:.2f}, {r.q3:.2f}]",f"{r.minimum:.2f}–{r.maximum:.2f}"])
table(heads,rows,[2.5,1,.6,1.85,1.05],9)
para(f"동일 환자 836쌍 Spearman ρ={corr.spearman_rho:.4f}; percentile bootstrap 95% CI {corr.lower95:.4f}–{corr.upper95:.4f} (2,000회). 이는 순위 연관성이며 절대 일치도 지표가 아니다.")
para("권장 해석: 임상 TMB는 임상 단위를 가진 주 지표로, Variant count는 별도 탐색 지표로 제시한다. 결측 임상 TMB를 count로 대체하거나 두 출처를 합쳐 하나의 TMB 열로 만들지 않는다. 패널 callable Mb와 필터 정보를 확인하면 정규화 TMB의 절대 차이·일치도 분석을 별도로 설계할 수 있다.")
refs(13,["27_TMB_distribution_summary.tsv","27_TMB_analysis_limits.tsv","27_TMB_MAF_rank_correlation.tsv"])
export_table("Table5_Manuscript_burden_summary.tsv",heads,rows)
figure(14,"동일 환자의 Clinical TMB–Variant count 순위 연관성","27_TMB_MAF_paired_scatter.png",
 "완전쌍 836명, Spearman 상관과 환자 단위 bootstrap 신뢰구간.",
 f"ρ={corr.spearman_rho:.3f}로 강한 순위 연관성이 있다. 값이 함께 증가하는 경향을 보이지만 같은 측정량이라는 증거는 아니다.",
 ["27_TMB_MAF_rank_correlation.tsv","27_TMB_distribution_summary.tsv"],
 "공통 assay/pipeline 유래 여부가 미확인이다. 서로 다른 단위에 identity line이나 Bland–Altman 한계를 적용하지 않았다.")
miss=read("27_TMB_missingness_MAF_count_comparison.tsv").iloc[0]
figure(15,"임상 TMB 관측군과 결측군의 Variant count","27_TMB_missingness_MAF_counts.png",
 "임상 TMB 가용성으로 나누어 MAF Variant count 분포를 비교했다.",
 f"군 간 전체 검정 p={fp(miss.p)}. 이 결과 하나만으로 결측이 무작위(MCAR/MAR)인지 판정할 수 없다.",
 ["27_TMB_missingness_MAF_count_comparison.tsv","27_TMB_missingness_clinical_comparisons.tsv"],
 "관측군과 결측군의 임상 구성도 함께 확인하고, 결측 175명을 자동 보간하지 않는다.")
figure(16,"두 수치 지표의 보정 생존 연관성과 동일 환자 민감도","29_TMB_source_outcome_sensitivity.png",
 "각 지표를 연속형으로 유지한 Cox 모형. 전체 관측군 및 동일 836쌍에서 출처별 paired SD 증가당 HR을 표시한다.",
 "동일 환자 OS 분석에서 Clinical TMB와 Variant count의 HR은 각각 약 1.03/SD이며 p=0.580, 0.585로 유의한 선형 연관성을 확인하지 못했다. 두 측정량의 동등성 또는 예후 효과가 없음을 증명한 결과는 아니다.",
 ["29_Source_specific_Cox_coefficients.tsv","29_Source_specific_Cox_diagnostics.tsv","29_Source_specific_Cox_PH_tests.tsv"],
 "동일 환자 OS 모형의 PH p는 두 지표 모두 약 0.03으로 일정 HR 해석에 주의한다. 재발 완전사례는 동일 TMB 쌍 중 807명이다. 임의 high/low cut-off는 사용하지 않았다.")
figure(17,"임상 변수와 유전체 지표의 탐색적 연관성","28_Clinical_genomic_associations.png",
 "연속형은 Kruskal–Wallis, 범주형은 Fisher 검정. 각 grouping별 BH q를 표시하고 다군 쌍별 검정은 Holm 보정했다.",
 "새 임상 항목을 포함한 연관성의 전체 지도를 제공한다. 개별 유의 조합은 표본 수·결측·효과 방향을 보조표에서 확인한다.",
 ["28_Clinical_genomic_association_tests.tsv","28_Clinical_pairwise_tests.tsv","28_Clinical_source_stratified_tables.tsv","27_Burden_pairwise_tests.tsv"],
 "치료 관련 항목의 연관성은 치료 효과가 아니다. 수많은 탐색 검정과 희소 범주로 인한 우연·불안정성을 고려한다.")
recur=read("22_Observed_recurrence_patterns.tsv")
figure(18,"기록된 재발 환자의 재발 양상","22_Observed_recurrence_patterns.png",
 "재발이 기록된 환자에서 local/distant/both 등 범주의 빈도와 비율을 계산했다.",
 "; ".join(f"{r.Recurrence_pattern}: {int(r.n)}명 ({r.percent:.1f}%)" for r in recur.itertuples())+".",
 ["22_Observed_recurrence_patterns.tsv","22_Recurrence_patterns_by_driver.tsv"],
 "분모는 재발 환자이며 전체 코호트의 누적발생률이나 특정 시점의 재발 위험이 아니다.")
figure(19,"Distant-only 재발에서 관찰된 장기 분포","22_Distant_only_recurrence_sites.png",
 "Distant-only 재발 사례에 정의된 Distant_pattern을 집계했다.",
 "관찰된 원격 재발의 장기 분포를 보여준다. Local+distant 동시 재발은 이 표의 분모에 포함되지 않는다.",
 ["22_Distant_only_recurrence_sites.tsv"],
 "재발 전 사망의 경쟁위험을 반영한 장기별 발생률과 구분한다.")
figure(20,"DNA 복구 관련 유전자의 retained variant 빈도","24_DNA_repair_gene_variants.png",
 "참고논문의 HRD 관련 gene list 및 MMR 유전자에서 retained MAF call이 있는 환자의 비율을 요약했다.",
 "후속 변이 검토 대상을 제시하는 탐색적 기술통계이다. 빈도가 곧 약제 감수성이나 기능적 HRD를 뜻하지 않는다.",
 ["24_DNA_repair_gene_variant_counts.tsv","24_Analysis_definitions_and_limits.tsv"],
 "생식세포/체세포 구분 검증, 병원성, LOH 및 기능적 HRD 자료가 없으므로 HRD-positive로 명명하지 않는다.")


# New requested gene-list oncoplots and clinical profiles.
repair_summary = read("31_Repair_gene_list_summary.tsv")
repair_freq = read("31_Repair_gene_patient_frequencies.tsv")
repair_clin = read("32_Repair_group_clinical_distributions.tsv")
repair_tests = read("32_Repair_group_clinical_tests.tsv")
repair_overlap = read("31_Repair_gene_list_overlap.tsv")
def repair_result(set):
    r = repair_summary[repair_summary.gene_list==set].iloc[0]
    ff = repair_freq[repair_freq.gene_list==set].sort_values("patients",ascending=False).head(3)
    top = "; ".join(f"{q.gene} {int(q.patients)}명 ({q.percent:.1f}%)" for q in ff.itertuples())
    return (f"전체 {int(r.total_patients):,}명 중 목록 내 변이 보유자는 {int(r.variant_positive_n)}명 "
      f"({r.percent:.1f}%)이고, 두 개 이상 유전자에서 변이가 있는 환자는 {int(r.multiple_genes_n)}명이다. "
      f"전체 코호트에서 빈도가 높은 유전자는 {top}이다. "
      "이 그림의 확대 분모는 변이 보유 환자이며, 전체 코호트 빈도와 구분한다.")
for num,set,title in [(21,"HRD","HRD 관련 유전자 변이의 환자별 분포"),
                      (22,"MMR","MMR 유전자 변이의 환자별 분포")]:
    figure(num,title,f"31_{set}_oncoplot_variant_positive_patients.png",
      f"McIntyre 2020 목록의 {'18' if set=='HRD' else '4'}개 유전자를 사용했다. "
      "열은 한 환자, 행은 유전자, 색은 MAF 변이 종류이며 임상 특성을 상단에 정렬했다. 한 셀의 복수 변이 종류는 나눠 표시했다.",
      repair_result(set),
      ["31_Repair_gene_list_definition.tsv","31_Repair_gene_patient_frequencies.tsv",
       "31_Repair_variant_class_counts.tsv","31_Repair_gene_list_summary.tsv"],
      "변이 미보유 환자까지 포함한 전체 분모 oncoplot은 "
      f"31_{set}_oncoplot_all_patients.png이다. 변이 종류·병원성·biallelic 소실·기능 장애는 서로 다른 개념이며, 이 그림으로 HRD-positive 또는 dMMR을 판정하지 않는다.")
newpage("Table6. HRD 관련 및 MMR 유전자별 변이 빈도")
heads=["Gene list","Gene","Patients, n (%)","Variant rows"]
rows=[[r.gene_list,r.gene,f"{int(r.patients)} ({r.percent:.1f})",str(int(r.variant_rows))]
      for r in repair_freq.itertuples()]
table(heads,rows,[1.2,1.6,2.5,1.7],9)
para("환자 비율의 분모는 1,011명이다. 같은 환자가 여러 유전자나 여러 변이 행에 기여할 수 있어 합산 환자 수와 변이 행 수는 다르다. 0은 이 MAF에 유지된 변이가 없다는 뜻이며 패널의 완전한 음성 판정이 아니다. FAM175A는 논문의 표기이며 ABRAXAS1 동의어를 함께 처리한다.",9.5)
refs(21,["31_Repair_gene_list_definition.tsv","31_Repair_gene_patient_frequencies.tsv","31_Repair_gene_list_summary.tsv"])
export_table("Table6_Manuscript_repair_gene_frequencies.tsv",heads,rows)
def clinical_result(set):
    grouping=set+"_variant_status"
    tests=repair_tests[repair_tests.grouping==grouping]
    sig=tests[tests.q_BH<.05]
    result="BH 보정 후 유의한 임상/수치 항목은 "+(
      "; ".join(f"{display(r.variable)} (q={fp(r.q_BH)})" for r in sig.itertuples()) or "없었다")+". "
    for var in ["TMB","MAF_variant_count"]:
        z=repair_clin[(repair_clin.grouping==grouping)&(repair_clin.variable==var)]
        result+=display(var)+" 중앙값은 "+" vs ".join(
            f"{r['median']:.2f}" for _,r in z.set_index("group").loc[["No retained variant","Variant detected"]].iterrows())+" (미검출군 vs 변이 보유군)였다. "
    return result
for num,set in [(23,"HRD"),(24,"MMR")]:
    figure(num,f"{'HRD 관련' if set=='HRD' else 'MMR'} 유전자군 변이 보유 여부에 따른 임상 특성",
      f"32_{set}_clinical_distributions.png",
      "각 유전자군에서 하나 이상의 retained variant가 있는 환자와 없는 환자를 비교했다. "
      "막대는 각 군 내 임상 범주의 비율이며 결측도 분모에 포함했다. 연속형은 Wilcoxon, 범주형은 Fisher, 항목군 내 BH 보정을 적용했다.",
      clinical_result(set),
      ["32_Repair_group_clinical_distributions.tsv","32_Repair_group_clinical_tests.tsv",
       "32_Repair_group_clinical_pairwise_tests.tsv","32_Per_gene_clinical_distributions.tsv"],
      "Variant count 자체가 많으면 목록 내 변이가 하나 이상 발견될 기회도 늘므로 count 차이는 부분적으로 정의에 내재된 연관성이다. "
      "TMB와의 연관성 역시 HRD 또는 dMMR 기능 장애의 증거가 아니다. 임상 분포의 비유의 차이는 두 군의 동등성을 증명하지 않는다.")
    newpage(f"Table7. 유전자군별 임상 특성 — {'A. HRD 관련 유전자군' if set=='HRD' else 'B. MMR 유전자군'}")
    variables=["Age","TMB","MAF_variant_count","Neoadjuvant","T_stage","N_stage","M_stage","Differentiation","R_status","MSI"]
    heads=["Characteristic","No retained variant","Variant detected","BH q"]
    rows=[]
    for var in variables:
        z=repair_clin[(repair_clin.grouping==set+"_variant_status")&(repair_clin.variable==var)]
        q=repair_tests[(repair_tests.grouping==set+"_variant_status")&(repair_tests.variable==var)].iloc[0].q_BH
        categories=sorted(z.category.unique())
        for j,cat in enumerate(categories):
            vals=[]
            for group in ["No retained variant","Variant detected"]:
                a=z[(z.group==group)&(z.category==cat)]
                if a.empty: vals.append("0 (0.0%)"); continue
                rr=a.iloc[0]
                vals.append(f"{rr['median']:.1f} [{rr.q1:.1f}, {rr.q3:.1f}]; m={int(rr.missing_n)}"
                  if cat=="Continuous" else f"{int(rr.n)} ({rr.percent:.1f}%)")
            rows.append([display(var)+(f": {cat}" if cat!="Continuous" else ""),*vals,fp(q) if j==0 else "—"])
    table(heads,rows,[2.35,1.95,1.95,.75],8)
    ss=repair_summary[repair_summary.gene_list==set].iloc[0]
    para(f"군별 전체 N: No retained variant {int(ss.no_retained_variant_n)}명, Variant detected {int(ss.variant_positive_n)}명. "
      "연속형은 median [IQR], m=결측. 비율은 결측 포함 각 군의 전체 N 기준. q는 각 행의 범주가 아니라 해당 변수 전체 검정의 BH 보정값이다. "
      "나머지 임상 항목과 각 유전자별 분포는 동일 Supplementary table에 있다.",9)
    refs(num,["32_Repair_group_clinical_distributions.tsv","32_Repair_group_clinical_tests.tsv","32_Per_gene_clinical_distributions.tsv"])
    export_table(f"Table7{set}_Manuscript_clinical_profiles.tsv",heads,rows)
figure(25,"HRD 관련 및 MMR 유전자군 변이의 중복 분포","31_HRD_MMR_variant_overlap.png",
  "두 유전자 목록에서 하나 이상 변이가 있는지를 조합해 서로 배타적인 네 군을 만들었다.",
  "; ".join(f"{r.Repair_pattern} {int(r.n)}명 ({r.percent:.1f}%)" for r in repair_overlap.itertuples())+
  ". 두 목록 동시 보유군은 단일 유전자군 보유군과 분리해 임상 분포 및 모든 쌍별 비교를 제공했다. "
  "그러나 이 중복은 같은 환자에서 목록 내 변이가 함께 발견되었다는 의미이며 두 DNA 복구 경로의 기능적 동시 결핍을 증명하지 않는다.",
  ["31_Repair_gene_list_overlap.tsv","32_Repair_group_clinical_distributions.tsv","32_Repair_group_clinical_pairwise_tests.tsv"],
  "임상 항목의 네 군 비교는 전체 검정과 Holm 보정 쌍별 검정을 구분한다. 환자가 두 막대에 중복 집계되지 않는다.")

newpage("원고 작성 시 해석과 공개 범위")
para("전체 log-rank p가 작다는 것은 모든 그룹이 서로 다르다는 뜻이 아니다. 특정 쌍의 차이는 해당 pairwise Holm p에서 판단한다. 유의하지 않은 결과는 동등성을 입증하지 않으며, 희소 그룹의 넓은 CI와 검정력을 고려해야 한다.")
para("다변량 Cox 결과는 조정된 관찰 연관성이다. 비례위험 검정이 유의한 경우 시간에 따라 효과가 달라질 수 있어 일정 HR만으로 요약하는 데 한계가 있다. 36개월 RMST는 미보정 보조 분석으로, 조정 Cox의 대체가 아니다.")
para("재발 관찰기간이 정의되지 않은 M1과 일부 결측 사례는 재발 분석에서만 제외했다. 검체·시퀀싱 시점이 없어 생존자 선택 편향을 배제할 수 없다. 새 임상 사전과 원자료 확인 후 논문 Methods에 표본 추출 경로를 명시해야 한다.")
para("첨부 가능한 Supplementary table은 환자 식별자를 포함하지 않는 집계표만 선정했다. 파일명에 _private가 있거나 patient-level인 자료는 공개용 보조표에서 제외했다. 집계표도 희소 셀 공개 정책과 연구 승인 범위를 연구자가 검토해야 한다.")
para("Fig/Table 번호는 이 보고서의 원고 구성안이다. Supplementary table의 앞 번호는 최초 연결된 Figure 묶음, 뒤 번호는 그 묶음 안의 순서이다. 동일 파일을 여러 Figure/Table에서 참조하면 최초 부여한 번호를 유지한다.")
doc.add_heading("재현 방법",1)
para("저장소 루트에서 Rscript --vanilla run_all.R --config=config/local.R 실행 후, python3 scripts/build_results_report.py --results outputs/updated_v3_20260928 을 실행한다. Python 의존성: python-docx, pandas, Pillow. 원본 그림은 결과 폴더 PNG를 사용하며 Word에는 읽기용 축소본을 삽입했다.",9.5)
doc.add_heading("참고문헌 및 방법 문서",1)
para("1. McIntyre CA, et al. Alterations in Driver Genes Are Predictive of Survival in Patients With Resected Pancreatic Ductal Adenocarcinoma. Cancer. 2020;126:3939–3949. doi:10.1002/cncr.33038.",9)
para("2. Campbell BA, et al. The Mutational Status of Driver Genes in Patients With Resected Pancreatic Ductal Adenocarcinoma Is Associated With Pathological Characteristics and Overall Survival. Ann Surg. 2025;282. doi:10.1097/SLA.0000000000006794.",9)
para("3. R survival::survdiff; stats::p.adjust (Holm). https://stat.ethz.ch/R-manual/R-devel/library/survival/html/survdiff.html ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html",8.5)
para("4. Dunn 양측 검정 및 ties 처리 검증 기준: rstatix::dunn_test. https://rpkgs.datanovia.com/rstatix/reference/dunn_test.html",8.5)
# Exhaustive supporting-table registry and all-figure inventory are machine-readable.
pd.DataFrame(supp).to_csv(out/"30_Supplementary_table_index.tsv",sep="\t",index=False)
figs=sorted(out.glob("*.png"))
pd.DataFrame({"figure_file":[p.name for p in figs],
              "role":["본문 또는 추가 탐색 그림; 기존 결과 파일명 유지"]*len(figs)}).to_csv(out/"30_Figure_inventory.tsv",sep="\t",index=False)
manifest=[]
for name in sorted({r["file"] for r in supp}|{p.name for p in figs}):
    manifest.append({"file":name,"sha256":hashlib.sha256((out/name).read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(out/"30_Report_source_manifest.tsv",sep="\t",index=False)
para("전체 Supplementary table 번호 목록: 30_Supplementary_table_index.tsv",9)
para("전체 그림 파일 목록(본문 외 추가 분석 포함): 30_Figure_inventory.tsv",9)
para("보고서 생성 시점 근거 파일 체크섬: 30_Report_source_manifest.tsv",9)
in_last_section = False
for p in doc.paragraphs:
    if p.text == "원고 작성 시 해석과 공개 범위": in_last_section = True
    if in_last_section and not p.style.name.startswith("Heading"):
        p.paragraph_format.line_spacing = 1.04
        p.paragraph_format.space_after = Pt(4)
        for r in p.runs:
            if r.font.size is None or r.font.size.pt > 9.5: r.font.size = Pt(9.5)
target.parent.mkdir(parents=True,exist_ok=True)
doc.save(target)
print(f"Saved {target}; embedded figures=25; editable tables={len(doc.tables)}; supporting tables={len(supp)}")
