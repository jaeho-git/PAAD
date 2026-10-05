# 논문용 분석과 상세 보고서

## 연구 구성

| 구분 | 내용 |
| --- | --- |
| Figure 1 | 코호트, oncoplot, 선택한 KRAS 아형의 분포 |
| Figure 2 | 주요 유전자와 병리 소견의 연관성 |
| Figure 3 | 선행 platinum 노출의 임상·유전체 맥락, OS 및 기록된 재발 |
| Figure 4 | 재발양상과 재발 후 생존 |
| Figure 5 | 주요 driver 및 KRAS 아형의 예후 연관성 |
| Table 1 | 전체 코호트의 임상·유전체 특성 |
| Table 2 | 메인 생존모형의 HR, 95% CI, 원 p와 보정 p |
| Supplementary Figure 1–5 | HRD 관련 유전자, MMR 유전자, Reported TMB, SBS 가능성 및 집단 합산 민감도 분석 |
| Input 파일 | 각 그림·표의 계산 근거인 집계 TSV. Supplementary Table로 호칭하지 않음 |

임상/MAF 자료의 일치도·출처 비교는 내부 QC이며 논문용 결과 및 보충자료에 포함하지 않습니다. 원자료 파일명은 `Input 1. ...`로 표시하고, 실제 결과표와 구분합니다.

## 분석 원칙

- M0와 M1을 모두 사용하며 중복 patient_id의 모든 행을 제외합니다. 현재 자료에서는 1,015행 중 4행을 제외한 1,011명입니다.
- 현재 원고는 업데이트된 엑셀의 `KRAS_subtype_MAF`을 직접 선택하여 사용합니다. 아형, 이분형 KRAS, 주요 변이 유전자 수에도 같은 정의를 적용합니다. 출처 간 일치도는 원고에 넣지 않습니다.
- Oncoplot의 유전자 행은 **KRAS를 포함하여 모두 MAF의 실제 변이 유형**을 표시합니다. KRAS는 메인 그림에서 항상 포함합니다. 상단 KRAS 아형 annotation은 선택한 단일 출처를 사용하며, 유전자 행과 annotation의 출처를 그림 각주에 표시합니다. 이 복원은 병리·생존분석의 KRAS 정의를 바꾸지 않으며 별도의 출처 일치도·비교 분석을 논문에 추가하지 않습니다.
- 상단 **Variant count**: 해당 환자의 모든 유지 변이 수. HRD/MMR 그림에서도 해당 목록에 한정하지 않습니다.
- 오른쪽 **Variants per gene**: 표시된 환자들에게서 해당 유전자의 변이 수. 왼쪽 **Patients n (%)**: 해당 변이 보유 환자 수·비율입니다.
- 임상정보 track은 `R/oncoplot_annotations.R`에서 명시적으로 관리하며 필수 컬럼이 빠지면 중단합니다. 메인 oncoplot에는 기존 14개 임상항목과 선택된 KRAS 아형, Reported TMB를 표시합니다. HRD/MMR에는 기존 11개 항목과 Reported TMB를 표시합니다. 임상정보의 결측은 회색으로 표시하며 0으로 대체하지 않습니다.
- Annotation 색상은 `R/annotation_palettes.R`에서 항목·범주별로 고정합니다. 서로 다른 범주에 같은 색상 코드를 재사용하지 않으며, 연속형 Age·BMI·CA19-9·CEA·종양 크기·TMB는 각각 별도 색상 척도를 사용합니다. 같은 항목은 메인·HRD·MMR 및 개별 oncoplot에서 같은 색을 사용합니다. 실제 색상 코드는 `Input_*_clinical_tracks.tsv`의 `color_mapping`에 기록합니다. 결측 회색은 음성·No와 구분합니다.
- Forest plot은 `R/forest_display.R`에서 비교 대상, HR·95% CI 그림, 통계값 열을 행별로 정렬합니다. p·q를 축 라벨에 붙이지 않습니다. 원 p값은 `p`, BH 보정값은 `q (BH)`, KRAS 쌍별 Holm 보정값은 `Holm p`로 구분합니다. 추정치와 검정 결과 자체는 바꾸지 않습니다.
- 메인 그림의 CA19-9·CEA 색상은 기존 그림과 같이 80백분위수에서 포화됩니다. 범례의 `>=`는 그 이상의 값을 같은 색으로 표시한다는 뜻이며 원래 수치는 변경하지 않습니다. 집계 범위와 결측 수는 `Input_*_clinical_tracks.tsv`에서 확인할 수 있습니다.
- **Figure 1a**는 상자와 화살표로 구성한 환자 선정 흐름도, **Figure 1b**는 변이·임상정보 oncoplot, **Figure 1c**는 KRAS 아형 분포입니다. 흐름도는 전체 기록에서 중복 ID의 모든 기록을 제외한 뒤, 최종 코호트에 M0와 M1이 모두 포함됨을 보여줍니다. `R/cohort_flow.R`의 동일 집계·배치 정보를 사용하여 R 그림과 PPTX의 편집 가능한 흐름도를 생성합니다. PNG·PDF 그림에는 패널 문자 없이 분석 제목만 표시합니다. 합본에도 a/b/c 문자를 덧붙이지 않습니다. PPTX는 메인 그림의 모든 패널을 개별 슬라이드로 배치하고, 편집 가능한 슬라이드 제목에 `Figure 1a`, `Figure 1b`, `Figure 2a`처럼 번호를 표시합니다. DOCX의 그림 번호·캡션과 그림 파일명은 식별을 위해 유지합니다.
- 같은 환자·유전체 사건의 중복 annotation은 제거하지만, 서로 다른 환자의 같은 사건은 각각 집계합니다. 검출되지 않음은 확인된 wild type과 다릅니다. 유전자별 검사 가능 범위는 별도로 확인해야 합니다.
- TMB는 임상자료의 **Reported TMB**를 사용합니다. 확인되지 않은 단위, TMB-high 임계값, Variant count에 의한 대체를 적용하지 않습니다.
- 전체 log-rank와 쌍별 log-rank를 구분하고 쌍별 p에 Holm 보정을 적용합니다. Cox 대비는 별도로 보고합니다. 전체 쌍별 결과는 TSV로 보존합니다.
- 사망 포함 DFS/RFS가 아닌 **기록된 재발까지의 시간**을 사용합니다. 비례위험 위반을 함께 보고하고 일정한 HR 해석을 제한합니다.
- HRD 관련/MMR 유전자 변이는 기능적 HRD/dMMR 또는 치료 반응을 의미하지 않습니다. LOH, CNV, purity, cellularity 분석은 포함하지 않습니다.
- `Recurrence_pattern=Distant only`는 `Distant_pattern`으로 Liver only, Lung only, Other distant의 3군으로 분류합니다. Peritoneum only, Multiple distant 및 나머지 원격부위는 Other distant입니다. 기존 3군 분석과 Local only·Local+distant를 포함한 확장 5군 분석을 모두 보존합니다.
- SBS 분석은 targeted-panel MAF의 SNV만 사용한 탐색적 가능성 평가입니다. 환자별 사전기준은 context-counted SNV≥10, reconstruction cosine≥0.90, target count≥5, proportion≥0.20, bootstrap stability≥0.80입니다. 현재 자료에서는 1,011명 중 193명만 SNV≥10이었고 최고 cosine이 0.899여서 통과자가 0명입니다. 따라서 platinum/HRD/MMR process별 positive rate·burden·proportion 9개 지표는 0이 아니라 추정 불가(NE)입니다. 집단 합산 constrained refit은 통합 process 6개 검정과 SBS31·35·3·6·14·15·20·21·26·44의 개별 20개 검정을 분리하고 각각 Holm 보정합니다. Standard pooled assignment와 forced constrained refit도 별도 감사표·그림으로 비교합니다. 이 결과들은 방법론적 민감도 분석일 뿐 환자별 signature 또는 치료 유발 효과의 근거가 아닙니다.

## R 분석 실행

저장소 최상위 폴더에서 실행합니다. 기존 `config/local.R`을 덮어쓰지 말고 아래 항목을 확인합니다.

```r
clinical_file = "data/raw/260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx"
analysis_profile = "manuscript"
kras_source = "workbook_maf" # Main!KRAS_subtype_MAF를 논문 분석에 사용
updated_cohort = "all_unique"
manuscript_subdir = "manuscript_extension_20261005"
signature_analysis = TRUE # 별도 SigProfiler 환경과 GRCh37 reference가 준비된 경우
clinical_extensions = FALSE
pretreatment_measurements = TRUE
```

```sh
Rscript --vanilla run_all.R --config=config/local.R
# 또는 동일 분석만 직접 실행
Rscript --vanilla scripts/manuscript_analysis.R --config=config/local.R
```

출력은 `output_dir/manuscript_subdir/` 안에 저장됩니다.

```text
figures/       PNG 및 PDF 그림
tables/        집계 결과 TSV
deliverables/  생성된 DOCX 및 PPTX
private/       환자별 중간자료, provenance, 렌더링·검증 자료
```

현재 분석은 `260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx`의 `Main` 시트를 사용합니다. 기존 결과 폴더를 자동 삭제하지 않습니다. 다른 출처나 조건을 실행할 때는 새 `manuscript_subdir`를 사용해 이력을 분리합니다. 모든 실제 자료와 출력은 Git에서 제외됩니다.

## Word와 PowerPoint 생성

Word 생성에는 Python의 `python-docx`와 `Pillow`가 필요합니다. 먼저 R 분석이 성공해야 합니다. 아래 경로는 실제 결과 폴더로 바꿉니다.

```sh
python3 scripts/build_updated_manuscript_report.py \
  --results outputs/updated_v3_20260928/manuscript_extension_20261005_kras_maf \
  --output outputs/updated_v3_20260928/manuscript_extension_20261005_kras_maf/deliverables/PAAD_Updated_Manuscript_Report_20261005_KR_Explained_v4.docx
```

그림별 분석 방법, 실제 수치, 임상·생물학적 해석, 한계를 구분하고 공통 원고 `private/document_build/content.json`도 생성합니다. 이번 원고에는 출처 확정 전 검토 문구를 넣지 않습니다.

PowerPoint는 Codex workspace의 `@oai/artifact-tool`과 presentation 검증 도구를 사용합니다. **일반 Node만 설치한 환경에서는 그대로 실행되지 않습니다.** 실행환경 경로를 확인한 후 다음 환경변수를 지정합니다.

```sh
export PAAD_RUNTIME_NODE_MODULES="/실제/runtime/node_modules"
export PAAD_RUNTIME_PYTHON="/실제/runtime/python3"
export PAAD_PRESENTATIONS_SKILL="/실제/presentations/skill/폴더"
node scripts/build_updated_manuscript_slides.mjs \
  outputs/updated_v3_20260928/manuscript_extension_20261005_kras_maf \
  outputs/updated_v3_20260928/manuscript_extension_20261005_kras_maf/deliverables/PAAD_Updated_Manuscript_Presentation_20261005_KR_Explained_v5.pptx
```

PPTX는 그림 슬라이드와 상세 해석 슬라이드를 구분합니다. 표·본문은 편집 가능한 개체이며, Pretendard 폰트가 필요합니다. 같은 최종 PPTX가 이미 있으면 실수로 덮어쓰지 않도록 중단하므로 마지막 인자로 새 파일명을 지정하십시오. 생성 뒤 DOCX 모든 페이지와 PPTX 모든 슬라이드의 렌더링 검토가 필요합니다. macOS headless Word 렌더링은 사용자 폰트 디렉터리를 포함한 fontconfig 설정이 필요할 수 있습니다.

## 검증

```sh
Rscript --vanilla tests/test_oncoplot_counts.R
Rscript --vanilla tests/test_oncoplot_annotations.R
Rscript --vanilla tests/test_cohort_flow.R
Rscript --vanilla tests/test_figure_titles.R
Rscript --vanilla tests/test_forest_display.R
Rscript --vanilla tests/test_pairwise.R
Rscript --vanilla tests/test_clinical_validation.R
Rscript --vanilla tests/test_source_genomics.R
Rscript --vanilla tests/run_tests.R
python3 tests/verify_manuscript_outputs.py --results outputs/updated_v3_20260928/manuscript_extension_20261005
```

마지막 검사는 단일 출처, 집계 일관성, 코호트 구성, 그림/표 구성, Input 파일 존재, 보고서에 내부 비교가 없는지를 점검합니다. 소프트웨어 검증은 임상적 타당성·외부 검증을 대체하지 않습니다. 투고 전에는 검사 범위·버전, 변이 병원성 및 연구윤리·연구설계 정보를 별도로 확정해야 합니다.

## 임상적 가치 확장 (2026-09-29)

상세 계획은 `docs/clinical_value_analysis_plan.md`에 기록했습니다. 기존 결과를 본 뒤 진행한 탐색적 확장으로, 사전등록 확증 분석으로 표현하지 않습니다.

- 측정 시점: 연구자 확인에 따라 임상 검사·분자정보는 진단 후 선행치료 전입니다. 수술 병리와 절제연은 별도 수술 정보입니다. OS 시작점은 수술로 유지합니다.
- 주 분석: 연령 자연 스플라인, 성별, T/N/M, 선행치료, LVI/PNI, 절제연 보정. 분화도의 비례위험 위반을 고려해 분화도별 기저위험을 허용합니다.
- KRAS HR 표에는 비교군과 기준군, 사건을 명시합니다. HR > 1은 비교군의 사망 순간위험률이 높다는 뜻이며 누적 사망확률비가 아닙니다.
- 같은 complete-case 환자에서 임상병리 단독 및 유전체 추가 모형을 비교합니다. 500회 환자 bootstrap으로 모형 전체를 재적합하고 낙관성 보정 C-index·Brier score와 OOB 보정도를 산출합니다. OOB 성능 차이의 중앙 95% 범위는 보정 효과의 신뢰구간이 아닙니다.
- 절대 효과는 동일 공변량 분포로 표준화한 12/24/36개월 생존확률과 차이로 제시합니다. 선행치료 상호작용, M0 및 CA19-9/ASA 추가 보정은 민감도 분석입니다.
- 추가 예측 개선이 작거나 없는 결과도 그대로 보고하며, 관찰적 연관성만으로 치료 변경을 권고하지 않습니다.
