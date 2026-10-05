# 확정 매핑을 적용한 새 임상자료 분석

> 현재 논문용 최종 실행은 `data/raw/260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx`의 `Main` 시트를 사용하며, 결과는 `outputs/updated_v3_20260928/manuscript_extension_20261005/`에 저장합니다. 이 문서의 18–32 번호 결과는 이전 확장 분석 이력이고, 현재 메인 Figure 1–5 구성은 [manuscript_analysis.md](manuscript_analysis.md)를 따릅니다.

## 실행

저장소 루트에서 실행합니다.

```sh
Rscript --vanilla run_all.R --config=config/local.R
# 추가 분석만 독립 실행
Rscript --vanilla scripts/source_comparisons.R --config=config/local.R
# 검증
Rscript --vanilla tests/test_curated_input.R
Rscript --vanilla tests/test_source_genomics.R
Rscript --vanilla tests/verify_curated_outputs.R --config=config/local.R
```

설정은 `clinical_schema="curated_v3"`, `clinical_sheet="Main"`, `updated_cohort="all_unique"`입니다.
Main과 Data_dictionary를 읽고, 사용자가 완성한 매핑의 의미를 코드에 명시적으로 반영했습니다.
매핑 파일의 빈 지시 칸은 해당 새 컬럼명을 채택한다는 뜻입니다. 임상 결측값은 별개이며 NA로 유지합니다.
매핑 XLSX와 원본 XLSX/MAF는 읽기만 합니다. 사전/매핑의 자유서술문을 실행 코드로 평가하지 않습니다.

## 대상 및 입력 우선순위

- 새 workbook 전체에서 **중복 patient_id의 모든 행만 제외**합니다. 임의로 첫 검체를 고르지 않습니다.
- **M0와 M1 모두 주 분석에 포함**합니다. OS Cox에는 M_stage를 보정하고 M0-only 결과를 보조 민감도로 제공합니다.
- 기록된 재발 endpoint가 정의되지 않은 M1은 재발 분석에서만 빠집니다. 이를 전 코호트 M0 제한으로 해석하면 안 됩니다.
- 각 분석은 필요한 항목의 complete cases를 사용하고 N·사건 수·결측 수를 기록합니다.
- 중복 포함 전체 검체 행은 입력 감사와 별도 기술통계에만 남습니다.
- 이전 임상자료와 TXT는 감사 비교용입니다. 새 임상값이나 코호트를 덮어쓰지 않습니다.

## 컬럼과 자료형

| 항목 | 적용 |
| --- | --- |
| Age, Tumor_size | 연속형. Age_group은 ≤55, >55 및 <70, ≥70 |
| ASA, T_stage, N_stage, LN_positive, M_stage | 숫자로 저장되어 있어도 범주형 |
| AJCC_stage | 새 값. Stage_Group은 I/II/III/IV |
| Preop_platinum_exposure | 표기 약어 Preop_platinum |
| Operation | 표기 약어 OP |
| Adjuvant_treatment / chemotherapy / radiotherapy | Adjuvant_Tx / Adjuvant_CT / Adjuvant_RT |
| KRAS_subtype | 임상 보고서 원문 유지 |
| KRAS_subtype_raw | MAF 유래 세부 변이 문자열 |
| TMB | 임상 mutations/Mb. 결측 자동 대체 없음 |
| MAF_variant_count | 유지된 MAF 행의 검체별 개수. mutations/Mb 아님 |

신규 스크립트는 새 컬럼명을 사용합니다. 기존 1–17 분석의 연결을 위해 NAC, Size, T/N/N_status,
Stage, RM, OS_m/survive, RFS_m/Recur 등의 호환 별칭을 R/data.R에 명시적으로 남겼습니다.
KRAS_report와 TMB_report도 같은 임상값의 호환 별칭일 뿐, 새 출처가 아닙니다.
LVI/PNI y/n은 Positive/Negative, 분화도는 대문자로 표시되며 의미와 결측은 보존됩니다.
Classification은 legacy 호환용 PDAC 상수이며 새 코호트 필터가 아닙니다. OS_d/RFS_d는 만들지 않습니다.
Adjuvant_RT 빈칸은 No가 아닙니다. MSI의 원문 범주도 임의로 합치지 않습니다.

## KRAS 두 출처 분석

두 출처 각각 빈도, 임상·병리 비교, OS와 기록된 재발 KM, 보정 Cox를 만듭니다.
기존과 같은 G12D/G12V/G12R/Other/no call 묶음 및 G12D 대 기타 변이 비교를 별도로 적용합니다.
생존분석의 희귀 변이는 Other로 묶되, 세부 원문은 별도 빈도표와 환자별 표에 모두 남습니다.

R/genomics.R에 공통 MAF 판정 규칙이 있습니다. GRCh37 chr12의 기존 G12/G13 좌표를 적용하고,
단백질 주석이 있으면 명시적인 표기를 우선합니다. 현 MAF의 Q61 등 애매한 좌표는
`Q61 (unresolved)`로 기록하며 특정 아미노산 변이라고 단정하지 않습니다.
다중 변이의 쉼표/세미콜론과 순서를 정규화한 비교용 열을 따로 만듭니다.

- 25: 축약 범주 교차표·heatmap, 상세 문자열 교차표·heatmap, 불일치 유형, 개별 검토 목록.
- 정확한 subtype 일치율은 두 출처가 변이를 보고하고 MAF 상세형이 판정 가능한 쌍만 분모로 합니다.
- 양쪽 no call은 별도 유형입니다. 보고서 WT/Not detected나 MAF no call은 생물학적 WT를 증명하지 않습니다.
- 26: 출처별 빈도, KM, log-rank/BH q, median/95% CI, risk table.
- 28/29: 출처별 임상 연관성 및 보정 생존모형.

## TMB 비교

임상 TMB와 MAF 변이 개수는 단위가 다른 지표입니다. 분석영역 Mb와 필터·검증된 assay 규칙이
없으므로 `TMB_MAF` mutations/Mb를 만들어내지 않습니다. 두 값의 비례관계로 분모를 역산하지 않습니다.

- 전체 관측군 및 **동일 paired 환자군**의 분포·결측·요약표.
- paired scatter와 Spearman rho, 환자 단위 2,000회 bootstrap 95% CI.
- paired ECDF, 출처별 표준화 heatmap(단위 변환 아님).
- 임상 TMB 결측군과 관측군의 임상 특성 및 MAF count 비교.
- 각 지표의 임상군별 boxplot, 연속형 변수와의 Spearman, BH q.
- OS 및 기록된 재발 Cox를 원 단위와 paired cohort에서 계산한 출처별 SD 단위로 제공.
- 전체 관측값 모형과 동일 paired cohort 모형을 구분해 표본 구성의 영향을 점검.

정규화된 같은 단위 측정치가 없어 Bland–Altman/절대 일치도는 산출하지 않습니다.
자동 보간·대체, 임의 TMB high/low cut-off도 적용하지 않습니다.
높은 상관관계는 교환 가능성 또는 독립 검증의 증거가 아닙니다.
두 지표의 assay/pipeline이 얼마나 공유되는지 확인해야 최종 지표 선택이 가능합니다.

## 참고논문과 추가 결과

| 근거/목적 | 결과 | 접두사 |
| --- | --- | --- |
| McIntyre / Campbell Table 1 | 전체 검체·주 코호트·Neoadjuvant별 임상표 | 18 |
| Campbell Table 2 | 4개 유전자와 병리의 Neoadjuvant 층화 generalized CMH | 19 |
| 두 논문 Cox 분석 | 단변량·보정 OS, 확장 모형, PH 진단, forest | 20 |
| 두 논문 KM / McIntyre Table 2 | 유전자, driver 0/1/2/3+, TP53 truncation, KRAS/TP53, G12D, median, 1/3/5년 OS | 21 |
| PH 가정 보조 검토 | 비보정 36개월 RMST 및 차이 | 21 |
| 재발 양상 | 관측 Local/Distant/Both, distant-only 세부 부위 | 22 |
| 출처 교차 점검 | KRAS 축약 범주와 TMB/MAF count 초기 요약 | 23 |
| McIntyre DNA repair | HRD 관련/MMR retained MAF variant 빈도 | 24 |
| 이번 확정 사항 | 상세 출처 비교·별도 분석·신규 변수 분포와 연관성 | 25–29 |

출처: McIntyre et al. Cancer 2020, DOI 10.1002/cncr.33038;
Campbell et al. Ann Surg 2025, DOI 10.1097/SLA.0000000000006794.
동일 논문의 완전한 재현이라고 주장하지 않습니다.

## 통계 및 제한

- 주 Cox 보정: Age/10, Sex, T_stage, Neoadjuvant, M_stage. M_stage가 상수인 재발 분석은 해당 항목 제외.
- 확장 문헌 민감도: N, LVI, PNI, R_status 추가. source-specific 결과에는 PH 및 연속형 비선형성(3-df natural spline 대 선형) 진단.
- BH는 명시된 각 검정 family 안에서 적용. 기존 그림의 nominal p는 탐색적으로 해석.
- 새 범주형 변수는 분포·유전자·KRAS 및 TMB 비교에 포함합니다. 매우 희소하거나 한 범주뿐인 검정은 NA와 사유를 표시.
- Adjuvant 치료 시점이 없으므로 치료효과를 주장하는 baseline Cox를 만들지 않습니다.
- 사망을 사건에 포함한 DFS/RFS, CIF/Gray/Fine–Gray는 미구현: 재발 전 사망의 관측 정의가 미확정.
- LOH/CNA, germline, pathogenicity/VUS, 기능적 HRD, Clavien-Dindo, 상세 치료·NGS 시점이 없어 해당 분석은 보류.
- TMB와 달리 Tumor_size·CA19-9/CEA 단위, AJCC 판수 등은 원자료 정의 확인 전 임의 지정하지 않음.
- NGS 시점 부재에 따른 선택/생존자 편향을 교정한 연구가 아님.

## 저장·재실행

결과는 설정의 `outputs/updated_v3_20260928/`에 저장합니다. 파일 이름이 같은 결과는 재실행 시 갱신하지만,
일반 run_all.R은 폴더를 자동 삭제하지 않습니다. 이번 요청에 한해 이전 폴더를 복구 가능하게
이동한 뒤 비운 폴더에 새 결과를 생성했습니다. 매핑 검토 폴더는 삭제하지 않았습니다.

실제 자료, 결과 그림/표, 환자별 private TSV/XLSX와 config/local.R은 Git 제외 대상입니다.
공유되는 것은 코드·문서·합성 예제입니다. 환자별 목록을 공개 저장소나 논문 부록에 올리면 안 됩니다.

## 추가: Variant count·모든 그룹 쌍·Word 보고서

표시명은 Variant count로 통일했습니다. 내부 열 MAF_variant_count와 기존 파일명은 추적·호환성을 위해 유지합니다.

전체 p가 유의한 그림만 골라 검정하지 않고 모든 그룹 쌍을 계산합니다. 유의하지 않은 쌍과 추정 불가 사유도 남깁니다.

| 분석 | 쌍별 방법 | 저장 |
| --- | --- | --- |
| KM | log-rank chi-square/df, n/사건수 | 14–17 개별 pairwise_logrank, 21_OS_pairwise_logrank, 26_KRAS_pairwise_logrank TSV |
| 연속형 다군 | pooled ranks/tie-corrected 양측 Dunn Z; 2군 Wilcoxon W | 4_Variant_count_pairwise, 27_Burden_pairwise_tests, 28_Clinical_pairwise_tests TSV |
| 범주형 | Fisher 2×2 exact; 더 많은 결과 범주에서는 10,000회 Monte Carlo | 4_Frequency_pairwise, 13_Driver_pathology_pairwise_Fisher, 28_Clinical_pairwise_tests TSV |
| 층화 병리 | 두 범주씩 neoadjuvant 층화 CMH | 19_Pathology_pairwise_CMH.tsv |
| 다범주 Cox 노출 | 동일 모형 계수·공분산의 Wald Z contrast | 20_Cox_pairwise_contrasts.tsv, 29_Source_Cox_pairwise_contrasts.tsv |
| 36개월 RMST | 군 간 평균 생존시간 차이와 SE, 양측 Wald | 21_RMST_pairwise_contrasts.tsv |

Holm family는 한 변수×grouping / 한 모형·노출 / 한 KM endpoint의 전체 쌍입니다. 연구 전체를 한꺼번에 보정한 값이 아닙니다. BH q는 기존 omnibus/모형 계수 family 보정으로 별개입니다. KM에는 Holm p<0.05인 쌍을 모두 표시하고 없으면 없다고 적습니다. 박스플롯은 쌍별 주석 또는 표 파일명을 제공합니다.

Cox HR과 RMST 차이는 group1 대 group2 방향입니다. CI는 개별 95% CI이며 동시 구간이 아닙니다. Fisher/CMH OR은 outcome_first_level 대 outcome_second_level의 odds를 group1/group2로 비교합니다. Fisher에는 chi-square 통계량이 없으므로 statistic/df가 NA인 것이 정상입니다.

### 검증·보고서 재생성

```sh
Rscript --vanilla tests/test_pairwise.R
Rscript --vanilla tests/verify_pairwise_outputs.R --config=config/local.R
python3 scripts/build_results_report.py --results outputs/updated_v3_20260928
```

Dunn 구현은 rstatix::dunn_test와 대조하므로 해당 검증에만 rstatix가 필요합니다. 보고서는 임상 특성→유전자/병리/OS→KRAS 두 출처→TMB/Variant count→재발·DNA 복구 순입니다. Table1–7은 현재 집계 자료로 새로 작성한 Word 편집 가능 표이며 같은 표시값의 Table*_Manuscript_*.tsv도 저장합니다. 30_Supplementary_table_index.tsv는 최초 연결된 Figure 번호 기준으로 중복 파일에 같은 보조표 번호를 유지합니다. 30_Figure_inventory.tsv에는 본문 외 탐색 그림까지 포함합니다. 환자별 private 자료는 공개 보조표 목록에 넣지 않습니다.

## HRD 관련·MMR oncoplot 및 임상 분포

McIntyre 2020 Methods(p3940)의 HRD 관련 목록 18개와 MMR 4개를 사용합니다. 정확한 목록과 출처는 `31_Repair_gene_list_definition.tsv`에 기록합니다. FAM175A의 현행 별칭 ABRAXAS1을 연결하며 다른 RAD51 계열 유전자를 임의로 추가하지 않습니다.

- 31: 전체 환자 및 각 유전자군 변이 보유 환자만의 oncoplot 4개, 유전자별 환자 수·MAF 행 수·변이 분류, 두 목록 중복 분포. oncoplot 행의 백분율은 **해당 그림에 표시된 환자 수**가 분모입니다. 집계 빈도 TSV는 전체 주 코호트가 분모입니다.
- 32: 각 유전자군 변이 검출/미검출 및 4개 중복 조합군의 임상 분포, 전체 검정과 모든 쌍의 Holm 보정, 개별 22개 유전자별 임상 기술통계. 개별 희귀 유전자의 검정 결과를 선별해 제시하지 않습니다.
- 연속형은 2군 Wilcoxon/다군 Kruskal–Wallis, 범주형은 Fisher(2×2 exact, 그 외 10,000회 Monte Carlo). 전체 검정 BH는 grouping별 임상 항목 family, 쌍별 Holm은 grouping×항목별 모든 쌍입니다.
- 기술통계의 범주 백분율은 결측을 포함한 군 전체 N, 검정은 해당 항목 complete cases입니다. 연속형은 관측값의 중앙값[IQR] 및 결측 N을 제공합니다.
- 유지된 MAF 변이를 사용하며 병원성/기능적 결핍을 판정하지 않습니다. 0개는 패널 검출 범위가 확인된 진정한 음성이 아닙니다. 유전자군 변이 보유와 Variant count의 연관은 정의상·검출 기회상 연관이 포함됩니다.
- LOH/CNV, tumor purity, tumor cellularity는 분석·모형·요약 비교·그림 주석에서 제외합니다. 입력 감사표와 원본 자료는 보존하므로 원래 컬럼이 남아 있을 수 있습니다.
- 보고서 Fig21–25와 Table6–7은 이 추가 결과를 사용합니다. 본문과 편집 가능한 표의 글꼴은 Pretendard입니다. 원본 PNG는 축소하지 않은 별도 논문용 검토 자료입니다.

```sh
Rscript --vanilla scripts/repair_gene_profiles.R --config=config/local.R
Rscript --vanilla tests/verify_repair_outputs.R --config=config/local.R
Rscript --vanilla tests/verify_pairwise_outputs.R --config=config/local.R
python3 scripts/build_results_report.py --results outputs/updated_v3_20260928
```
