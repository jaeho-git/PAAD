# 분석 스크립트 안내

## 현재 논문용 실행 경로

새 임상자료에서는 `run_all.R` → `manuscript_analysis.R`이 기본입니다. 현재 기준 입력은 `data/raw/260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx`입니다. 메인 Figure 5개와 Table 2개, 선택된 보충 Figure를 논문 순서로 생성합니다. Supplementary Table 시리즈는 만들지 않고 그림·표의 계산 근거를 `Input_*.tsv`로 저장합니다.

- `manuscript_analysis.R`: 단일 KRAS 출처, 임상 Reported TMB, 병리 및 생존분석, HRD/MMR oncoplot의 기본 분석.
- `R/manuscript_story_extensions.R`: 선행 platinum, 재발양상, 명확한 병리 비교 그림, 메인 Table 1–2를 추가하고 5개 Figure 논문 흐름을 완성.
- `R/recurrence_signature_extensions.R`: Distant only를 Liver only·Lung only·Other distant로 확장한 재발 후 생존분석, 5군 재발 분석, 선행 platinum별 HRD/MMR retained-variant 비교를 추가.
- `mutational_signature_analysis.py`: MAF의 SNV 문맥을 점검하고 SigProfiler 기반 환자별 품질평가, 요청된 positive rate·burden·proportion의 산출 가능성, 통합 process와 개별 SBS의 집단 합산 민감도 분석, standard assignment 대 constrained refit 감사를 수행. 품질기준 미달 환자를 음성으로 바꾸지 않으며 환자별 결과는 `private/`에만 저장.
- `build_updated_manuscript_report.py`: 집계 결과와 그림을 이용해 상세 해석이 포함된 Pretendard DOCX를 생성.
- `build_updated_manuscript_slides.mjs`: 확정된 5개 Figure 흐름으로 Nature-inspired PPTX를 생성.
- `R/oncoplot_counts.R`: 환자별 변이, 유전자별 변이, 변이 보유 환자를 구분하는 집계 함수.

설정, 입력·결과 파일과 문서 생성 명령은 [논문용 분석 안내](../docs/manuscript_analysis.md)에 있습니다. 아래 기존 파일들은 삭제하지 않고 재현 이력으로 보존했습니다. 특히 `source_comparisons.R`은 내부 QC이며 논문 또는 보충자료에 포함하지 않습니다. 이전 `build_results_report.py`는 새 논문용 보고서 생성기가 아닙니다.

## 기존 분석 이력

이 디렉터리에는 기존 분석 파일 7개와 새 임상자료 감사·문헌·출처·HRD/MMR 확장 파일 4개가 있습니다. 파일명은 분석 내용을 나타내며, 각 파일 안에 해당 분석의 함수와 실행 순서가 함께 있습니다.

스크립트 이름의 순번은 제거했지만 생성 결과의 앞 번호는 기존 figure·표와의 대응을 보존하기 위해 유지했습니다.

공통 파일인 `R/config.R`, `R/data.R`, `R/plot_helpers.R`에는 설정 검증, 입력 전처리, 안전한 출력 경로와 반복되는 그림 보조 기능만 둡니다. 통계 검정, 표 구성, 유전자 선택, survival group 정의와 실제 출력 호출은 아래 분석 파일에서 확인할 수 있습니다.

## 빠르게 찾기

| 알고 싶은 분석 | 먼저 열 파일 | 생성 결과 |
| --- | --- | --- |
| 새/구 임상 변수·식별자·대상자 비교 | `clinical_update_audit.R` | Output 00: 비공개 비교표 |
| 전체 및 T/N 정보 보유군 oncoplot | `oncoplots.R` | Output 1, 1-2 |
| 임상 특성 요약표 | `clinical_summary.R` | Output 2 |
| 임상 그룹별 변이빈도·변이개수 비교 | `clinical_group_comparisons.R` | Output 4 |
| Top mutation 및 변이개수 heatmap | `mutation_heatmaps.R` | Output 5–7 |
| 유전자 간 somatic interaction | `somatic_interactions.R` | Output 9 |
| Stage별 환자 수 | `stage_distribution.R` | Output 10 |
| KRAS·driver·clinicopathology·생존 | `driver_kras_survival.R` | Output 11–17 |
| 논문 참고 확장표, 보정 Cox/CMH, 생존·재발 양상, KRAS 교차검증 | `literature_extensions.R` | Output 18–24 |
| 임상/MAF KRAS 별도 분석·비교, TMB paired/결측 비교, 신규 임상변수 분포, 출처별 생존모형 | `source_comparisons.R` | Output 25–29 |
| HRD 관련·MMR oncoplot, 유전자별·유전자군별 임상 분포 및 쌍별 비교 | `repair_gene_profiles.R` | Output 31–32 |
| Pretendard Word 결과보고서와 원고용 표 | `build_results_report.py`, `report_interpretations.py` | Fig1–25, Table1–7, Output 30 |

원본 v19 실행부에는 Output 3과 Output 8이 없습니다. 존재하지 않는 분석을 추정해 추가하지 않았습니다.

## 실행 방법

모든 명령은 저장소 루트에서 실행합니다. 실제 로컬 입력은 Git에서 제외되는 `config/local.R`을 사용합니다.

```sh
Rscript --vanilla scripts/oncoplots.R --config=config/local.R
Rscript --vanilla scripts/clinical_summary.R --config=config/local.R
Rscript --vanilla scripts/clinical_group_comparisons.R --config=config/local.R
Rscript --vanilla scripts/mutation_heatmaps.R --config=config/local.R
Rscript --vanilla scripts/somatic_interactions.R --config=config/local.R
Rscript --vanilla scripts/stage_distribution.R --config=config/local.R
Rscript --vanilla scripts/driver_kras_survival.R --config=config/local.R
```

전체 분석은 같은 순서로 한 번에 실행할 수 있습니다.

```sh
Rscript --vanilla run_all.R --config=config/local.R
```

저장소에 포함된 완전 합성 입력으로 실행 흐름만 확인하려면 설정을 바꿉니다.

```sh
Rscript --vanilla scripts/clinical_summary.R --config=config/config.synthetic.R
Rscript --vanilla run_all.R --config=config/config.synthetic.R
```

각 분석 파일은 새 R session에서 독립적으로 실행할 수 있으며, 다른 분석 파일이 만든 전역 객체나 중간 결과 파일을 요구하지 않습니다. 결과는 설정의 `output_dir`에 생성됩니다. `clinical_update_audit.R`, `literature_extensions.R`, `source_comparisons.R`, `repair_gene_profiles.R`은 `curated_v3`에서만 실행하고 legacy 설정에서는 자동 생략합니다. 합성 설정은 기존 42개 결과 계약을 유지합니다. [새 분석 방법](../docs/updated_analysis.md)을 먼저 확인하십시오.

현재 논문용 재분석은 `data/raw/260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx`의 `Main` 시트와 `data/raw/PDAC_oncopanel.maf`를 사용하며 결과는 `outputs/updated_v3_20260928/manuscript_extension_20261005_kras_maf/`에 저장합니다. 임상 원자료 1,015행 중 중복 ID 2개에 해당하는 4행을 제외해 M0 972명과 M1 39명, 총 1,011명을 포함합니다.

## 코드를 읽는 순서

개별 파일은 위에서 아래로 다음 흐름을 따릅니다.

1. 저장소 위치와 `--config` 인자를 확인합니다.
2. 공통 설정·입력·그림 보조 파일을 불러옵니다.
3. 그 분석에만 필요한 함수와 통계 규칙을 정의합니다.
4. 설정과 세 입력 파일을 읽습니다.
5. 분석을 실행하고 이름이 고정된 결과를 `output_dir`에 저장합니다.

정확한 원본 코드 구간, 입력 조건, 출력 파일명과 검증 범위는 [`docs/analysis_map.md`](../docs/analysis_map.md)를 참조하십시오. 합성 예제는 소프트웨어 실행 확인용이며 과학적 결과 또는 통계적 타당성을 검증하는 자료가 아닙니다.

`R/genomics.R`에는 여러 그림에서 동일하게 쓰는 KRAS 변이 판정 규칙만 있습니다. 실제 그림·검정·저장 코드는 위 스크립트에 있습니다. `source_comparisons.R`도 원본 입력을 직접 읽으므로 다른 스크립트 결과를 먼저 만들 필요가 없습니다.

## 쌍별 통계 및 보고서

R/pairwise_tests.R는 log-rank, Dunn, Fisher/CMH, Cox contrast의 작은 통계 함수 모음입니다. 실제 비교 변수 선택·그림·저장은 각 스크립트에 그대로 있습니다. curated_v3 실행에서 pairwise TSV가 추가되며 legacy 합성 분석의 42개 결과 계약은 유지됩니다.

build_results_report.py는 완료된 R 결과를 읽어 Word 보고서를 만듭니다. run_all.R 성공 후 실행하십시오. Python 패키지: python-docx, pandas, Pillow.

```sh
python3 scripts/build_results_report.py --results outputs/updated_v3_20260928
```

NE는 추정 불가, Not in model은 해당 모형에 미포함을 뜻합니다.

HRD/MMR 분석만 갱신하려면 `Rscript --vanilla scripts/repair_gene_profiles.R --config=config/local.R`을 실행합니다. 유전자 목록/별칭은 `R/repair_gene_sets.R`, 그림의 명확한 표시명은 `R/data.R::figure_label()`에 있습니다. 원자료/매핑 컬럼명과 그림용 표시명은 구분합니다. 제외 항목(LOH, CNV, purity, cellularity)은 주석에도 넣지 않습니다.

보고서 생성 후에는 Word/PDF 렌더링으로 페이지 배치를 확인해야 합니다. 보고서는 원고 구성·검토용이며 결과 해석과 연구 설계의 한계 검토를 대체하지 않습니다.
