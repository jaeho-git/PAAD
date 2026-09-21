# 분석 대응표

이 문서는 기준 통합 스크립트 `PDAC_plot_n527_modified_v19_publication_update.R`의 분석을 새 실행 파일과 연결한다. 원본의 줄 번호는 조사한 v19 파일(SHA-256 `c88ceb7476c6a7767f30f6036ed811c74b7fb50d4f4393cd53ef94b87c98b073`) 기준이다. 파일명의 `n527`은 원본 대상 목록을 가리킬 뿐이며, 아래 분석의 최종 환자 또는 검체 수가 527임을 뜻하지 않는다.

모든 명령은 저장소 루트에서 실행한다. 실제 자료는 `data/raw/`에만 로컬로 복제되며 Git에서 제외된다. `config/local.R`을 사용하는 실행은 그 실제 자료를 읽고, 결과를 설정된 `outputs/` 하위 경로에 새로 만든다. 기존 결과 폴더의 파일은 복사하거나 덮어쓰지 않는다.

| 분석명 | 원본 파일·구간 | 새 실행 파일 | 입력·선행 처리 | 출력 파일/폴더 | 실행 명령 | 검증 상태 |
| --- | --- | --- | --- | --- | --- | --- |
| 공통 입력 및 PDAC cohort 구성 | v19 19–163행: 경로, XLSX/CP949 TSV/MAF 입력, join, 범주 정리, `Classification` 필터 | `R/config.R`, `R/data.R` | 임상 XLSX, 대상 목록 TXT, MAF; barcode join 후 `Exception` 제외 및 `PDAC` 문자열 필터 | 독립 파일 없음; 이후 분석에서 사용하는 메모리 객체 | 각 실행 파일에서 공통 호출 | 실제 입력 load 완료: filter 후 임상 509행, retained MAF 4,224행 |
| Main 및 T/N 유효 표본 oncoplot (Output 1, 1-2) | 166–305행의 색상/oncoprint 함수, 1182–1204행 실행부 | `scripts/oncoplots.R` | 공통 입력; 두 번째 그림은 T와 N이 모두 결측이 아닌 barcode | `1_Oncoplot_Main_PDAC.png`, `1_2_Oncoplot_TN_Filtered_PDAC.png` | `Rscript --vanilla scripts/oncoplots.R --config=config/local.R` | 실제 자료 실행 완료(top 5 검증 설정)·figure 동등성 미확인 |
| 임상 요약표 (Output 2) | 1206–1217행 | `scripts/clinical_summary.R` | 공통 입력; annotation 변수 | `2_Clinical_Table_PDAC.xlsx` | `Rscript --vanilla scripts/clinical_summary.R --config=config/local.R` | 실제 자료 실행 및 기존 12개 sheet 값 동등성 확인 |
| 변이 빈도 및 원본 명칭의 “TMB” 비교 (Output 4) | 345–475행, 1219–1222행 | `scripts/clinical_group_comparisons.R` | 공통 입력; `Differentiation`, `NAC`, `T`, `N`, `N_status`, `Stage`, `Stage_Group` 중 존재하는 변수 | 변수별 `4_FreqPlot_PDAC_<변수>.png`, `4_TMB_BoxPlot_PDAC_<변수>.png` | `Rscript --vanilla scripts/clinical_group_comparisons.R --config=config/local.R` | 실제 자료 실행 완료(top 5 검증 설정)·수치/figure 동등성 미확인 |
| 변이 binary heatmap Top 20/5 및 변이 개수 heatmap (Output 5–7) | 308–342행, 478–514행, 1224–1249행 | `scripts/mutation_heatmaps.R` | 공통 입력; MAF에서 빈도 상위 유전자와 sample별 retained MAF row count | `5_Heatmap_Top20_PDAC.png`, `6_Heatmap_Top5_PDAC.png`, `7_TMB_Heatmap_PDAC.png` | `Rscript --vanilla scripts/mutation_heatmaps.R --config=config/local.R` | 실제 자료 실행 완료(top 5 검증 설정)·figure 동등성 미확인 |
| Somatic interaction (Output 9) | 1251–1257행 | `scripts/somatic_interactions.R` | 공통 입력; `maftools::somaticInteractions`, 상위 25 유전자 | `9_Somatic_Interactions_PDAC.png`, `9_Somatic_Interactions_Results_PDAC.xlsx` | `Rscript --vanilla scripts/somatic_interactions.R --config=config/local.R` | 실제 자료 실행 완료·수치/figure 동등성 미확인 |
| Stage별 수 (Output 10) | 1259–1274행 | `scripts/stage_distribution.R` | 공통 입력; 결측 Stage 제외, 코드에 정의된 level 유지 | `10_Stage_Counts_PDAC.png` | `Rscript --vanilla scripts/stage_distribution.R --config=config/local.R` | 실제 자료 실행 완료·figure 동등성 미확인 |
| KRAS subtype, 4개 driver 요약, clinicopathology association, OS/RFS (Output 11–17) | 518–1171행, 1276–1277행 | `scripts/driver_kras_survival.R` | 공통 입력; `KRAS`, `TP53`, `SMAD4`, `CDKN2A`; OS는 `survive`, RFS는 `Recur`를 event로 사용 | `11_Driver_KRAS_patient_level_summary.xlsx`, `11_KRAS_subtype_pie.png`, `12_Driver_mutation_status_pie.png`, `13_Driver_gene_mutation_clinicopath_table.{xlsx,tsv}`, `13_Driver_gene_mutation_clinicopath_table_display.xlsx`, `14_17_PDAC_survival_input_patient_level.xlsx`, `14_*`, `15_*`, `15_2_*`, `16_*`, `17_*`, `17_2_*` PNG와 `*_number_at_risk.tsv` | `Rscript --vanilla scripts/driver_kras_survival.R --config=config/local.R` | 실제 자료 실행·일부 표/TSV 동등성 확인; simulated p-value와 figure는 미확인 |
| 전체 분석 | v19 1173–1279행 실행부 전체 | `run_all.R` | 위 공통 입력과 일곱 분석 파일 | 설정된 `output_dir` 바로 아래 위 파일 전체 | `Rscript --vanilla run_all.R --config=config/local.R` | 합성 전체 실행 완료(42개 파일); 실제 production 설정 전체 실행은 미수행 |

## 실행 순서와 독립성

각 `scripts/*.R`은 새 R session에서 설정과 입력을 직접 준비한다. 다른 분석 스크립트가 만든 전역 객체나 중간 파일을 선행 조건으로 사용하지 않는다. 각 분석의 함수와 실행 호출은 해당 파일에 함께 있어, 그 파일 하나에서 분석 흐름을 읽을 수 있다. `R/config.R`, `R/data.R`, `R/plot_helpers.R`에는 공통 기반 기능만 둔다.

`run_all.R`은 아래와 같은 명시적 순서로 일곱 분석 파일을 실행한다. 전체 실행을 위한 별도의 분석 구현은 두지 않으므로 개별 실행과 전체 실행이 같은 코드를 사용한다.

1. `oncoplots.R`
2. `clinical_summary.R`
3. `clinical_group_comparisons.R`
4. `mutation_heatmaps.R`
5. `somatic_interactions.R`
6. `stage_distribution.R`
7. `driver_kras_survival.R`

합성 예제 실행은 `config/config.synthetic.R`을 사용한다. `config/config.example.R`은 `data/raw/`의 실제 로컬 입력을 가리키는 `config/local.R` 작성 예시다.

```sh
Rscript --vanilla scripts/clinical_summary.R --config=config/config.synthetic.R
Rscript --vanilla run_all.R --config=config/config.synthetic.R
```

합성 결과는 입출력 연결과 실행 가능성을 확인하는 용도이며 실제 연구 결과의 재현 근거가 아니다.

실제 자료 검증 산출물은 Git 제외 대상 `outputs/code_validation/`에 격리했다. 속도와 부작용을 제한하기 위해 oncoplot, comparison, heatmap은 `top_n=5` 검증 설정을 사용했으며 production 기본값인 top 20·1000 dpi 전체 실행은 하지 않았다. 기존 결과와 확인한 동등성은 임상 요약 workbook 12개 sheet, 509×35 driver/KRAS 및 survival input 표, 세 개 risk TSV, driver clinicopathology 표의 p-value 이외 필드다. RNG 조건을 재현하지 못한 simulated Fisher p-value와 모든 PNG의 pixel/figure 동등성은 확인하지 않았다.

## 원본 번호에 없는 항목

v19 실행부에는 Output 3과 Output 8이 없다. 존재하지 않는 분석을 추정해 새 placeholder를 만들지 않았다. `01_Results_copy` 등 다른 시점의 결과 이름은 v19 실행부의 근거가 아니므로 이 표에 통합하지 않았다.
