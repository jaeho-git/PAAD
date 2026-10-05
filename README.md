# PAAD 분석 코드

PAAD 환자의 targeted exome sequencing 자료와 임상정보를 이용한 PDAC 분석 코드를 정리한 저장소입니다. 기준 통합 스크립트 `PDAC_plot_n527_modified_v19_publication_update.R`의 분석 의미와 출력 형식을 유지하면서, 공동연구자가 분석 내용을 파일명과 코드만 보고 따라갈 수 있도록 분석별로 나누었습니다.

`scripts/`의 각 파일은 단순한 실행용 wrapper가 아닙니다. 해당 분석의 함수, 통계 규칙, 표·그림 생성 코드와 실제 실행 순서가 한 파일 안에 함께 있습니다. 따라서 분석을 이해하기 위해 별도의 초기화 파일이나 실행 중계 파일을 찾아갈 필요가 없습니다.

이 저장소는 공개 가능한 코드·문서·합성 예제만 Git으로 관리합니다. 실제 환자자료와 새 분석 결과는 로컬에만 보관합니다.

## 업데이트 임상파일 분석

`Main`과 `Data_dictionary` 시트가 있는 새 임상파일은 `clinical_schema = "curated_v3"`로 읽습니다. 이 모드에서는 새 파일의 임상·병리·생존 값을 사용하고, 예전 대상 TXT는 비교 감사에만 사용합니다. 설정 예시는 `config/config.curated.example.R`입니다. **기존 `config/local.R`이 있으면 덮어쓰지 말고 필요한 항목을 확인하십시오.**

새 임상자료의 기본 실행은 **논문용 분석**입니다. `run_all.R`이 `scripts/manuscript_analysis.R`을 실행하여 메인 Figure 1–5, Table 1–2와 선택된 보충자료를 만듭니다. M0·M1을 모두 포함하고 중복 환자 ID의 모든 행만 제외합니다. 현재 로컬 분석은 `data/raw/260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx`의 `Main` 시트를 사용하며, 1,015행 중 중복 환자 ID에 해당하는 4행을 제외한 1,011명을 분석합니다. 임상/MAF 출처 간 일치도 비교는 이 경로에서 생성하지 않습니다. 현재 논문용 KRAS 아형은 엑셀의 `KRAS_subtype_MAF`를 직접 사용하도록 `config/local.R`의 `kras_source = "workbook_maf"`로 고정했습니다. `"clinical"` 또는 원시 MAF 재계산값인 `"maf"`도 검증용 선택지로 남아 있지만 한 실행에서 서로 섞이지 않습니다.

TMB는 임상정보의 **Reported TMB**를 탐색적으로 사용합니다. 단위·검사법을 확인하지 못한 상태에서는 mut/Mb 또는 TMB-high로 표시하지 않습니다. Variant count는 oncoplot의 기술적 주석이며 TMB 결측을 대체하지 않습니다. 기록된 재발은 사망 포함 DFS/RFS와 구분합니다. 새 결과 구성·재생성 방법은 [논문용 분석 및 문서 안내](docs/manuscript_analysis.md)를 먼저 읽으십시오. [이전 확장 분석 안내](docs/updated_analysis.md)는 내부 검토·기존 결과의 이력입니다.

```sh
Rscript --vanilla run_all.R --config=config/local.R
```

논문용 결과는 `output_dir/manuscript_subdir/`에 저장됩니다. 기본 하위 폴더명은 `manuscript`입니다. 과거 legacy 분석은 `output_dir`에 직접 저장됩니다. 실제 데이터·환자별 자료·산출물은 계속 Git에서 제외됩니다.

## 핵심 저장 원칙

- 원본 프로젝트는 읽기 전용으로 취급합니다.
- 원본 `00_Data`의 파일은 로컬 실행용 사본을 `data/raw/`에 두며, 이 경로 전체를 Git에서 제외합니다.
- `data/example/`에는 실제 record를 추출·익명화한 자료가 아닌 완전 합성 예제만 포함합니다.
- 기존 `01_Results`와 `01_Results_copy`는 이 저장소로 복사하거나 백업하지 않습니다.
- `outputs/`는 새 코드를 실행해서 생성되는 figure, 표와 로그만 저장하며 Git에서 제외합니다.
- 개인 설정은 `config/local.R`에 두고 Git에서 제외합니다.
- `.gitignore`는 데이터 확장자를 일괄 제외하지 않고 실제 데이터·결과 경로를 제외하므로, 합성 XLSX·TXT·MAF는 Git으로 관리할 수 있습니다.

## 저장소 구조

```text
PAAD/
├── R/                         # 설정·입력 전처리·공통 그림 보조 기능만 포함
├── scripts/                   # manuscript_analysis.R + 문서 생성기 + 기존 분석 이력
├── config/
│   ├── config.example.R       # 실제 로컬 입력용 공개 설정 예시
│   ├── config.curated.example.R # 새 임상파일용 설정 예시
│   ├── config.synthetic.R     # 합성 예제 smoke test 설정
│   └── local.R                # 현재 컴퓨터의 실제 실행 설정; Git 제외
├── data/
│   ├── README.md              # 입력 형식과 합성 자료 설명
│   ├── example/               # 공개 가능한 합성 입력과 generator
│   └── raw/                   # 실제 입력 사본; Git 제외
├── docs/
│   ├── analysis_map.md        # 원본 코드·새 실행 파일·출력 대응표
│   ├── data_dictionary.md     # 입력 schema와 파생 변수
│   └── refactoring_notes.md   # 변경, 검증 범위와 판단 보류 사항
├── outputs/                   # 새 실행 결과; Git 제외
├── local_notes/               # 비공개 조사·검증 기록; Git 제외
├── tests/run_tests.R          # 합성 자료 기반 통합 점검
├── run_all.R                  # curated_v3: 논문용 분석 / legacy: 기존 분석
└── PAAD.Rproj
```

## 실행환경

검증한 환경은 R 4.5.3입니다. 분석에 필요한 package는 다음과 같습니다.

- CRAN: `readxl`, `writexl`, `tidyverse`, `circlize`, `RColorBrewer`, `survival`, `cowplot`, `ragg`, `jsonlite`, `png`, `broom`
- Bioconductor: `maftools`, `ComplexHeatmap`

스크립트는 package를 자동 설치하거나 사용자의 global library를 변경하지 않습니다. 누락된 package가 있으면 실행을 중단하고 이름을 표시합니다. 새 환경에서는 R console에서 다음과 같이 준비할 수 있습니다.

```r
install.packages(c(
  "readxl", "writexl", "tidyverse",
  "circlize", "RColorBrewer", "survival", "cowplot", "ragg", "jsonlite", "png", "broom"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("maftools", "ComplexHeatmap"))
```

검증 시 사용한 주요 버전은 `maftools 2.26.0`, `readxl 1.4.5`, `writexl 1.5.4`, `tidyverse 2.0.0`, `ComplexHeatmap 2.26.1`, `circlize 0.4.18`, `RColorBrewer 1.1-3`, `survival 3.8-6`입니다. 별도 `renv.lock`을 이용한 다른 컴퓨터에서의 재현성은 아직 검증하지 않았습니다.

## 먼저 합성 예제로 확인하기

저장소 루트에서 실행합니다.

```sh
Rscript --vanilla tests/run_tests.R
```

이 검사는 다음을 포함합니다.

- 필수 파일과 R 문법 확인
- `.gitignore` 및 출력 경로 안전장치 확인
- 합성 XLSX, CP949·CRLF TXT와 MAF 형식 확인
- PDAC cohort 구성과 driver 요약 확인
- 합성 입력을 이용한 7개 분석의 개별 실행과 전체 실행
- 실행별로 격리된 `outputs/test_runtime_<process-id>/`에 예상 산출물 42개가 생성되는지 확인한 뒤 성공 시 테스트 산출물 정리

전체 합성 분석만 별도로 실행하려면 다음 명령을 사용합니다.

```sh
Rscript --vanilla run_all.R --config=config/config.synthetic.R
```

합성 입력은 24개 가상 record로 구성되며 v19 filter 후 20개 PDAC record와 300개 MAF row가 분석에 남습니다. 모든 식별자는 `SYN-` 접두사를 사용합니다. 이는 파일 연결과 실행 가능성을 확인하는 smoke test이며 실제 연구 결과의 재현 또는 통계적 검증이 아닙니다.

합성 입력을 결정적으로 다시 만들거나 현재 파일을 점검할 수도 있습니다.

```sh
python3 data/example/generate_examples.py
python3 data/example/generate_examples.py --check
```

자세한 형식은 [data/README.md](data/README.md)를 참조합니다.

## 실제 로컬 데이터 준비

새 clone에서는 원본 입력의 로컬 사본을 다음 위치에 준비합니다.

```text
data/raw/PDAC_clinical_info.xlsx
data/raw/result_527명.txt
data/raw/PDAC_oncopanel.maf
```

`분석용데이터1_수작업_eOJH_sdk0116.xlsx`도 현재 로컬 `data/raw/`에 복제되어 있지만, 기준 v19 분석에서는 사용하지 않습니다.

공개 설정을 로컬 설정으로 복사한 뒤 필요한 경로만 수정합니다.

```sh
cp config/config.example.R config/local.R
```

현재 작업 환경에는 위 실제 입력 사본과 `config/local.R`이 이미 준비되어 있습니다. 두 항목 모두 Git에서 제외됩니다.

## 전체 분석 실행

```sh
Rscript --vanilla run_all.R --config=config/local.R
```

논문용 분석 결과는 `output_dir/manuscript_subdir/`에 생성됩니다. 새 임상 설정 예시에서는 `outputs/updated_v3_run/manuscript/`입니다. 현재 로컬 설정의 재분석 결과는 `outputs/updated_v3_20260928/manuscript_extension_20261005_kras_maf/`에 있습니다. 설정 검증은 출력 위치가 이 저장소의 `outputs/` 또는 그 하위인지 확인합니다. `run_all.R`은 R 분석까지만 수행하며 Word·PowerPoint 생성은 아래 연결 문서의 별도 명령을 사용합니다.

재분석 결과를 이전 실행과 명확히 구분하려면 `config/local.R`의 `output_dir`을 실행별 새 하위 폴더로 지정하십시오. 같은 폴더를 다시 사용하면 같은 이름의 파일은 갱신되지만, 입력 조건 부족으로 이번 실행에서 생략된 분석의 과거 파일은 자동 삭제하지 않습니다.

## 분석별 실행

각 실행 파일은 새로운 R session에서 독립적으로 입력과 설정을 준비합니다. 파일을 열면 상단의 목적·입력·출력 설명에 이어 분석 함수와 명시적인 실행부를 같은 곳에서 볼 수 있습니다. 전체 파일별 안내는 [scripts/README.md](scripts/README.md)를 참조하십시오.

스크립트 파일명에서는 `01_`, `02_` 같은 순번을 없앴습니다. 결과 파일의 앞 번호는 기존 figure·표와의 대응을 유지하기 위해 그대로 두었습니다.

| 분석 | 실행 파일 | 주요 출력 |
| --- | --- | --- |
| Main 및 T/N-filtered oncoplot | `scripts/oncoplots.R` | `1_*.png` |
| 임상 요약표 | `scripts/clinical_summary.R` | `2_Clinical_Table_PDAC.xlsx` |
| 임상 그룹별 변이빈도·변이개수 비교 | `scripts/clinical_group_comparisons.R` | `4_FreqPlot_*.png`, `4_TMB_BoxPlot_*.png` |
| Top 20/5 mutation 및 변이개수 heatmap | `scripts/mutation_heatmaps.R` | `5_*.png`, `6_*.png`, `7_*.png` |
| Somatic interaction | `scripts/somatic_interactions.R` | `9_*.png`, `9_*.xlsx` |
| Stage count | `scripts/stage_distribution.R` | `10_Stage_Counts_PDAC.png` |
| KRAS·driver·clinicopathology·survival | `scripts/driver_kras_survival.R` | `11_*`–`17_*` |

예를 들어 실제 로컬 데이터로 임상 요약표만 생성하려면 다음과 같이 실행합니다.

```sh
Rscript --vanilla scripts/clinical_summary.R --config=config/local.R
```

원본 구간, 정확한 출력 파일과 현재 검증 상태는 [docs/analysis_map.md](docs/analysis_map.md)에 정리되어 있습니다.

## 주요 설정

`config/local.R`과 `config/config.synthetic.R`에서 다음 항목을 관리합니다.

| 항목 | 의미 |
| --- | --- |
| `clinical_file`, `target_patients_file`, `maf_file` | 세 입력 파일 |
| `output_dir` | 반드시 저장소의 `outputs/` 하위 |
| `target_encoding` | cohort TXT 인코딩; 현재 `CP949` |
| `clinical_schema` | `legacy` 또는 `curated_v3`; 입력 값의 출처와 매핑 선택 |
| `clinical_sheet` | 새 파일은 `Main`; `Data_dictionary`도 필수 |
| `updated_cohort` | 새 입력 분석 대상: 기본 `all_unique` (M0+M1, 중복 ID 행 전체 제외). `localized_unique`는 별도 요청 시에만 쓰는 과거 M0 제한 설정 |
| `sex_column_policy` | `legacy`, `target_patients`, `clinical_korean` 중 선택 |
| `random_seed` | `NULL`이면 v19의 비결정적 simulated Fisher 동작 유지; 정수이면 재현 가능 |
| `top_n` | 기본 상위 유전자 수 20; 낮은 값은 빠른 코드 검증용이며 기존 figure 제목·파일명은 유지 |
| `plot_dpi` | production 기본값 1000; smoke test는 150 |

기존 `legacy` 입력에서는 `sex_column_policy`로 성별의 출처를 선택합니다. `curated_v3`에서는 사전에 정의된 새 `Sex`의 F/M 값을 직접 사용하며 이 옵션은 적용하지 않습니다.

## 기존 v19 리팩터링 검증 기록

- 기준 v19 스크립트 전체 1,279행의 입력·분석·출력을 조사했습니다.
- 합성 입력으로 7개 분석 파일의 개별 실행과 `run_all.R` 전체 실행이 각각 완료되어 예상 산출물 42개를 확인했습니다.
- 실제 로컬 입력은 v19 filter 후 임상 509행, retained MAF 4,224행으로 로드됐습니다.
- 실제 자료로 각 분석을 격리된 `outputs/code_validation/`에서 실행했습니다. 계산량이 큰 그림 분석은 top 5 검증 설정을 사용했습니다.
- 기존 결과와 임상 요약 workbook 12개 sheet, 509×35 driver/KRAS 표, 509×35 survival input 표, risk TSV 3개가 일치했습니다.
- Driver-clinicopathology 표는 구조와 p-value 이외 필드가 일치했습니다.
- simulated Fisher p-value의 기존 RNG 상태와 PNG의 pixel/figure 동등성은 확인하지 않았습니다.
- production 설정의 top 20·1000 DPI 전체 실제 데이터 `run_all.R`은 실행하지 않았습니다.

상세 변경과 한계는 [docs/refactoring_notes.md](docs/refactoring_notes.md)를 참조합니다.

2026-09-28 업데이트에서는 새 schema 단위 검사와 기존 합성 회귀검사를 통과하고, top 20·1000 DPI 설정으로 실제 `run_all.R` 전체를 실행했습니다. 확장 그림은 파일 크기와 가독성을 위해 최대 400 DPI로 저장합니다. 구체적인 환자 수·비교 결과는 Git 제외 결과 폴더의 감사표에 저장합니다.

## 해석 시 주의사항

- 기존 `TMB` 파일명은 유지하지만 축·제목은 Variant count로 표시합니다. 새 실제 TMB(mutations/Mb)는 별도로 분석합니다.
- 일부 기존 omnibus Fisher p는 nominal 값입니다. 새 curated 분석에는 전 그룹 쌍별 검정과 Holm 보정을 추가했으며, 전체 연관성 표의 BH q와 구분합니다.
- legacy 입력의 event 정의는 기존 사전을 확인해야 합니다. 새 입력에서는 y/n을 사망/재발 event 1/0으로 명시적으로 변환하며, 재발은 사망 포함 DFS/RFS와 구분합니다.
- KRAS subtype의 좌표 fallback은 MAF의 genome build와 일치하는지 확인해야 합니다.
- 환자와 검체 단위, 중복 barcode 대표행, staging edition과 biomarker 검출한계 처리에는 연구자 판단이 필요합니다.

확인된 입력 schema와 파생 변수는 [docs/data_dictionary.md](docs/data_dictionary.md)에, 전체 판단 보류 목록은 [docs/refactoring_notes.md](docs/refactoring_notes.md)에 기록했습니다.

## Git 공개 전 확인

`data/raw/`, `config/local.R`, `outputs/`, `local_notes/`는 기본적으로 Git에서 제외됩니다. 합성 예제는 제외되지 않습니다.

```sh
git check-ignore -v data/raw/PDAC_clinical_info.xlsx
git check-ignore -v config/local.R
git check-ignore -v outputs/local_run/1_Oncoplot_Main_PDAC.png
git status --short
```

실제 자료, 생성 결과 또는 개인 경로가 Git 변경 목록에 나타나지 않는지 확인한 뒤 사용자가 staging·commit·push 여부를 결정합니다.

## 쌍별 비교와 원고용 보고서

curated_v3 분석의 다군 비교는 전체 검정 외에 모든 쌍의 통계량·원 p·Holm p를 TSV로 저장합니다. KM에는 전체 log-rank p와 Holm p<0.05인 그룹 쌍을 표시합니다. Variant count는 유지된 MAF 행 수이며 임상 TMB(mutations/Mb)와 별개입니다.

전체 R 분석이 성공한 뒤 원고 구성용 Word 보고서를 별도로 생성합니다.

```sh
python3 -m pip install python-docx pandas Pillow
python3 scripts/build_results_report.py --results outputs/updated_v3_20260928
```

Figure 25개, 편집 가능한 Table 7개(일부 패널 분할), 분석·해석·주의사항, Supplementary table 파일명을 포함합니다. 결과 폴더에 DOCX와 Table1–7 TSV, 보조표 번호·그림 목록이 생성됩니다. 환자 단위 자료는 공개용 보조표에서 제외됩니다. 결과 폴더는 계속 Git 제외 대상입니다. [상세 방법](docs/updated_analysis.md)을 참고하십시오.

## HRD 관련·MMR 유전자군 추가 분석

`scripts/repair_gene_profiles.R`는 McIntyre 2020의 HRD 관련 18개 및 MMR 4개 유전자로 전체/변이 보유군 oncoplot, 임상 분포와 쌍별 통계표를 만듭니다. `run_all.R`의 마지막 단계이며 결과 접두사는 31–32입니다. 유전자 목록은 `R/repair_gene_sets.R`에서 확인할 수 있습니다.

목록 내 변이는 기능적 HRD/dMMR 진단이 아닙니다. LOH, CNV, tumor purity, tumor cellularity는 분석·비교·그림 주석에서 제외하며 원본 값은 보존합니다. 보고서는 Pretendard 본문·표를 사용하며, 읽는 컴퓨터에도 해당 폰트를 설치하는 것이 좋습니다. 그림용 명칭은 `R/data.R`의 `figure_label()`에서 한 번에 확인할 수 있습니다.
