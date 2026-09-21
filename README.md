# PAAD 분석 코드

PAAD 환자의 targeted exome sequencing 자료와 임상정보를 이용한 PDAC 분석 코드를 정리한 저장소입니다. 기준 통합 스크립트 `PDAC_plot_n527_modified_v19_publication_update.R`의 분석 의미와 출력 형식을 유지하면서, 공통 처리와 분석별 실행 파일을 분리했습니다.

이 저장소는 공개 가능한 코드·문서·합성 예제만 Git으로 관리합니다. 실제 환자자료와 새 분석 결과는 로컬에만 보관합니다.

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
├── R/                         # 입력 처리, plot, driver/survival, workflow 함수
├── scripts/                   # 분석별 실행 진입점 7개
├── config/
│   ├── config.example.R       # 실제 로컬 입력용 공개 설정 예시
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
├── run_all.R                  # 전체 분석 진입점
└── PAAD.Rproj
```

## 실행환경

검증한 환경은 R 4.5.3입니다. 분석에 필요한 package는 다음과 같습니다.

- CRAN: `readxl`, `writexl`, `tidyverse`, `circlize`, `RColorBrewer`, `survival`
- Bioconductor: `maftools`, `ComplexHeatmap`

스크립트는 package를 자동 설치하거나 사용자의 global library를 변경하지 않습니다. 누락된 package가 있으면 실행을 중단하고 이름을 표시합니다. 새 환경에서는 R console에서 다음과 같이 준비할 수 있습니다.

```r
install.packages(c(
  "readxl", "writexl", "tidyverse",
  "circlize", "RColorBrewer", "survival"
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
- 합성 입력을 이용한 전체 workflow 실행
- `outputs/synthetic_run/`에 예상 산출물 42개가 생성되는지 확인

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

`분석용데이터1_수작업_eOJH_sdk0116.xlsx`도 현재 로컬 `data/raw/`에 복제되어 있지만, 기준 v19 workflow에서는 사용하지 않습니다.

공개 설정을 로컬 설정으로 복사한 뒤 필요한 경로만 수정합니다.

```sh
cp config/config.example.R config/local.R
```

현재 작업 환경에는 위 실제 입력 사본과 `config/local.R`이 이미 준비되어 있습니다. 두 항목 모두 Git에서 제외됩니다.

## 전체 분석 실행

```sh
Rscript --vanilla run_all.R --config=config/local.R
```

모든 결과는 기본적으로 `outputs/local_run/`에 새로 생성됩니다. 설정 검증은 `output_dir`이 이 저장소의 `outputs/` 또는 그 하위인지 확인하며, 원본 프로젝트나 임의의 외부 경로에는 결과를 쓰지 않습니다.

## 분석별 실행

각 실행 파일은 새로운 R session에서 독립적으로 입력과 설정을 준비합니다.

| 분석 | 실행 파일 | 주요 출력 |
| --- | --- | --- |
| Main 및 T/N-filtered oncoplot | `scripts/01_oncoplots.R` | `1_*.png` |
| 임상 요약표 | `scripts/02_clinical_summary.R` | `2_Clinical_Table_PDAC.xlsx` |
| 임상 그룹별 변이빈도·변이개수 비교 | `scripts/03_comparison_plots.R` | `4_FreqPlot_*.png`, `4_TMB_BoxPlot_*.png` |
| Top 20/5 mutation 및 변이개수 heatmap | `scripts/04_heatmaps.R` | `5_*.png`, `6_*.png`, `7_*.png` |
| Somatic interaction | `scripts/05_somatic_interactions.R` | `9_*.png`, `9_*.xlsx` |
| Stage count | `scripts/06_stage_counts.R` | `10_Stage_Counts_PDAC.png` |
| KRAS·driver·clinicopathology·survival | `scripts/07_driver_kras_survival.R` | `11_*`–`17_*` |

예를 들어 실제 로컬 데이터로 임상 요약표만 생성하려면 다음과 같이 실행합니다.

```sh
Rscript --vanilla scripts/02_clinical_summary.R --config=config/local.R
```

원본 구간, 정확한 출력 파일과 현재 검증 상태는 [docs/analysis_map.md](docs/analysis_map.md)에 정리되어 있습니다.

## 주요 설정

`config/local.R`과 `config/config.synthetic.R`에서 다음 항목을 관리합니다.

| 항목 | 의미 |
| --- | --- |
| `clinical_file`, `target_patients_file`, `maf_file` | 세 입력 파일 |
| `output_dir` | 반드시 저장소의 `outputs/` 하위 |
| `target_encoding` | cohort TXT 인코딩; 현재 `CP949` |
| `sex_column_policy` | `legacy`, `target_patients`, `clinical_korean` 중 선택 |
| `random_seed` | `NULL`이면 v19의 비결정적 simulated Fisher 동작 유지; 정수이면 재현 가능 |
| `top_n` | 기본 상위 유전자 수 20 |
| `plot_dpi` | production 기본값 1000; smoke test는 150 |

기본 실제 데이터 설정은 `sex_column_policy="legacy"`입니다. 원본 clinical workbook에는 `Sex`가 아니라 `성별코드`, cohort TXT에는 `sex`가 있으므로 v19과 동일하게 성별 annotation을 자동 대체하지 않습니다. 어느 source와 coding을 사용할지 확인한 뒤에만 정책을 변경해야 합니다.

## 검증 결과

- 기준 v19 스크립트 전체 1,279행의 입력·분석·출력을 조사했습니다.
- 합성 입력의 전체 workflow가 새 R session에서 완료되어 42개 파일을 생성했습니다.
- 실제 로컬 입력은 v19 filter 후 임상 509행, retained MAF 4,224행으로 로드됐습니다.
- 실제 자료로 각 workflow를 격리된 `outputs/code_validation/`에서 실행했습니다. 무거운 figure workflow는 top 5 검증 설정을 사용했습니다.
- 기존 결과와 임상 요약 workbook 12개 sheet, 509×35 driver/KRAS 표, 509×35 survival input 표, risk TSV 3개가 일치했습니다.
- Driver-clinicopathology 표는 구조와 p-value 이외 필드가 일치했습니다.
- simulated Fisher p-value의 기존 RNG 상태와 PNG의 pixel/figure 동등성은 확인하지 않았습니다.
- production 설정의 top 20·1000 DPI 전체 실제 데이터 `run_all.R`은 실행하지 않았습니다.

상세 변경과 한계는 [docs/refactoring_notes.md](docs/refactoring_notes.md)를 참조합니다.

## 해석 시 주의사항

- 기존 출력에서 `TMB`로 표시된 값은 panel 크기로 정규화한 mutations/Mb가 아니라 sample별 retained MAF row 수입니다. 이번 구조 정리에서는 계산과 기존 파일명을 임의로 바꾸지 않았습니다.
- 변이빈도 및 driver 연관성 분석의 일부 Fisher test는 simulation을 사용하고 multiple-testing correction이 적용되지 않았습니다.
- OS는 `survive`, RFS는 `Recur`를 event로 사용하며 `1=event`, `0=censored`로 해석합니다. 실제 임상 사전 확인이 필요합니다.
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
