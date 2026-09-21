# 분석 스크립트 안내

이 디렉터리에는 연구자가 직접 읽고 실행할 수 있는 분석 파일 일곱 개가 있습니다. 파일명은 실행 순번보다 분석 내용을 우선해 정했습니다. 각 파일 안에 해당 분석의 함수와 실행 순서가 함께 있으므로, 분석 정의를 확인하기 위해 별도의 실행 중계 파일을 찾아갈 필요가 없습니다.

스크립트 이름의 순번은 제거했지만 생성 결과의 앞 번호는 기존 figure·표와의 대응을 보존하기 위해 유지했습니다.

공통 파일인 `R/config.R`, `R/data.R`, `R/plot_helpers.R`에는 설정 검증, 입력 전처리, 안전한 출력 경로와 반복되는 그림 보조 기능만 둡니다. 통계 검정, 표 구성, 유전자 선택, survival group 정의와 실제 출력 호출은 아래 분석 파일에서 확인할 수 있습니다.

## 빠르게 찾기

| 알고 싶은 분석 | 먼저 열 파일 | 생성 결과 |
| --- | --- | --- |
| 전체 및 T/N 정보 보유군 oncoplot | `oncoplots.R` | Output 1, 1-2 |
| 임상 특성 요약표 | `clinical_summary.R` | Output 2 |
| 임상 그룹별 변이빈도·변이개수 비교 | `clinical_group_comparisons.R` | Output 4 |
| Top mutation 및 변이개수 heatmap | `mutation_heatmaps.R` | Output 5–7 |
| 유전자 간 somatic interaction | `somatic_interactions.R` | Output 9 |
| Stage별 환자 수 | `stage_distribution.R` | Output 10 |
| KRAS·driver·clinicopathology·생존 | `driver_kras_survival.R` | Output 11–17 |

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

각 분석 파일은 새 R session에서 독립적으로 실행할 수 있으며, 다른 분석 파일이 만든 전역 객체나 중간 결과 파일을 요구하지 않습니다. 결과는 설정의 `output_dir`에 생성됩니다. `config/local.R`의 기본값은 `outputs/local_run/`, 합성 설정의 기본값은 `outputs/synthetic_run/`입니다. 두 경로 모두 Git에서 제외됩니다.

## 코드를 읽는 순서

개별 파일은 위에서 아래로 다음 흐름을 따릅니다.

1. 저장소 위치와 `--config` 인자를 확인합니다.
2. 공통 설정·입력·그림 보조 파일을 불러옵니다.
3. 그 분석에만 필요한 함수와 통계 규칙을 정의합니다.
4. 설정과 세 입력 파일을 읽습니다.
5. 분석을 실행하고 이름이 고정된 결과를 `output_dir`에 저장합니다.

정확한 원본 코드 구간, 입력 조건, 출력 파일명과 검증 범위는 [`docs/analysis_map.md`](../docs/analysis_map.md)를 참조하십시오. 합성 예제는 소프트웨어 실행 확인용이며 과학적 결과 또는 통계적 타당성을 검증하는 자료가 아닙니다.
