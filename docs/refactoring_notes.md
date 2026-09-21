# 리팩터링 및 검증 기록

## 범위

기준은 `PDAC_plot_n527_modified_v19_publication_update.R`이다. 통합 스크립트의 cohort 정의와 분석 로직을 공통 모듈과 일곱 개 실행 파일로 분리했으며, 원본 프로젝트는 조사 중 읽기 전용으로 유지했다.

사용자 요청에 따라 다음 저장 정책을 적용했다.

- 실제 `00_Data` 입력은 대상 저장소의 `data/raw/`에 로컬 실행용으로 복제하되 Git에서 제외한다.
- 실제 환자 값, 식별자, 날짜, 개인 절대경로는 공개 코드·문서·설정 예시에 포함하지 않는다.
- `data/example/`에는 원본 header와 파일 형식을 재현한 완전 합성 fixture만 둔다.
- 기존 `01_Results`와 `01_Results_copy`는 `outputs/`로 복사하거나 백업하지 않는다.
- `outputs/`는 새 코드 실행으로 생성되는 결과만 저장하며 Git에서 제외한다.
- Git staging, commit, push와 `.git` metadata 변경은 수행하지 않는다.

## 구조 정리에서 적용한 변경

아래 항목은 분석 정의를 바꾸기 위한 수정이 아니라 독립 실행과 공개 저장소 운영을 위한 구조 변경이다.

| 변경 | 구현 | 결과 의미에 대한 처리 |
| --- | --- | --- |
| 개인 절대경로 제거 | `config/config.example.R` 및 Git 제외 대상 `config/local.R` | 입력 선택은 설정으로 이동; 로직은 유지 |
| 안전한 출력 경로 | `output_dir`을 저장소 `outputs/` 또는 하위 경로로 제한 | 기존 결과를 덮어쓰지 않음 |
| 공통 module 분리 | `R/config.R`, `R/data.R`, plot/driver/workflow module | 개별/전체 실행이 같은 함수 사용 |
| 분석별 entrypoint | `scripts/01_*.R`부터 `07_*.R`, `run_all.R` | v19 실행부의 순서와 출력명을 보존 |
| 이전 session 의존 제거 | 각 entrypoint가 bootstrap, 설정, 입력을 명시적으로 로드 | 새 R session에서 실행 가능 |
| package 명시적 사전 점검 | 필요한 package를 점검하고 누락 시 종료 | 실행 중 사용자의 global library를 변경하지 않음 |
| 실제/합성 자료 분리 | `data/raw/`는 Git 제외, `data/example/`은 합성 fixture | 합성 실행을 실제 연구 검증으로 간주하지 않음 |
| Sex source ambiguity 표면화 | `sex_column_policy`에 `legacy`, `target_patients`, `clinical_korean` 제공 | 기본값은 v19처럼 자동 대체하지 않음; 대체는 사용자 명시 필요 |
| 모의 검정 seed 선택지 | `random_seed=NULL` 기본, 정수 설정 가능 | 기본은 v19 비결정성을 보존; reproducibility 선택 가능 |
| 용어 경고 | 실행 시 원본 `TMB`가 변이 행 수임을 안내 | 계산식을 조용히 변경하지 않음 |
| 명시적 출력 인자 | comparison·driver 함수에 `output_dir`/파일 경로를 전달 | 전역 `output_dir` 의존성을 제거하되 기본 출력명은 유지 |
| Somatic interaction device 정리 | plot device 안에서 `somaticInteractions()`를 한 번 실행해 표와 그림에 같은 반환값 사용 | v19의 중복 실행과 stray `Rplots.pdf` 가능성을 제거; figure 동등성은 별도 미확인 |
| 상수 변이개수 heatmap 보호 | 모든 sample의 변이 행 수가 같을 때 color scale 범위를 소폭 확장 | 값 자체는 바꾸지 않고 v19이 실패하던 퇴화 입력에서만 렌더링 가능하게 함 |

## 보존한 v19 분석 정의

- 대상 목록으로 임상 행을 제한한 뒤 `Classification != "Exception"` 및 `grepl("PDAC", ...)`을 적용한다.
- 대상 목록의 생존 변수와 `LVI/PNI/RM`이 임상 workbook의 동명 컬럼보다 우선한다.
- 임상 범주·수치 변환, annotation 순서, top gene 선택, driver gene 네 개, survival grouping을 유지한다.
- figure 크기, 해상도, 색상과 기존 output filename을 유지한다.
- v19에 없는 Output 3·8을 임의로 만들지 않는다.

## 사용자 판단이 필요한 분석·재현성 이슈

다음 사항은 결과를 바꿀 수 있으므로 이번 구조 정리에서 임의 수정하지 않았다.

1. **`TMB` 정의:** sample별 retained MAF row 수이며 panel territory로 나누지 않는다. 정규화된 mutations/Mb로 해석하거나 표기하려면 panel size와 variant counting rule을 정한 뒤 별도 변경해야 한다.
2. **다중검정:** subgroup별 상위 유전자 Fisher 검정과 driver-clinicopathology 검정에 multiple-testing correction이 없다. 보정 방식과 family 정의가 필요하다.
3. **모의 Fisher 재현성:** v19는 simulated Fisher test에 seed를 지정하지 않는다. 기본 설정은 이를 보존한다. 고정 결과가 필요하면 `random_seed`를 정하되 원본과의 직접 비교 조건을 함께 기록해야 한다.
4. **KRAS 좌표 fallback:** amino-acid annotation이 없을 때 특정 genomic coordinate와 allele을 사용한다. MAF의 genome build, strand 표현, panel annotation과의 일치를 확인해야 한다.
5. **Sex 컬럼:** 조사한 clinical workbook에는 `Sex` 대신 `성별코드`, 대상 목록에는 `sex`가 있으나 v19는 둘을 `Sex`로 매핑하지 않는다. 어느 source와 coding을 쓸지 연구자가 선택해야 한다.
6. **생존 event:** 코드는 OS에 `survive`, RFS에 `Recur`를 사용하고 1을 event로 처리한다. 실제 임상 사전과 censoring 정의를 확인해야 한다.
7. **환자/검체 단위:** barcode별 첫 행을 남기고 이를 figure에서 patient N으로 표현한다. 한 환자에 여러 검체가 존재할 수 있다면 환자-검체 mapping과 대표 검체 규칙을 정해야 한다.
8. **중복 barcode:** 대상 목록과 임상표에서 첫 행을 보존한다. 중복이 있을 때 어떤 행이 대표인지 명시적 기준이 없다.
9. **범주 변환:** `Stage_Group`은 I–III 일부만 만들며 T 변환은 1–3만 `T` prefix로 바꾼다. 실제 자료의 T0/T4, Stage 0/IV와 staging edition을 확인해야 한다.
10. **검출한계 문자열:** CA 19-9와 CEA에서 부등호를 제거해 경계 숫자로 사용한다. censored measurement 처리로 적합한지 검토가 필요하다.
11. **모델 범위:** survival output은 unadjusted Kaplan–Meier/log-rank 비교다. confounder-adjusted model로 해석하지 않아야 한다.

## 검증 기록

| 점검 | 상태 | 범위와 한계 |
| --- | --- | --- |
| 기준 코드 분석 목록과 입출력 추적 | 완료 | v19 전체 1,279행의 section, input reader, output call을 정적으로 확인 |
| 공개 경로의 개인 절대경로 검사 | 완료 | 공개 후보 code/config/docs/tests를 대상으로 검사; 실제 값은 기록하지 않음 |
| 합성 fixture 자체 검사 | 완료 | Python 3.9에서 `generate_examples.py --check`로 schema·encoding·hash 확인 |
| R module parse/source | 완료 | 새 R session에서 모든 module과 entrypoint 문법·참조 점검 |
| 합성 입력 load/preprocess | 완료 | CP949 TXT, XLSX, MAF join 및 PDAC filter 확인 |
| 합성 예제 개별 실행 | 완료 | 임상 요약, Stage count, 핵심 driver 집계와 CLI 진입점을 확인 |
| 합성 예제 전체 실행 | 완료 | 24개 합성 입력행 중 v19 filter 후 20개 PDAC 행과 304개 MAF 행 중 300개 retained 행으로 전체 workflow를 실행해 42개 파일 생성; 과학적 결과 검증이 아님 |
| 실제 `data/raw/` load | 완료 | v19 filter 후 임상 509행과 retained MAF 4,224행 확인 |
| 실제 자료 workflow 실행 | 부분 완료 | 각 workflow를 격리된 `outputs/code_validation/`에 실행; oncoplot/comparison/heatmap은 top 5 검증 설정, production top 20·1000 dpi `run_all.R`은 미실행 |
| 원본 대비 표 수치 동등성 | 부분 확인 | 임상 요약 12개 sheet, 509×35 driver/KRAS·survival input 표, risk TSV 3개 값 동일; driver clinicopathology 표는 p-value 이외 필드 동일 |
| 원본 대비 p-value·figure 동등성 | 미확인 | simulated Fisher의 기존 RNG 조건을 재현하지 못했고 PNG pixel/figure 비교는 수행하지 않음 |

구체적인 실행 로그와 생성 결과는 Git 제외 대상 `outputs/`에만 둔다. 위 “완료” 상태는 최종 검증 명령이 성공한 뒤 유지하며, 실패가 확인되면 상태와 원인을 갱신한다.

검증 workstation의 R은 4.5.3이었다. 별도의 `renv.lock`으로 package version portability를 검증한 것은 아니다. 합성 전체 실행은 exit code 0이었으며, 동률 때문에 Wilcoxon exact p-value를 계산할 수 없다는 비치명적 warning이 두 차례 발생했다. 이 경우 R이 사용한 근사 p-value를 실제 연구 결과로 해석하기 전 확인해야 한다.

## 공개 전 확인사항

- `data/raw/`, `config/local.R`, `outputs/`, `local_notes/`가 Git에서 제외되는지 다시 확인한다.
- 공개 후보 text file에 개인 경로·실제 식별자가 없는지 확인한다.
- 합성 fixture가 실제 record를 변형한 것이 아니라 generator의 고정 가상 값만으로 만들어졌는지 확인한다.
- 사용자 판단 이슈를 해결해 분석 정의를 바꾸는 경우 v19 호환 결과와 새 결과를 별도로 표시한다.
