# 입력 데이터 준비

이 저장소는 실제 환자별 임상자료와 변이자료를 Git에 저장하지 않습니다. 로컬 실제 입력과 공개 가능한 합성 예제를 분리합니다.

## 디렉터리 구분

| 경로 | 용도 | Git 관리 |
| --- | --- | --- |
| `data/raw/` | 로컬에서 분석에 사용하는 실제 입력 파일의 복제 위치 | 제외 |
| `data/processed/` | 로컬 전처리 자료 | 제외 |
| `data/private/` | 그 밖의 비공개 자료 | 제외 |
| `data/example/` | 입력 형식 확인과 smoke test를 위한 완전 합성 예제 | 포함 |

`data/raw/`의 자료는 로컬 실행을 위한 사본입니다. 원본과 마찬가지로 환자·검체 정보가 포함될 수 있으므로 공유하거나 Git에 추가하지 않습니다. 분석 결과는 `outputs/`에 새로 생성하며, 기존 결과 폴더를 이 위치로 복사하지 않습니다.

## 합성 예제

`data/example/`의 식별자는 모두 `SYN-`으로 시작합니다. 행과 값은 원자료를 복사·익명화·변형한 것이 아니라 generator에 명시한 규칙으로 처음부터 생성했습니다. 이 자료는 파일 파싱과 실행 흐름을 점검하기 위한 것이며, 연구 결과 재현이나 과학적·통계적 검증에는 사용할 수 없습니다.

| 파일 | 형식 | 내용 |
| --- | --- | --- |
| `PDAC_clinical_info.example.xlsx` | XLSX, 시트 `명단` | 실제 임상 입력과 같은 32개 헤더, 첫 행 고정, 자동 필터, 24개 합성 행 |
| `result_527명.example.txt` | CP949, 탭 구분, CRLF | 실제 코호트 입력과 같은 58개 헤더, 24개 합성 행 |
| `PDAC_oncopanel.example.maf` | UTF-8, 탭 구분, CRLF | 실제 MAF 입력과 같은 37개 헤더, 합성 변이 행 |

MAF의 첫 두 줄은 `#`으로 시작하는 공개용 주석입니다. `synthetic_example=true`는 실제 연구자료가 아님을 표시합니다. 예제의 유전자명과 MAF 범주명은 파서 동작을 확인하기 위한 공개 도메인 용어입니다. KRAS G12D/G12V/G12R 및 Q61 fallback smoke test에는 분석 코드가 참조하는 공개 GRCh37 좌표를 사용합니다. 이 좌표는 환자자료가 아닙니다. 그 밖의 `SYN_GENE_*` 좌표, 검체 식별자와 측정값은 합성 값입니다.

## 스키마

### Clinical XLSX

- 파일 단위: 검체와 연결된 임상 record
- 조인 키: `Tumor_Sample_Barcode`
- 시트: `명단`
- 헤더: `임상리포트ID_NGS`, `Tumor_Sample_Barcode`, `성별코드`, `검사나이`, `접수일자`, `검사코드`, `검사코드2NM`, `검사결과본문내용`, `암종`, `Classification`, `특이사항`, `Differentiation`, `위치`, `NAC`, `Size`, `LV`, `PN`, `RM`, `positive LN`, `Total LN`, `T`, `N`, `N status`, `M`, `Stage`, `키`, `몸무게`, `BMI`, `CA 19-9`, `CEA`, `수술_수술일자`, `수술_수술나이(년)`
- 확인된 형식: `Size`, 림프절 수, TNM 숫자, 신체계측치와 수술 관련 숫자는 numeric으로 표현될 수 있습니다. `CA 19-9`와 `CEA`는 부등호가 포함될 수 있어 원본 입력에서 text입니다.
- 미확인: 각 임상 변수의 정확한 수집 정의, 단위, 허용 범주, 결측 표현은 제공된 코드와 파일 형식만으로 확정하지 않았습니다.

### Cohort TXT

- 파일 단위: 환자·검체 연결 record
- 조인 키: `Tumor_Sample_Barcode`; 보조 식별 필드 `ID_NGS`
- 인코딩과 줄바꿈: CP949, tab-separated, CRLF
- 주요 생존 필드: `OS_m`, `OS_d`, `RFS_m`, `RFS_d`, `survive`, `Recur`
- 코드에서 확인된 event 처리: 분석 코드는 `Surv(time, survive)`와 `Surv(time, Recur)`를 사용하고 두 event 변수를 `0` 또는 `1`로 제한합니다. 따라서 현재 pipeline은 `1 = event`, `0 = censored`로 해석합니다.
- 미확인: 원자료 수집 단계에서 `survive`라는 변수명에 부여한 임상적 정의와 censoring 기준은 별도 codebook으로 확인되지 않았습니다. 공개 전 실제 codebook과 대조해야 합니다.
- 기타 58개 필드의 전체 헤더는 generator의 `COHORT_HEADERS`에 고정되어 있습니다. 변수별 임상 단위와 범주 정의는 확인되지 않은 항목을 추정하지 않습니다.

### MAF

- 파일 단위: 한 행당 변이 record
- 검체 키: `Tumor_Sample_Barcode`
- 인코딩: UTF-8, tab-separated
- 필수 핵심 필드: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Variant_Classification`, `Variant_Type`, `Reference_Allele`, `Tumor_Seq_Allele1`, `Tumor_Seq_Allele2`, `Tumor_Sample_Barcode`
- 추가 입력 필드: matched-normal allele, validation 상태, sequencing metadata, UUID, `t_depth`, `t_ref_count`, `t_alt_count`를 포함해 실제 입력의 37개 헤더를 유지합니다.
- KRAS fallback 좌표는 공개 GRCh37 reference locus이며 나머지 좌표와 모든 측정값은 합성입니다. 어떤 값도 생물학적 결과로 해석할 수 없습니다.
- 20개 PDAC 합성 검체는 driver mutation count 0-4에 각각 4개씩 배치하고, KRAS subtype `G12D`, `G12V`, `G12R`, `Other KRAS`, `WT`에 각각 4개씩 배치합니다. 추가 4개 행은 non-PDAC 또는 `Exception` filter 경로를 점검합니다.
- `SYN_GENE_*`를 포함해 25개 이상의 gene을 반복 배치하여 interaction 단계가 빈 입력 때문에 중단되지 않게 했습니다. 이는 통계적 검증 설계가 아닙니다.

## 재생성 및 점검

generator는 Python 표준 라이브러리만 사용하며 `PAAD_tmp`, `data/raw/`, 환경변수 또는 외부 자료를 읽지 않습니다.

```bash
python3 data/example/generate_examples.py
python3 data/example/generate_examples.py --check
```

첫 명령은 세 예제 파일을 결정적으로 다시 만듭니다. 두 번째 명령은 기존 파일이 generator의 예상 byte와 정확히 같은지, XLSX 시트·헤더·고정행·필터, TXT의 CP949·CRLF, MAF의 주석·헤더 계약을 만족하는지 확인합니다. 각 실행은 파일별 SHA-256을 출력합니다.

## 실제 입력 파일

로컬 `data/raw/`에는 다음 이름으로 실제 입력을 준비합니다.

```text
data/raw/PDAC_clinical_info.xlsx
data/raw/result_527명.txt
data/raw/PDAC_oncopanel.maf
data/raw/분석용데이터1_수작업_eOJH_sdk0116.xlsx
```

`분석용데이터1_수작업_eOJH_sdk0116.xlsx`는 조사·추적을 위해 보관하는 선행 수작업 workbook입니다. 현재 v19 기준 실행의 직접 입력은 아닙니다. 공개 저장소에 포함하지 않으며 내부 값은 공개 문서에서 설명하지 않습니다.

실제 경로 선택은 공개 가능한 예제 설정과 Git에서 제외되는 로컬 설정을 구분해 관리합니다. 실제 자료를 합성 예제 파일 위에 덮어쓰지 않습니다.
