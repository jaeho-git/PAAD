# 입력자료 사전

**새 `curated_v3` 임상자료:** 이 문서의 아래 내용은 기존 v19 입력의 정의입니다. 새 workbook의 `Data_dictionary` 기반 매핑과 생존 endpoint, 데이터 우선순위는 [updated_analysis.md](updated_analysis.md)를 우선합니다. 새 자료에서 대상 TXT는 감사 비교용이며 임상 값을 제공하지 않습니다.

이 문서는 공개 가능한 입력 schema만 설명한다. 실제 환자·검체 값, 식별자, 날짜, 개인 경로는 포함하지 않는다. 실제 입력은 로컬의 `data/raw/`에 복제하여 사용하며 Git에서 제외한다. `data/example/`의 파일은 같은 형식을 설명하기 위해 처음부터 만든 합성 자료이며 Git에 포함할 수 있다.

## 파일 수준 정의

| 설정 항목 | 합성 예제 | 형식 | 분석에서의 역할 |
| --- | --- | --- | --- |
| `clinical_file` | `data/example/PDAC_clinical_info.example.xlsx` | Excel workbook, sheet `명단` | 기본 임상·병리 annotation |
| `target_patients_file` | `data/example/result_527명.example.txt` | tab-delimited text, CP949, CRLF | 분석 대상 barcode, 생존 및 병리 변수를 제공 |
| `maf_file` | `data/example/PDAC_oncopanel.example.maf` | tab-delimited MAF, UTF-8 | 체세포 변이와 sample barcode 제공 |

세 파일의 연결 key는 `Tumor_Sample_Barcode`다. v19 호환 순서상 최초 cohort 포함 검사는 두 입력을 읽은 그대로 비교하고, 대상 목록의 barcode는 join 직전에, 임상 barcode는 join 뒤에 앞뒤 공백을 제거한다. 따라서 입력의 공백 차이는 포함 또는 join 결과에 영향을 줄 수 있으므로 실제 자료에서 별도 점검해야 한다. 대상 목록에 같은 barcode가 여러 번 있으면 현재 동작은 첫 행을 남긴다. 임상표도 마지막에 barcode별 첫 행을 남기므로, 중복 행 선택 규칙의 과학적 적합성은 실제 자료에서 확인해야 한다.

## 임상 XLSX

원본에서 확인한 worksheet 이름은 `명단`이며 32개 header를 가진다. 새 코드는 그중 아래 컬럼을 직접 사용한다. 합성 workbook은 원본 header 이름과 순서를 유지하지만 실제 값을 포함하지 않는다.

| 컬럼 | 기대 형식 | 필수 | 처리·허용값 | 단위·주의사항 |
| --- | --- | --- | --- | --- |
| `Tumor_Sample_Barcode` | 문자 | 예 | 세 파일의 join key; 공백 제거 | 한 행을 환자 단위로 볼지 검체 단위로 볼지는 실제 자료에서 확인 필요 |
| `Classification` | 문자 | 예 | `Exception` 제외 후 대소문자 구분 없이 `PDAC`를 포함하는 행만 유지 | 정확한 질환 분류 정의 미확인 |
| `Differentiation` | 문자 | 예 | 문자열에 `wd`, `md`, `pd`가 있으면 각각 `WD`, `MD`, `PD`로 정규화 | 병리 등급의 원 임상 정의 미확인 |
| `NAC` | 문자 | 예 | `y`→`Yes`, `n`→`No`; 그 밖의 값은 유지 | 치료 정의·시점 미확인 |
| `Size` | 숫자로 변환 가능 | 예 | `as.numeric()` 변환 | 단위 미확인 |
| `T` | 문자/숫자 | 예 | `1`, `2`, `3`을 `T1`, `T2`, `T3`로 변환; 나머지는 문자열 유지 | staging edition은 코드만으로 확정하지 않음 |
| `N` | 문자/숫자 | 예 | `0`, `1`, `2`를 `N0`, `N1`, `N2`로 변환 | staging edition은 코드만으로 확정하지 않음 |
| `N status` | 문자/숫자 | 아니요 | 로드 후 `N_status`; `0`→`N0`, `1`→`N1`; 없으면 전부 결측 | 원 변수 정의 미확인 |
| `Stage` | 문자 | 예 | 원 값을 유지; `IA/IB`→group I, `IIA/IIB`→II, `III`→III | 그 밖의 stage는 `Stage_Group`에서 결측 처리 |
| `BMI` | 숫자로 변환 가능 | 예 | `as.numeric()` 변환 | 단위 미확인 |
| `CA 19-9` | 문자/숫자 | 예 | 로드 후 `CA19_9`; `<`, `>`, `=`, 공백 제거 뒤 숫자 변환 | 측정 단위·검출한계 처리의 적합성 미확인 |
| `CEA` | 문자/숫자 | 예 | `<`, `>`, `=`, 공백 제거 뒤 숫자 변환 | 측정 단위·검출한계 처리의 적합성 미확인 |
| `성별코드` | 문자 | 정책 의존 | v19 호환 기본값(`legacy`)에서는 `Sex`로 자동 대체하지 않음; `clinical_korean` 정책에서만 명시적으로 매핑 | 코드 체계 미확인 |

XLSX의 나머지 header는 형식 재현을 위해 합성 예제에 존재할 수 있지만 현재 분석에는 전달하지 않는다.

## 대상 목록 TXT

파일은 CP949로 읽는다. 아래는 v19에서 실제로 join 또는 분석에 사용하는 컬럼이다.

| 컬럼 | 기대 형식 | 필수 | 처리·코딩 |
| --- | --- | --- | --- |
| `Tumor_Sample_Barcode` | 문자 | 예 | join key, 공백 제거 |
| `ID_NGS` | 문자 | 아니요 | join 단계에는 포함되나 공개 출력용 최종 임상 컬럼에서는 제거 |
| `age` | 문자/숫자 | 아니요 | 숫자와 점 이외 문자를 제거하여 `Age`로 변환; ≤55, 56–69, ≥70 세 그룹 생성 |
| `OS_m`, `RFS_m` | 숫자로 변환 가능 | 생존분석 조건부 | 코드가 month-scale 변수로 해석; 실제 임상 정의 미확인 |
| `OS_d`, `RFS_d` | 숫자로 변환 가능 | 생존분석 조건부 | month 변수가 전부 결측이면 day-scale fallback; 실제 임상 정의 미확인 |
| `survive` | 0/1 숫자 | OS 조건부 | 현재 코드는 `1=event`, `0=censor`로 사용; 컬럼명과 임상 event 정의는 사용자 확인 필요 |
| `Recur` | 0/1 숫자 | RFS 조건부 | 현재 코드는 `1=event`, `0=censor`로 사용; 임상 정의 확인 필요 |
| `LVI`, `PNI`, `RM` | 문자/숫자 | 아니요 | `1/positive/pos/yes/y`→`Positive`; `0/negative/neg/no/n`→`Negative`; 그 밖의 값 유지 |
| `sex` | 문자 | 정책 의존 | `sex_column_policy="target_patients"`일 때만 `Sex`로 매핑; coding 미확인 |

동명 컬럼이 임상 XLSX에도 있으면 `ID_NGS`, 생존 변수, `LVI/PNI/RM`, age 관련 컬럼은 임상표에서 제거한 뒤 대상 목록 값을 우선 join한다. 이는 v19의 기존 우선순위를 보존한 것이다.

결측값은 R reader가 인식하는 빈 셀/빈 field와 `NA`로 처리된다. 원자료에서 쓰인 모든 결측 표기 규칙은 확인되지 않았으므로 추가 문자열 sentinel이 있다면 로컬 검증이 필요하다.

## MAF

MAF는 `maftools::read.maf()`로 읽으며 임상 cohort의 barcode로 `subsetMaf()`를 수행한다. 분석에 직접 필요한 핵심 컬럼은 다음과 같다.

| 컬럼 | 기대 형식 | 역할 |
| --- | --- | --- |
| `Hugo_Symbol` | 문자 | 유전자별 빈도, 상위 유전자, 4개 driver 상태 |
| `Chromosome` | 문자 | MAF 기본 locus |
| `Start_Position`, `End_Position` | 정수형 좌표 | MAF 기본 locus; protein annotation이 없을 때 KRAS subtype fallback에도 사용 |
| `Reference_Allele`, `Tumor_Seq_Allele2` | 문자 | allele 및 KRAS 좌표 fallback |
| `Variant_Classification` | 문자 | oncoplot mutation type; maftools의 retained variant 결정 |
| `Variant_Type` | 문자 | MAF 기본 변이 유형 |
| `Tumor_Sample_Barcode` | 문자 | 임상표와 연결되는 sample key |
| `NCBI_Build` | 문자 | 좌표계 metadata; KRAS fallback 좌표와 일치하는 build인지 확인 필요 |
| `HGVSp_Short`, `Protein_Change`, `Amino_Acid_Change`, `AAChange`, `HGVSp` | 문자 | 하나 이상 있으면 KRAS amino-acid subtype 판정에 우선 사용; 모두 없으면 좌표/allele fallback |

합성 MAF에는 reader 호환을 위한 표준 확장 컬럼이 더 포함된다. 현재 코드의 “TMB”는 `maftools` 객체에 남은 MAF 행을 sample별로 센 값이며 panel 크기로 나눈 mutations/Mb가 아니다.

## 파생 변수

| 파생 변수 | 정의 |
| --- | --- |
| `Age` | 대상 목록의 `age`에서 숫자/점을 제외한 문자를 제거하고 numeric 변환 |
| `Age_group` | `Age <=55`, `Age 56-69`, `Age >=70` |
| `Stage_Group` | IA/IB→I, IIA/IIB→II, III→III, 그 밖은 결측 |
| `KRAS_subtype` | 단일 subtype이 G12D/G12V/G12R이면 해당 값; 그 밖의 KRAS 변이는 `Other KRAS`; KRAS 변이가 없으면 `WT` |
| `driver_mutation_count` | KRAS, TP53, SMAD4, CDKN2A 중 변이가 있는 유전자 수(0–4) |
| `driver_mutation_count_collapsed` | 0, 1–2, 3–4 mutations |
| 원본 명칭 `TMB` | sample별 retained MAF row 수; 정규화된 종양변이부담 아님 |

## 합성 예제의 보장 범위

`data/example/generate_examples.py`는 외부 자료를 읽지 않고 결정론적으로 세 예제를 만든다. 다음 명령은 schema, encoding, line ending 및 기대 hash를 확인한다.

```sh
python3 data/example/generate_examples.py --check
```

이 검사는 파일 형식과 합성 fixture의 안정성만 보장한다. 실제 cohort 선정, 임상 정의, 통계 결과 또는 원본 결과 동등성을 검증하지 않는다.
