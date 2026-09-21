#!/usr/bin/env python3
"""Generate deterministic, fully synthetic PAAD example inputs.

This generator is intentionally self-contained. It never reads PAAD_tmp,
data/raw, environment variables, or any external data source.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape


OUTPUT_DIR = Path(__file__).resolve().parent
CLINICAL_NAME = "PDAC_clinical_info.example.xlsx"
COHORT_NAME = "result_527명.example.txt"
MAF_NAME = "PDAC_oncopanel.example.maf"

CLINICAL_HEADERS = [
    "임상리포트ID_NGS",
    "Tumor_Sample_Barcode",
    "성별코드",
    "검사나이",
    "접수일자",
    "검사코드",
    "검사코드2NM",
    "검사결과본문내용",
    "암종",
    "Classification",
    "특이사항",
    "Differentiation",
    "위치",
    "NAC",
    "Size",
    "LV",
    "PN",
    "RM",
    "positive LN",
    "Total LN",
    "T",
    "N",
    "N status",
    "M",
    "Stage",
    "키",
    "몸무게",
    "BMI",
    "CA 19-9",
    "CEA",
    "수술_수술일자",
    "수술_수술나이(년)",
]

COHORT_HEADERS = [
    "ID_NGS",
    "Tumor_Sample_Barcode",
    "ID",
    "OPdate",
    "OS_m",
    "RFS_m",
    "survive",
    "expire_date",
    "OS_d",
    "Recur",
    "Recur_date",
    "RFS_d",
    "recur_carte2",
    "Recur_treat",
    "sex",
    "age",
    "BMI",
    "NACT2",
    "PreOpCA",
    "CEA",
    "preopSUVmax",
    "preHb",
    "postHb",
    "opname",
    "Open_Lap",
    "OP_name",
    "Op_vessel",
    "Op_organ",
    "intraop_Tf_1",
    "Optime",
    "discharge.date",
    "POD",
    "size",
    "diff",
    "CAP_new",
    "T_8th",
    "N_8th",
    "M_8th",
    "TNM8th",
    "TNM1.4",
    "Positive_LN",
    "Total_LN",
    "LVI",
    "PNI",
    "RM",
    "AdjCTx",
    "AdjReg",
    "AdjReg_name",
    "AdjRTx",
    "nintydayMort",
    "postopCx_1",
    "POPF",
    "POPF_Gr",
    "POPF_new",
    "Clavien",
    "Clavien_new",
    "SevereCx",
    "Clavien_expl",
]

MAF_HEADERS = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "Match_Norm_Seq_Allele1",
    "Match_Norm_Seq_Allele2",
    "Tumor_Validation_Allele1",
    "Tumor_Validation_Allele2",
    "Match_Norm_Validation_Allele1",
    "Match_Norm_Validation_Allele2",
    "Verification_Status",
    "Validation_Status",
    "Mutation_Status",
    "Sequencing_Phase",
    "Sequence_Source",
    "Validation_Method",
    "Score",
    "BAM_File",
    "Sequencer",
    "Tumor_Sample_UUID",
    "Matched_Norm_Sample_UUID",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
]


def synthetic_patients() -> list[dict[str, object]]:
    diffs = ["WD", "MD", "PD"]
    stages = ["IA", "IB", "IIA", "IIB", "III"]
    return [
        {
            "index": i,
            "report_id": f"SYN-CLIN-{i:03d}",
            "sample_id": f"SYN-TUMOR-{i:03d}",
            "normal_id": f"SYN-NORMAL-{i:03d}",
            "sex": "F" if i % 2 == 0 else "M",
            "age": 45 + ((i * 7) % 35),
            "diff": diffs[(i - 1) % len(diffs)],
            "nac": "y" if i % 3 == 0 else "n",
            "stage": stages[(i - 1) % len(stages)],
            "t": 1 + ((i - 1) % 3),
            "n": (i - 1) % 3,
            "classification": (
                "PDAC_SYNTHETIC"
                if i <= 20
                else ("OTHER_SYNTHETIC" if i <= 22 else "Exception")
            ),
        }
        for i in range(1, 25)
    ]


def clinical_rows() -> list[list[object]]:
    rows: list[list[object]] = []
    for p in synthetic_patients():
        i = int(p["index"])
        height = 160 + i * 2
        weight = 50 + i * 3
        bmi = round(weight / ((height / 100) ** 2), 2)
        rows.append(
            [
                p["report_id"],
                p["sample_id"],
                p["sex"],
                str(p["age"]),
                f"2099-01-{i:02d}",
                "SYN-TEST-CODE",
                "SYN-TEST-NAME",
                "SYNTHETIC EXAMPLE ONLY",
                "SYNTHETIC PANCREAS",
                p["classification"],
                "SYNTHETIC",
                p["diff"],
                "SYNTHETIC_SITE",
                p["nac"],
                round(1.5 + i * 0.4, 1),
                "0" if i % 2 else "1",
                "1" if i in (3, 4, 5) else "0",
                "1" if i == 5 else "0",
                int(p["n"]),
                int(p["n"]) + 12,
                int(p["t"]),
                int(p["n"]),
                0 if int(p["n"]) == 0 else 1,
                0,
                p["stage"],
                height,
                weight,
                bmi,
                str(20 + i * 11),
                str(round(1.1 + i * 0.7, 1)),
                20990100 + i,
                int(p["age"]),
            ]
        )
    return rows


def cohort_rows() -> list[list[object]]:
    rows: list[list[object]] = []
    for p in synthetic_patients():
        i = int(p["index"])
        event_os = 1 if (i - 1) % 4 in (1, 3) else 0
        event_rfs = 1 if (i - 1) % 4 in (2, 3) else 0
        row = {
            "ID_NGS": p["report_id"],
            "Tumor_Sample_Barcode": p["sample_id"],
            "ID": f"SYN-PATIENT-{i:03d}",
            "OPdate": f"2099-01-{i:02d}",
            "OS_m": 10 + i * 4,
            "RFS_m": 6 + i * 3,
            "survive": event_os,
            "expire_date": f"2099-12-{i:02d}" if event_os else "",
            "OS_d": (10 + i * 4) * 30,
            "Recur": event_rfs,
            "Recur_date": f"2099-08-{i:02d}" if event_rfs else "",
            "RFS_d": (6 + i * 3) * 30,
            "recur_carte2": "SYN_RECURRENCE" if event_rfs else "SYN_NONE",
            "Recur_treat": "SYN_TREATMENT" if event_rfs else "SYN_NONE",
            "sex": p["sex"],
            "age": p["age"],
            "BMI": round(19.0 + i * 0.8, 1),
            "NACT2": 1 if p["nac"] == "y" else 0,
            "PreOpCA": 20 + i * 11,
            "CEA": round(1.1 + i * 0.7, 1),
            "preopSUVmax": round(2.0 + i * 0.5, 1),
            "preHb": round(14.5 - i * 0.2, 1),
            "postHb": round(12.0 - i * 0.1, 1),
            "opname": "SYN_OPERATION",
            "Open_Lap": "SYN_OPEN" if i % 2 else "SYN_LAP",
            "OP_name": "SYN_PROCEDURE",
            "Op_vessel": "SYN_NONE",
            "Op_organ": "SYN_NONE",
            "intraop_Tf_1": 0,
            "Optime": 200 + i * 10,
            "discharge.date": f"2099-02-{i + 3:02d}",
            "POD": 8 + i,
            "size": round(1.5 + i * 0.4, 1),
            "diff": p["diff"],
            "CAP_new": "SYN_NA",
            "T_8th": p["t"],
            "N_8th": p["n"],
            "M_8th": 0,
            "TNM8th": p["stage"],
            "TNM1.4": p["stage"],
            "Positive_LN": p["n"],
            "Total_LN": int(p["n"]) + 12,
            "LVI": 1 if i % 2 == 0 else 0,
            "PNI": 1 if i in (3, 4, 5) else 0,
            "RM": 1 if i == 5 else 0,
            "AdjCTx": 1 if i != 1 else 0,
            "AdjReg": "SYN_REGIMEN" if i != 1 else "SYN_NONE",
            "AdjReg_name": "SYN_REGIMEN_A" if i != 1 else "SYN_NONE",
            "AdjRTx": 0,
            "nintydayMort": 0,
            "postopCx_1": 1 if i == 4 else 0,
            "POPF": 1 if i == 4 else 0,
            "POPF_Gr": "SYN_B" if i == 4 else "SYN_NONE",
            "POPF_new": 1 if i == 4 else 0,
            "Clavien": 2 if i == 4 else 0,
            "Clavien_new": "SYN_MINOR" if i == 4 else "SYN_NONE",
            "SevereCx": 0,
            "Clavien_expl": "합성예시" if i == 4 else "SYN_NONE",
        }
        rows.append([row[h] for h in COHORT_HEADERS])
    return rows


def maf_rows() -> list[list[object]]:
    driver_genes = ["KRAS", "TP53", "SMAD4", "CDKN2A"]
    generic_genes = [f"SYN_GENE_{i:02d}" for i in range(1, 27)]
    rows: list[list[object]] = []
    patients = synthetic_patients()

    def add_record(
        gene: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        patient_index: int,
        classification: str = "Missense_Mutation",
    ) -> None:
        record_index = len(rows) + 1
        p = patients[patient_index - 1]
        record = {
            "Hugo_Symbol": gene,
            "Entrez_Gene_Id": 0,
            "Center": "SYN_CENTER",
            "NCBI_Build": "GRCh37",
            "Chromosome": chrom,
            "Start_Position": pos,
            "End_Position": pos,
            "Strand": "+",
            "Variant_Classification": classification,
            "Variant_Type": "SNP",
            "Reference_Allele": ref,
            "Tumor_Seq_Allele1": ref,
            "Tumor_Seq_Allele2": alt,
            "dbSNP_RS": "SYN_NONE",
            "dbSNP_Val_Status": "SYN_NOT_TESTED",
            "Tumor_Sample_Barcode": p["sample_id"],
            "Matched_Norm_Sample_Barcode": p["normal_id"],
            "Match_Norm_Seq_Allele1": ref,
            "Match_Norm_Seq_Allele2": ref,
            "Tumor_Validation_Allele1": ref,
            "Tumor_Validation_Allele2": alt,
            "Match_Norm_Validation_Allele1": ref,
            "Match_Norm_Validation_Allele2": ref,
            "Verification_Status": "SYN_NOT_VERIFIED",
            "Validation_Status": "SYN_UNVALIDATED",
            "Mutation_Status": "Somatic",
            "Sequencing_Phase": "SYN_PHASE",
            "Sequence_Source": "SYNTHETIC",
            "Validation_Method": "SYN_NONE",
            "Score": 100 + record_index,
            "BAM_File": "SYN_NONE",
            "Sequencer": "SYN_PLATFORM",
            "Tumor_Sample_UUID": f"SYN-TUMOR-UUID-{patient_index:03d}",
            "Matched_Norm_Sample_UUID": f"SYN-NORMAL-UUID-{patient_index:03d}",
            "t_depth": 80 + (record_index % 20),
            "t_ref_count": 50 + (record_index % 10),
            "t_alt_count": 20 + (record_index % 8),
        }
        rows.append([record[h] for h in MAF_HEADERS])

    # Twenty PDAC examples form five driver-count groups (0-4), four samples each.
    # The KRAS coordinates are public GRCh37 reference loci used by the pipeline's
    # fallback subtype logic; all sample IDs and non-KRAS positions are synthetic.
    kras_by_driver_count = {
        1: ("12", 25398284, "C", "T"),  # G12D
        2: ("12", 25398284, "C", "A"),  # G12V
        3: ("12", 25398285, "C", "G"),  # G12R
        4: ("12", 25380276, "T", "G"),  # Other KRAS (Q61 fallback range)
    }
    for patient_index in range(1, 21):
        driver_count = (patient_index - 1) // 4
        for driver_offset, gene in enumerate(driver_genes[:driver_count]):
            if gene == "KRAS":
                chrom, pos, ref, alt = kras_by_driver_count[driver_count]
            else:
                chrom = str(13 + driver_offset)
                pos = 200000 + patient_index * 100 + driver_offset
                ref, alt = ("A", "G") if driver_offset % 2 else ("G", "T")
            add_record(gene, chrom, pos, ref, alt, patient_index)

        for gene_index, gene in enumerate(generic_genes, 1):
            if (patient_index + gene_index) % 4 in (0, 1):
                add_record(
                    gene,
                    str(1 + (gene_index % 10)),
                    500000 + gene_index * 1000 + patient_index,
                    "A" if gene_index % 2 else "C",
                    "G" if gene_index % 2 else "T",
                    patient_index,
                )

    # Variants present only in records that downstream cohort filters must remove.
    for patient_index in range(21, 25):
        add_record(
            f"SYN_FILTER_ONLY_{patient_index}",
            "20",
            900000 + patient_index,
            "A",
            "T",
            patient_index,
        )
    return rows


def excel_column(index: int) -> str:
    result = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        result = chr(65 + remainder) + result
    return result


def xlsx_cell(ref: str, value: object, style: int = 0) -> str:
    style_attr = f' s="{style}"' if style else ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f'<c r="{ref}"{style_attr}><v>{value}</v></c>'
    text = escape(str(value))
    preserve = ' xml:space="preserve"' if text.startswith(" ") or text.endswith(" ") else ""
    return f'<c r="{ref}" t="inlineStr"{style_attr}><is><t{preserve}>{text}</t></is></c>'


def make_sheet_xml() -> bytes:
    all_rows = [CLINICAL_HEADERS] + clinical_rows()
    row_xml: list[str] = []
    for row_index, row in enumerate(all_rows, 1):
        cells: list[str] = []
        for col_index, value in enumerate(row, 1):
            style = 1 if row_index == 1 else (2 if col_index == 28 else (3 if col_index in (29, 30) else 0))
            cells.append(xlsx_cell(f"{excel_column(col_index)}{row_index}", value, style))
        row_xml.append(f'<row r="{row_index}">{"".join(cells)}</row>')
    sheet_data_xml = "".join(row_xml).replace(
        '<row r="1">', '<row r="1" ht="32" customHeight="1">'
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<dimension ref="A1:AF25"/>'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        '<sheetFormatPr defaultRowHeight="15"/>'
        '<cols><col min="1" max="1" width="18" customWidth="1"/>'
        '<col min="2" max="2" width="22" customWidth="1"/>'
        '<col min="3" max="5" width="12" customWidth="1"/>'
        '<col min="6" max="7" width="16" customWidth="1"/>'
        '<col min="8" max="8" width="24" customWidth="1"/>'
        '<col min="9" max="9" width="22" customWidth="1"/>'
        '<col min="10" max="13" width="18" customWidth="1"/>'
        '<col min="14" max="18" width="11" customWidth="1"/>'
        '<col min="19" max="20" width="12" customWidth="1"/>'
        '<col min="21" max="28" width="11" customWidth="1"/>'
        '<col min="29" max="30" width="12" customWidth="1"/>'
        '<col min="31" max="32" width="18" customWidth="1"/></cols>'
        f'<sheetData>{sheet_data_xml}</sheetData>'
        '<autoFilter ref="A1:AF25"/>'
        '</worksheet>'
    ).encode("utf-8")


def make_xlsx() -> bytes:
    files = {
        "[Content_Types].xml": b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/><Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/><Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/></Types>''',
        "_rels/.rels": b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/><Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/></Relationships>''',
        "docProps/app.xml": b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>PAAD synthetic example generator</Application></Properties>''',
        "docProps/core.xml": b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>Fully synthetic PAAD example input</dc:title><dc:creator>PAAD project</dc:creator><dcterms:created xsi:type="dcterms:W3CDTF">2000-01-01T00:00:00Z</dcterms:created><dcterms:modified xsi:type="dcterms:W3CDTF">2000-01-01T00:00:00Z</dcterms:modified></cp:coreProperties>''',
        "xl/workbook.xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="명단" sheetId="1" r:id="rId1"/></sheets><calcPr calcId="0" fullCalcOnLoad="1"/></workbook>'''.encode("utf-8"),
        "xl/_rels/workbook.xml.rels": b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>''',
        "xl/styles.xml": b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><numFmts count="1"><numFmt numFmtId="164" formatCode="0.00"/></numFmts><fonts count="2"><font><sz val="10"/><name val="Arial"/></font><font><b/><sz val="10"/><color rgb="FF000000"/><name val="Arial"/></font></fonts><fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FFFFFF00"/><bgColor indexed="64"/></patternFill></fill></fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="4"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment horizontal="center" vertical="center" wrapText="1"/></xf><xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/><xf numFmtId="49" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/></cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>''',
        "xl/worksheets/sheet1.xml": make_sheet_xml(),
    }
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, content in files.items():
            info = zipfile.ZipInfo(name, date_time=(2000, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, content, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return output.getvalue()


def tabular_bytes(headers: list[str], rows: list[list[object]], encoding: str, comments: list[str] | None = None) -> bytes:
    lines = list(comments or [])
    lines.append("\t".join(headers))
    for row in rows:
        if len(row) != len(headers):
            raise ValueError("row width does not match header width")
        values = ["" if value is None else str(value) for value in row]
        if any("\t" in value or "\r" in value or "\n" in value for value in values):
            raise ValueError("tabular values must not contain tabs or line breaks")
        lines.append("\t".join(values))
    return ("\r\n".join(lines) + "\r\n").encode(encoding)


def build_outputs() -> dict[str, bytes]:
    return {
        CLINICAL_NAME: make_xlsx(),
        COHORT_NAME: tabular_bytes(COHORT_HEADERS, cohort_rows(), "cp949"),
        MAF_NAME: tabular_bytes(
            MAF_HEADERS,
            maf_rows(),
            "utf-8",
            comments=[
                "#version 2.4",
                "#synthetic_example=true; contains no patient or study result data",
            ],
        ),
    }


def xlsx_headers(data: bytes) -> tuple[str, list[str], str | None, str | None]:
    namespaces = {
        "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet_name = workbook.find("m:sheets/m:sheet", namespaces).attrib["name"]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    first_row = sheet.find("m:sheetData/m:row", namespaces)
    headers = []
    for cell in first_row.findall("m:c", namespaces):
        text = cell.find("m:is/m:t", namespaces)
        headers.append(text.text if text is not None else "")
    pane = sheet.find("m:sheetViews/m:sheetView/m:pane", namespaces)
    auto_filter = sheet.find("m:autoFilter", namespaces)
    return sheet_name, headers, pane.attrib.get("topLeftCell") if pane is not None else None, auto_filter.attrib.get("ref") if auto_filter is not None else None


def validate_outputs(outputs: dict[str, bytes]) -> None:
    sheet_name, headers, top_left, filter_ref = xlsx_headers(outputs[CLINICAL_NAME])
    if sheet_name != "명단" or headers != CLINICAL_HEADERS:
        raise ValueError("clinical workbook sheet or header mismatch")
    if top_left != "A2" or filter_ref != "A1:AF25":
        raise ValueError("clinical workbook freeze pane or auto-filter mismatch")

    cohort = outputs[COHORT_NAME]
    if not cohort.endswith(b"\r\n") or cohort.count(b"\r\n") != 25:
        raise ValueError("cohort TXT must use CRLF for all 25 lines")
    cohort_text = cohort.decode("cp949")
    if cohort_text.splitlines()[0].split("\t") != COHORT_HEADERS:
        raise ValueError("cohort TXT header mismatch")
    if "합성예시" not in cohort_text:
        raise ValueError("cohort TXT must contain the CP949 encoding sentinel")
    try:
        cohort.decode("utf-8")
    except UnicodeDecodeError:
        pass
    else:
        raise ValueError("cohort TXT unexpectedly decodes as UTF-8")

    maf = outputs[MAF_NAME]
    if not maf.startswith(b"#version 2.4\r\n#synthetic_example=true;"):
        raise ValueError("MAF synthetic comment is missing")
    maf_lines = maf.decode("utf-8").splitlines()
    maf_header = next(line for line in maf_lines if line and not line.startswith("#"))
    if maf_header.split("\t") != MAF_HEADERS:
        raise ValueError("MAF header mismatch")

    combined_text = cohort_text + "\n" + maf.decode("utf-8")
    identifiers = re.findall(r"SYN-(?:TUMOR|NORMAL|CLIN|PATIENT)-[A-Z-]*\d{3}", combined_text)
    if not identifiers or any(not identifier.startswith("SYN-") for identifier in identifiers):
        raise ValueError("synthetic identifier contract failed")
    if "/Users/" in combined_text or "PAAD_tmp" in combined_text:
        raise ValueError("local or source path leaked into generated data")

    clinical = clinical_rows()
    if len(clinical) != 24:
        raise ValueError("clinical fixture must contain 24 synthetic records")
    classifications = [row[CLINICAL_HEADERS.index("Classification")] for row in clinical]
    if classifications.count("PDAC_SYNTHETIC") != 20:
        raise ValueError("clinical fixture must contain 20 PDAC smoke-test records")
    if classifications.count("OTHER_SYNTHETIC") != 2 or classifications.count("Exception") != 2:
        raise ValueError("clinical fixture filter-control records are incomplete")

    cohort_records = cohort_rows()
    sample_index = COHORT_HEADERS.index("Tumor_Sample_Barcode")
    survive_index = COHORT_HEADERS.index("survive")
    recur_index = COHORT_HEADERS.index("Recur")
    if {row[sample_index] for row in cohort_records} != {
        row[CLINICAL_HEADERS.index("Tumor_Sample_Barcode")] for row in clinical
    }:
        raise ValueError("clinical and cohort sample keys do not match")
    for start in range(0, 20, 4):
        block = cohort_records[start : start + 4]
        if {row[survive_index] for row in block} != {0, 1}:
            raise ValueError("each driver-count block needs both OS event states")
        if {row[recur_index] for row in block} != {0, 1}:
            raise ValueError("each driver-count block needs both RFS event states")

    maf_records = maf_rows()
    maf_sample_index = MAF_HEADERS.index("Tumor_Sample_Barcode")
    gene_index = MAF_HEADERS.index("Hugo_Symbol")
    position_index = MAF_HEADERS.index("Start_Position")
    ref_index = MAF_HEADERS.index("Reference_Allele")
    alt_index = MAF_HEADERS.index("Tumor_Seq_Allele2")
    driver_genes = {"KRAS", "TP53", "SMAD4", "CDKN2A"}
    for patient_index in range(1, 21):
        sample = f"SYN-TUMOR-{patient_index:03d}"
        observed = {
            row[gene_index]
            for row in maf_records
            if row[maf_sample_index] == sample and row[gene_index] in driver_genes
        }
        if len(observed) != (patient_index - 1) // 4:
            raise ValueError("driver-count smoke-test groups are not balanced")
    pdac_genes = {
        row[gene_index]
        for row in maf_records
        if int(str(row[maf_sample_index]).rsplit("-", 1)[1]) <= 20
    }
    if len(pdac_genes) < 25:
        raise ValueError("MAF needs at least 25 genes in the PDAC smoke-test cohort")
    kras_loci = {
        (25398284, "C", "T"),
        (25398284, "C", "A"),
        (25398285, "C", "G"),
        (25380276, "T", "G"),
    }
    locus_counts = {
        locus: sum(
            1
            for row in maf_records
            if row[gene_index] == "KRAS"
            and (row[position_index], row[ref_index], row[alt_index]) == locus
        )
        for locus in kras_loci
    }
    if set(locus_counts.values()) != {4}:
        raise ValueError("KRAS fallback subtype loci must each have four records")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare existing files with deterministic generator output without writing",
    )
    args = parser.parse_args()

    outputs = build_outputs()
    validate_outputs(outputs)

    if args.check:
        failures = []
        for name, expected in outputs.items():
            path = OUTPUT_DIR / name
            if not path.is_file():
                failures.append(f"missing: {name}")
            elif path.read_bytes() != expected:
                failures.append(f"content differs: {name}")
        if failures:
            print("\n".join(failures), file=sys.stderr)
            return 1
    else:
        for name, content in outputs.items():
            (OUTPUT_DIR / name).write_bytes(content)

    for name, content in outputs.items():
        print(f"{hashlib.sha256(content).hexdigest()}  {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
