"""
Reader for the Moorfields MYO7A patient table (data/moorfield.xlsx, sheet "Genetic results"):
two header rows (the group "Variant 1", "Variant 2", "Variant 3" over the columns "c.DNA" and
"Protein"), one row per patient, then footnotes in the first column only.

All patients have Usher syndrome type 1. Those marked '*' after the patient ID have a mild
phenotype and are given the auditory (non-syndromic) side; those marked '⁺⁺' have sector
retinitis pigmentosa and stay on the Usher side, with the note in the diagnosis text. The
table gives no age, so the age is always unknown.

A patient is homozygous if the 'Homozygous (Y/N)' column says so, or if the two c.DNA changes
are identical; the second allele is then left empty, as patient_record.py expects.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional

from openpyxl import load_workbook
from openpyxl.utils.exceptions import InvalidFileException

from patient_record import InputError, PatientRecord, Phenotype, check_input_path

SHEET = "Genetic results"
GENE = "MYO7A"  # the table caption: "Variants in MYO7A identified in the Usher Syndrome Type 1 patients"
PATIENT_COLUMN = "Patient ID"
HOMOZYGOUS_COLUMN = "Homozygous (Y/N)"
CDNA = "c.DNA"
PROTEIN = "Protein"
VARIANTS = ["Variant 1", "Variant 2", "Variant 3"]
REQUIRED_COLUMNS = [PATIENT_COLUMN, HOMOZYGOUS_COLUMN] + [f"{v}: {c}" for v in VARIANTS[:2] for c in (CDNA, PROTEIN)]
HEADER_ROWS = 2

YES = {"y", "yes"}
NO = {"n", "no"}
MILD_MARKER = "*"
SECTOR_RP_MARKER = "⁺⁺"
DIAGNOSIS = "Usher syndrome 1"
MILD_DIAGNOSIS = f"{DIAGNOSIS}, mild phenotype"
SECTOR_RP_DIAGNOSIS = f"{DIAGNOSIS}, sector retinitis pigmentosa"


def check_xlsx_path(raw_path: str) -> Path:
    return check_input_path(raw_path, ".xlsx")


def text(value: object) -> str:
    """A cell value as text, with non-breaking and repeated spaces collapsed and the ends stripped."""
    return " ".join(str(value).split()) if value is not None else ""


def column_names(group_row: List[str], name_row: List[str]) -> List[str]:
    """'Variant 1: c.DNA' etc. under a group header, the plain name elsewhere."""
    names: List[str] = []
    group = ""
    for group_cell, name in zip(group_row, name_row):
        if group_cell:
            group = group_cell
        elif name not in (CDNA, PROTEIN):
            group = ""
        names.append(f"{group}: {name}" if group else name)
    return names


def read_rows(path: Path) -> List[Dict[str, str]]:
    try:
        workbook = load_workbook(path, read_only=True, data_only=True)
    except (InvalidFileException, KeyError, ValueError) as exc:
        error = f"'{path}' is not a valid .xlsx file: {exc}"
        raise InputError(error) from exc
    try:
        if SHEET not in workbook.sheetnames:
            error = f"'{path}' has no sheet '{SHEET}' (sheets: {', '.join(workbook.sheetnames)})"
            raise InputError(error)
        rows = [[text(value) for value in row] for row in workbook[SHEET].iter_rows(values_only=True)]
    finally:
        workbook.close()
    if len(rows) <= HEADER_ROWS:
        error = f"'{path}': sheet '{SHEET}' has no data rows"
        raise InputError(error)

    names = column_names(rows[0], rows[1])
    missing = [column for column in REQUIRED_COLUMNS if column not in names]
    if missing:
        error = f"'{path}' lacks the column(s): {', '.join(missing)}"
        raise InputError(error)
    # footnotes and the caption fill the first column only
    patient_rows = [row for row in rows[HEADER_ROWS:] if any(row[1:])]
    return [dict(zip(names, row)) for row in patient_rows]


def protein_variant(row: Dict[str, str], variant: str) -> str:
    """The protein change of a variant column group, with the spaces removed ('p. (Asp2191ValfsTer4)')."""
    return re.sub(r"\s+", "", row.get(f"{variant}: {PROTEIN}", ""))


def allele(row: Dict[str, str], variant: str) -> str:
    """'c.635G>A, p.(Arg212His)', in the style of the Sloan-Heggen allele cells; "" if not given."""
    parts = [row.get(f"{variant}: {CDNA}", ""), protein_variant(row, variant)]
    return ", ".join(part for part in parts if part)


def is_homozygous(row: Dict[str, str]) -> bool:
    value = row.get(HOMOZYGOUS_COLUMN, "").lower()
    if value not in YES | NO:
        error = f"patient {row.get(PATIENT_COLUMN)}: unexpected '{HOMOZYGOUS_COLUMN}' value '{value}'"
        raise InputError(error)
    cdna1, cdna2 = (row.get(f"{variant}: {CDNA}", "") for variant in VARIANTS[:2])
    return value in YES or (cdna1 != "" and cdna1 == cdna2)


def diagnosis_and_phenotype(raw_id: str) -> tuple[str, Optional[Phenotype]]:
    if raw_id.endswith(MILD_MARKER):
        return MILD_DIAGNOSIS, Phenotype.NONSYNDROMIC
    if raw_id.endswith(SECTOR_RP_MARKER):
        return SECTOR_RP_DIAGNOSIS, Phenotype.USHER_1B
    return DIAGNOSIS, Phenotype.USHER_1B


def patient_record(row: Dict[str, str]) -> PatientRecord:
    raw_id = row.get(PATIENT_COLUMN, "?")
    patient = raw_id.rstrip(MILD_MARKER + SECTOR_RP_MARKER).strip()
    diagnosis, phenotype = diagnosis_and_phenotype(raw_id)
    homozygous = is_homozygous(row)
    allele2 = "" if homozygous else allele(row, VARIANTS[1])
    other_variants = "; ".join(filter(None, (allele(row, variant) for variant in VARIANTS[2:])))
    return PatientRecord(patient=patient, gene=GENE, allele1=allele(row, VARIANTS[0]), allele2=allele2,
                         homozygous=homozygous, other_variants=other_variants, diagnosis=diagnosis,
                         phenotype=phenotype, age="")


def read_patients(path: Path) -> List[PatientRecord]:
    """The patient records of an .xlsx file that passed check_xlsx_path()."""
    return [patient_record(row) for row in read_rows(path)]
