"""
Reader for the patient table in the tab-separated format written by 01_pdf_table_parser.py
from the Sloan-Heggen supplementary table: one row per patient, with the columns below.
Empty cells may also be written as '_' or '-'; they are returned as "".
"""

import csv
from pathlib import Path
from typing import Dict, List

from patient_record import InputError, PatientRecord, Phenotype, check_input_path

ALLELE1_COLUMN = "Allele #1: HGVS variant"
ALLELE2_COLUMN = "Allele #2: HGVS variant"
ZYGOSITY1_COLUMN = "Allele #1: zygosity"
OTHER_COLUMN = "Other variants: HGVS variant"
DIAGNOSIS_COLUMN = "Provided Diagnosis"
PATIENT_COLUMN = "Patient #"
GENE_COLUMN = "Gene"
AGE_COLUMN = "Age (range)"
REQUIRED_COLUMNS = [ALLELE1_COLUMN, ALLELE2_COLUMN, ZYGOSITY1_COLUMN, OTHER_COLUMN, DIAGNOSIS_COLUMN, AGE_COLUMN]

HOMOZYGOUS = "hom"
PHENOTYPES = {"autosomal recessive non-syndromic hearing loss": Phenotype.NONSYNDROMIC,
              "usher syndrome 1b": Phenotype.USHER_1B}
EMPTY_CELL = {"", "_", "-"}


def check_tsv_path(raw_path: str) -> Path:
    return check_input_path(raw_path, ".tsv")


def read_rows(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
            if missing:
                error = f"'{path}' lacks the column(s): {', '.join(missing)}"
                raise InputError(error)
            return [{key: (value or "").strip() for key, value in row.items() if key} for row in reader]
    except UnicodeDecodeError as exc:
        error = f"'{path}' is not UTF-8 text: {exc}"
        raise InputError(error) from exc
    except csv.Error as exc:
        error = f"'{path}' is not a valid TSV file: {exc}"
        raise InputError(error) from exc


def cell(row: Dict[str, str], column: str) -> str:
    value = row.get(column, "")
    return "" if value in EMPTY_CELL else value


def patient_record(row: Dict[str, str]) -> PatientRecord:
    diagnosis = cell(row, DIAGNOSIS_COLUMN)
    return PatientRecord(patient=row.get(PATIENT_COLUMN, "?"), gene=cell(row, GENE_COLUMN),
                         allele1=cell(row, ALLELE1_COLUMN), allele2=cell(row, ALLELE2_COLUMN),
                         homozygous=cell(row, ZYGOSITY1_COLUMN).lower() == HOMOZYGOUS,
                         other_variants=cell(row, OTHER_COLUMN), diagnosis=diagnosis,
                         phenotype=PHENOTYPES.get(diagnosis.lower()), age=cell(row, AGE_COLUMN))


def read_patients(path: Path) -> List[PatientRecord]:
    """The patient records of a .tsv file that passed check_tsv_path()."""
    return [patient_record(row) for row in read_rows(path)]
