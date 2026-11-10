"""
Reader for the MYO7A patients of Ambrosio et al., Doc Ophthalmol (2021) 143:39-51, Table 1
(data/ambrosio_myo7a.csv; see doc/ambrosio_table1_extraction.md). Despite the .csv extension the
file is tab-separated: one row per patient, with the columns below.

All the MYO7A patients of the table have Usher syndrome type 1. The protein changes are given in
one-letter amino acid code (p.G214R, p.R1701*, p.A1277Qfs*122); they are rewritten in three-letter
code (p.Gly214Arg, p.Arg1701Ter, p.Ala1277GlnfsTer122), as 03_allele_pair_visualization.py expects,
and a synonymous change (p.E968E) as p.Glu968=. The age at the first examination, in years, is
turned into the age ranges of the Sloan-Heggen table ("0-10", "11-20", ...).

A patient is homozygous if the two c.DNA changes are identical; the second allele is then left
empty, as patient_record.py expects.
"""

import csv
import math
import re
from pathlib import Path
from typing import Dict, List

from patient_record import InputError, PatientRecord, Phenotype, check_input_path

PATIENT_COLUMN = "patient"
GENE_COLUMN = "gene"
CDNA1_COLUMN = "mutation1_cdna"
PROTEIN1_COLUMN = "mutation1_protein"
CDNA2_COLUMN = "mutation2_cdna"
PROTEIN2_COLUMN = "mutation2_protein"
AGE_COLUMN = "age_first_exam_years"
REQUIRED_COLUMNS = [PATIENT_COLUMN, GENE_COLUMN, CDNA1_COLUMN, PROTEIN1_COLUMN, CDNA2_COLUMN, PROTEIN2_COLUMN,
                    AGE_COLUMN]

DIAGNOSIS = "Usher syndrome 1"
AGE_RANGE_WIDTH = 10

THREE_LETTER = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
                "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
                "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Ter"}
ALREADY_THREE_LETTER = re.compile(r"^p\.\(?[A-Z][a-z]{2}\d")
SYNONYMOUS = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)\1$")


def check_csv_path(raw_path: str) -> Path:
    return check_input_path(raw_path, ".csv")


def read_rows(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
            if missing:
                error = f"'{path}' lacks the column(s): {', '.join(missing)} (is it tab-separated?)"
                raise InputError(error)
            return [{key: (value or "").strip() for key, value in row.items() if key} for row in reader]
    except UnicodeDecodeError as exc:
        error = f"'{path}' is not UTF-8 text: {exc}"
        raise InputError(error) from exc
    except csv.Error as exc:
        error = f"'{path}' is not a valid tab-separated file: {exc}"
        raise InputError(error) from exc


def three_letter_protein(change: str) -> str:
    """'p.G214R' -> 'p.Gly214Arg', 'p.A1277Qfs*122' -> 'p.Ala1277GlnfsTer122', 'p.E968E' -> 'p.Glu968='."""
    if not change or ALREADY_THREE_LETTER.match(change):
        return change
    if not change.startswith("p."):
        error = f"protein change '{change}' does not start with 'p.'"
        raise InputError(error)

    def residue(match: re.Match) -> str:
        code = match.group(0)
        if code not in THREE_LETTER:
            error = f"protein change '{change}': unknown amino acid code '{code}'"
            raise InputError(error)
        return THREE_LETTER[code]

    # the keywords (fs, del, delins, ins, dup) are lower case, so every capital letter or '*' is a residue
    converted = "p." + re.sub(r"[A-Z*]", residue, change[2:])
    synonymous = SYNONYMOUS.match(converted)
    return f"p.{synonymous.group(1)}{synonymous.group(2)}=" if synonymous else converted


def allele(cdna: str, protein: str) -> str:
    """'c.640G>A, p.Gly214Arg', in the style of the Sloan-Heggen allele cells; "" if not given."""
    return ", ".join(part for part in (cdna, three_letter_protein(protein)) if part)


def age_range(raw_age: str, patient: str) -> str:
    """'0-10' for an age of up to 10 years, '11-20' for above 10 up to 20, etc.; "" if not given."""
    if not raw_age:
        return ""
    try:
        age = float(raw_age)
    except ValueError as exc:
        error = f"patient {patient}: '{AGE_COLUMN}' value '{raw_age}' is not a number"
        raise InputError(error) from exc
    if age < 0:
        error = f"patient {patient}: negative '{AGE_COLUMN}' value {age:.2f}"
        raise InputError(error)
    upper = max(1, math.ceil(age / AGE_RANGE_WIDTH)) * AGE_RANGE_WIDTH
    lower = 0 if upper == AGE_RANGE_WIDTH else upper - AGE_RANGE_WIDTH + 1
    return f"{lower}-{upper}"


def patient_record(row: Dict[str, str]) -> PatientRecord:
    patient = row.get(PATIENT_COLUMN, "?")
    cdna1, cdna2 = row.get(CDNA1_COLUMN, ""), row.get(CDNA2_COLUMN, "")
    homozygous = cdna1 != "" and cdna1 == cdna2
    allele1 = allele(cdna1, row.get(PROTEIN1_COLUMN, ""))
    allele2 = "" if homozygous else allele(cdna2, row.get(PROTEIN2_COLUMN, ""))
    return PatientRecord(patient=patient, gene=row.get(GENE_COLUMN, ""), allele1=allele1, allele2=allele2,
                         homozygous=homozygous, other_variants="", diagnosis=DIAGNOSIS, phenotype=Phenotype.USHER_1B,
                         age=age_range(row.get(AGE_COLUMN, ""), patient))


def read_patients(path: Path) -> List[PatientRecord]:
    """The patient records of a .csv file that passed check_csv_path()."""
    return [patient_record(row) for row in read_rows(path)]
