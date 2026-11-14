"""
The input-format independent patient record. Each patient table reader (e.g. sloan_heggen_tsv.py)
turns its own input into a list of PatientRecord, which 03_allele_pair_visualization.py and
05_combined_visualization.py then filter and draw; records from different readers can be
concatenated and drawn in the same figure.
"""

import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional


class InputError(Exception):
    """Raised when an input does not satisfy the reader's assumptions."""


class Phenotype(Enum):
    NONSYNDROMIC = "autosomal recessive non-syndromic hearing loss"  # drawn above the protein line
    USHER_1B = "Usher syndrome 1B"  # drawn below the protein line


@dataclass
class PatientRecord:
    patient: str  # patient id, as given in the source
    gene: str  # "" if not given
    allele1: str  # text with the HGVS protein-level (p.) variant(s) of allele #1, "" if not given
    allele2: str  # the same for allele #2; "" if not given (as is normal for a homozygous allele #1)
    homozygous: bool  # allele #1 is homozygous
    other_variants: str  # any further variants of the patient, "" if none
    diagnosis: str  # as given, for labels and notes
    phenotype: Optional[Phenotype]  # None for any diagnosis other than the two drawn
    age: str  # age range in the form "0-10", "" if unknown


def check_input_path(raw_path: str, suffix: str) -> Path:
    """The path of an existing, readable, non-empty file with the given extension (e.g. ".tsv")."""
    path = Path(raw_path).expanduser()
    if not path.exists():
        error = f"input file '{path}' does not exist"
        raise InputError(error)
    if not path.is_file():
        error = f"input path '{path}' is not a regular file"
        raise InputError(error)
    if path.suffix.lower() != suffix:
        error = f"input file '{path}' does not have a {suffix} extension"
        raise InputError(error)
    if not os.access(path, os.R_OK):
        error = f"input file '{path}' is not readable"
        raise InputError(error)
    if path.stat().st_size == 0:
        error = f"input file '{path}' is empty"
        raise InputError(error)
    return path
