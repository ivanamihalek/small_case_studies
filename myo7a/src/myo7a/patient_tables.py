"""
Read patient tables of any supported format into one list of PatientRecord (patient_record.py),
choosing the reader by the file extension: .tsv for sloan_heggen_tsv.py, .xlsx for moorfields_xlsx.py,
.csv for ambrosio_csv.py. With prefix_ids, each patient id is prefixed with the first three characters of
its file name (e.g. "amb 1" for ambrosio_myo7a.csv), to tell apart the patients of different tables.
output_path() names an output file after the same prefixes of all tables (e.g. data/amb_slo_moo.svg).
"""

from dataclasses import replace
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import ambrosio_csv
import moorfields_xlsx
import sloan_heggen_tsv
from patient_record import InputError, PatientRecord, check_input_path

READERS: Dict[str, Callable[[Path], List[PatientRecord]]] = {".tsv": sloan_heggen_tsv.read_patients,
                                                             ".xlsx": moorfields_xlsx.read_patients,
                                                             ".csv": ambrosio_csv.read_patients}
PREFIX_LENGTH = 3


def table_prefix(path: Path) -> str:
    return path.name[:PREFIX_LENGTH]


def output_path(paths: List[Path], suffix: str) -> Path:
    """In the first table's directory: the tables' prefixes, in the order given, joined by '_', plus suffix."""
    prefixes = list(dict.fromkeys(table_prefix(path) for path in paths))
    return paths[0].parent / f"{'_'.join(prefixes)}{suffix}"


def read_patient_tables(raw_paths: List[str], prefix_ids: bool = False) -> Tuple[List[Path], List[PatientRecord]]:
    """The checked paths and the patient records of all tables, in the order given."""
    paths: List[Path] = []
    records: List[PatientRecord] = []
    for raw_path in raw_paths:
        suffix = Path(raw_path).suffix.lower()
        if suffix not in READERS:
            error = f"input file '{raw_path}': unsupported extension (supported: {', '.join(READERS)})"
            raise InputError(error)
        path = check_input_path(raw_path, suffix)
        paths.append(path)
        table_records = READERS[suffix](path)
        if prefix_ids:
            prefix = table_prefix(path)
            table_records = [replace(record, patient=f"{prefix} {record.patient}") for record in table_records]
        records.extend(table_records)
    return paths, records
