#!/usr/bin/env python3
"""
Find all variants (VarID) shared by more than one affected individual (IID with Pheno == Affected),
then extract the rows of the affected carriers of these variants, sorted by the number of
affected carriers (descending), then by VarID and IID. Columns with the number of affected and
total (any phenotype) carriers of the variant are appended to each row.
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import openpyxl

IID_COLUMN = "IID"
PHENO_COLUMN = "Pheno"
VARID_COLUMN = "VarID"
AFFECTED_VALUE = "Affected"
N_AFFECTED_COLUMN = "N_affected_carriers"
N_CARRIERS_COLUMN = "N_carriers"


def read_table(input_path: Path, sheet_name: Optional[str]) -> Tuple[List[Any], List[Tuple[Any, ...]]]:
    workbook = openpyxl.load_workbook(input_path, read_only=True, data_only=True)
    worksheet = workbook[sheet_name] if sheet_name else workbook.worksheets[0]
    rows = [row for row in worksheet.iter_rows(values_only=True) if any(cell is not None for cell in row)]
    workbook.close()
    if not rows:
        error_msg = f"No data found in {input_path}"
        raise ValueError(error_msg)
    return list(rows[0]), rows[1:]


def column_index(header: List[Any], column_name: str) -> int:
    if column_name not in header:
        error_msg = f"Column '{column_name}' not found in the header: {header}"
        raise ValueError(error_msg)
    return header.index(column_name)


def is_affected(pheno: Any) -> bool:
    return pheno is not None and str(pheno).strip() == AFFECTED_VALUE


def select_shared_variant_rows(header: List[Any], data_rows: List[Tuple[Any, ...]],
                               min_affected: int) -> Tuple[List[Tuple[Any, ...]], Dict[Any, int]]:
    iid_idx = column_index(header, IID_COLUMN)
    pheno_idx = column_index(header, PHENO_COLUMN)
    varid_idx = column_index(header, VARID_COLUMN)
    affected_carriers: Dict[Any, Set[Any]] = defaultdict(set)
    all_carriers: Dict[Any, Set[Any]] = defaultdict(set)
    for row in data_rows:
        if row[varid_idx] is None:
            continue
        all_carriers[row[varid_idx]].add(row[iid_idx])
        if is_affected(row[pheno_idx]):
            affected_carriers[row[varid_idx]].add(row[iid_idx])
    n_affected = {varid: len(iids) for varid, iids in affected_carriers.items() if len(iids) >= min_affected}
    # only the rows of affected carriers are kept; N_carriers still counts carriers of any phenotype
    selected_rows = [row + (n_affected[row[varid_idx]], len(all_carriers[row[varid_idx]]))
                     for row in data_rows if row[varid_idx] in n_affected and is_affected(row[pheno_idx])]
    selected_rows.sort(key=lambda row: (-n_affected[row[varid_idx]], str(row[varid_idx]), str(row[iid_idx])))
    return selected_rows, n_affected


def write_table(output_path: Path, header: List[Any], rows: List[Tuple[Any, ...]], sheet_title: str) -> None:
    workbook = openpyxl.Workbook()
    worksheet = workbook.active
    worksheet.title = sheet_title
    worksheet.append(header)
    for row in rows:
        worksheet.append(list(row))
    worksheet.freeze_panes = "A2"
    workbook.save(output_path)


def parse_args() -> argparse.Namespace:
    tables_dir = Path(__file__).resolve().parent.parent / "tables"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, default=tables_dir / "sample.xlsx", help="input xlsx table")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="output xlsx table (default: <input>.shared_variants.xlsx)")
    parser.add_argument("-m", "--min-affected", type=int, default=2,
                        help="minimum number of affected individuals sharing a variant (default: 2)")
    parser.add_argument("-s", "--sheet", default=None, help="sheet name (default: first sheet)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_path = args.output or args.input.with_name(f"{args.input.stem}.shared_variants.xlsx")

    header, data_rows = read_table(args.input, args.sheet)
    shared_rows, n_affected = select_shared_variant_rows(header, data_rows, args.min_affected)

    summary_msg = (f"{args.input.name}: {len(data_rows)} rows, {len(n_affected)} variants shared by "
                   f"at least {args.min_affected} {AFFECTED_VALUE.lower()} individuals, {len(shared_rows)} rows extracted")
    print(summary_msg)
    if not shared_rows:
        warning_msg = "No shared variants found - no output written."
        print(warning_msg, file=sys.stderr)
        return 1

    write_table(output_path, header + [N_AFFECTED_COLUMN, N_CARRIERS_COLUMN], shared_rows, "shared_variants")
    output_msg = f"written to {output_path}"
    print(output_msg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
