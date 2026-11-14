#!/usr/bin/env python3
"""
Find all individuals (IID) carrying a variant in a given gene (HGNC column),
plus all other individuals sharing any variant (VarID, in any gene) with them,
then extract all rows (variants in any gene) for these individuals, sorted by IID.
"""
import argparse
import sys
from pathlib import Path
from typing import Any, List, Optional, Set, Tuple

import openpyxl

IID_COLUMN = "IID"
GENE_COLUMN = "HGNC"
VARID_COLUMN = "VarID"


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


def select_carrier_rows(header: List[Any], data_rows: List[Tuple[Any, ...]],
                        gene: str) -> Tuple[List[Tuple[Any, ...]], Set[Any], Set[Any]]:
    iid_idx = column_index(header, IID_COLUMN)
    gene_idx = column_index(header, GENE_COLUMN)
    varid_idx = column_index(header, VARID_COLUMN)
    carriers = {row[iid_idx] for row in data_rows if row[gene_idx] is not None and str(row[gene_idx]).strip() == gene}
    # all variants (in any gene) found in the carriers, and other individuals sharing any of them
    carrier_varids = {row[varid_idx] for row in data_rows if row[iid_idx] in carriers and row[varid_idx] is not None}
    varid_sharers = {row[iid_idx] for row in data_rows if row[varid_idx] in carrier_varids} - carriers
    selected = carriers | varid_sharers
    selected_rows = [row for row in data_rows if row[iid_idx] in selected]
    # sort is stable, so within an individual the original row order is kept
    return sorted(selected_rows, key=lambda row: str(row[iid_idx])), carriers, varid_sharers


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
                        help="output xlsx table (default: <input>.<gene>_carriers.xlsx)")
    parser.add_argument("-g", "--gene", default="CLPTM1", help="HGNC gene symbol (exact match)")
    parser.add_argument("-s", "--sheet", default=None, help="sheet name (default: first sheet)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_path = args.output or args.input.with_name(f"{args.input.stem}.{args.gene}_carriers.xlsx")

    header, data_rows = read_table(args.input, args.sheet)
    carrier_rows, carriers, varid_sharers = select_carrier_rows(header, data_rows, args.gene)

    summary_msg = (f"{args.input.name}: {len(data_rows)} rows, {len(carriers)} individuals with {args.gene}, "
                   f"{len(varid_sharers)} other individuals sharing a {VARID_COLUMN} with them, "
                   f"{len(carrier_rows)} rows extracted")
    print(summary_msg)
    if not carrier_rows:
        warning_msg = f"No individuals with {args.gene} in the {GENE_COLUMN} column - no output written."
        print(warning_msg, file=sys.stderr)
        return 1

    write_table(output_path, header, carrier_rows, f"{args.gene}_carriers")
    output_msg = f"written to {output_path}"
    print(output_msg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
