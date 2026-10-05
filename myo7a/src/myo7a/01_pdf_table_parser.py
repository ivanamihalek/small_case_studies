#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = ["pdfplumber>=0.11"]
# ///
"""
Extract the rows mentioning a given gene from a numbered table in a PDF.

Usage:
    uv run 01_pdf_table_parser.py <pdf> <output_dir> <gene> <table_id> [--force]

Example:
    uv run 01_pdf_table_parser.py data/paper_supp.pdf out/ MYO7A S3

The output is a tab-separated file <output_dir>/<pdf stem>.tsv holding the
table header followed by every row whose gene column holds the gene name as a
whole token (case-insensitive).

The table is located by its caption ("Table S3", "Supplementary Table S3", ...)
at the start of a text line, and its column layout must be known in advance
(TABLE_HEADERS). Extraction starts on the first page with the caption and the
expected header, and continues over subsequent pages as long as they repeat
the header.
"""

import argparse
import csv
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

PDF_MAGIC = b"%PDF-"
# PDF spec allows junk before the header; readers tolerate it within the first 1 KB.
PDF_MAGIC_SEARCH_BYTES = 1024
# characters that, adjacent to the gene name, make it part of a longer identifier
# (MYO7A must not match MYO7AB or MYO7A-AS1)
GENE_TOKEN_CHARS = r"A-Za-z0-9_\-"
TABLE_ID_PATTERN = re.compile(r"^[A-Za-z]{0,3}\d+[A-Za-z]?$")
GENE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")
ROW_NUMBER_PATTERN = re.compile(r"^\d+$")

# geometry, in PDF points
LINE_TOLERANCE = 0.5  # characters whose tops differ by less than this are on the same line
FRAGMENT_GAP = 1.0  # a wider horizontal gap between characters separates two cells
CLUSTER_TOLERANCE = 1.5  # fragment left edges this close belong to the same column
HEADER_LINE_GAP = 6.0  # max vertical distance between consecutive header lines

# Header of each known table. A column maps to a word of its label as printed in the PDF;
# repeated words ("HGVS" under each allele) are matched left to right.
# TODO: generalize - these are the tables of sloan-heggen.suppl1.pdf
TABLE_HEADERS = {
    "S3": {
        "Patient #": "Patient",
        "OtoSCOPE version": "OtoSCOPE",
        "Age (range)": "Age",
        "Gene": "Gene",
        "Allele #1": {"chromosomal location": "chromosomal", "HGVS variant": "HGVS", "zygosity": "zygosity"},
        "Allele #2": {"chromosomal location": "chromosomal", "HGVS variant": "HGVS", "zygosity": "zygosity"},
        "Other variants": {"chromosomal location": "chromosomal", "HGVS variant": "HGVS", "zygosity": "zygosity"},
        "Provided Diagnosis": "Provided",
        "Clinical information": {
            "sex (if pertinent)": "sex",
            "Reported inheritance": "Reported",
            "Onset": "Onset",
            "Minimum severity": "Minimum",
            "Symmetry/laterality": "Symmetry/",
            "Physical exam": "Physical",
        },
    },
}
GENE_COLUMN = "Gene"


class InputError(Exception):
    """Raised when an input does not satisfy the script's assumptions."""


def die(message: str, exit_code: int = 1) -> None:
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(exit_code)


def warn(message: str) -> None:
    print(f"Warning: {message}", file=sys.stderr)


#######################################
# input checks
def check_pdf_path(raw_path: str) -> Path:
    path = Path(raw_path).expanduser()
    if not path.exists():
        raise InputError(f"input file '{path}' does not exist")
    if not path.is_file():
        raise InputError(f"input path '{path}' is not a regular file")
    if path.suffix.lower() != ".pdf":
        raise InputError(f"input file '{path}' does not have a .pdf extension")
    if not os.access(path, os.R_OK):
        raise InputError(f"input file '{path}' is not readable")
    if path.stat().st_size == 0:
        raise InputError(f"input file '{path}' is empty")
    try:
        with path.open("rb") as handle:
            head = handle.read(PDF_MAGIC_SEARCH_BYTES)
    except OSError as exc:
        raise InputError(f"cannot read input file '{path}': {exc}") from exc
    if PDF_MAGIC not in head:
        raise InputError(f"input file '{path}' does not look like a PDF (no '%PDF-' header)")
    return path


def check_output_dir(raw_path: str) -> Path:
    path = Path(raw_path).expanduser()
    if not path.exists():
        raise InputError(f"output directory '{path}' does not exist")
    if not path.is_dir():
        raise InputError(f"output path '{path}' is not a directory")
    if not os.access(path, os.W_OK | os.X_OK):
        raise InputError(f"output directory '{path}' is not writable")
    return path


def check_gene_name(raw_name: str) -> str:
    name = raw_name.strip()
    if not name:
        raise InputError("gene name is empty")
    if not GENE_NAME_PATTERN.match(name):
        raise InputError(f"gene name '{raw_name}' contains unexpected characters")
    return name


def check_table_id(raw_id: str) -> str:
    # tolerate "Table S3" or "S3." as well as the bare "S3"
    table_id = re.sub(r"^\s*(supplementary\s+)?table\s*", "", raw_id, flags=re.IGNORECASE)
    table_id = table_id.strip().rstrip(".:")
    if not table_id:
        raise InputError(f"table id '{raw_id}' is empty")
    if not TABLE_ID_PATTERN.match(table_id):
        raise InputError(f"table id '{raw_id}' does not look like a table id (e.g. '3', 'S3', 'S3a')")
    return table_id


def check_table_layout(table_id: str) -> dict:
    header = next((h for key, h in TABLE_HEADERS.items() if key.lower() == table_id.lower()), None)
    if header is None:
        raise InputError(f"no header layout defined for table '{table_id}' (known: {', '.join(TABLE_HEADERS)})")
    if GENE_COLUMN not in header:
        raise InputError(f"header layout of table '{table_id}' has no '{GENE_COLUMN}' column")
    return header


#######################################
# table location and extraction
#
# The tables are borderless, multi-line cells are vertically centred on their row, and
# column positions shift from page to page. Therefore a table is rebuilt from character
# coordinates: on each page the column left edges are found by clustering the left edges
# of the text fragments, each column is identified through its header label, each
# numbered first-column entry starts a row, and every text line goes to the nearest row.


@dataclass
class Fragment:
    """A run of characters on one line with no gap wider than FRAGMENT_GAP (a cell's text on that line)."""
    x0: float
    x1: float
    top: float
    bottom: float
    text: str


@dataclass
class Column:
    name: str
    header_word: str
    header_x0: float = 0.0
    start: float = 0.0


def caption_regex(table_id: str) -> re.Pattern:
    """Matches a table caption at the start of a line."""
    return re.compile(
        rf"^\s*(?:supplementa(?:ry|l)\s+)?table\s*({re.escape(table_id)})(?![A-Za-z0-9])",
        flags=re.IGNORECASE | re.MULTILINE,
    )


def page_text(page) -> str:
    try:
        return page.extract_text() or ""
    except Exception as exc:  # pdfminer can fail on malformed content streams
        warn(f"could not extract text from page {page.page_number}: {exc}")
        return ""


def flatten_header(header: dict) -> list[Column]:
    columns = []
    for name, value in header.items():
        if isinstance(value, dict):
            columns += [Column(f"{name}: {sub_name}", word) for sub_name, word in value.items()]
        else:
            columns.append(Column(name, value))
    return columns


def header_vocabulary(header: dict) -> set[str]:
    words = set()
    for name, value in header.items():
        words.update(name.lower().split())
        if isinstance(value, dict):
            for sub_name, word in value.items():
                words.update(sub_name.lower().split())
                words.add(word.lower())
        else:
            words.add(value.lower())
    return words


def page_lines(page) -> list[list[Fragment]]:
    """Text lines of a page (top to bottom), each a left-to-right list of fragments."""
    chars = sorted((c for c in page.chars if c.get("upright", True)), key=lambda c: (c["top"], c["x0"]))
    grouped: list[list[dict]] = []
    for char in chars:
        if grouped and char["top"] - grouped[-1][0]["top"] <= LINE_TOLERANCE:
            grouped[-1].append(char)
        else:
            grouped.append([char])

    lines = []
    for line_chars in grouped:
        line_chars.sort(key=lambda c: c["x0"])
        fragments: list[Fragment] = []
        current: list[dict] = []
        for char in line_chars:
            if current and char["x0"] - current[-1]["x1"] > FRAGMENT_GAP:
                fragments.append(make_fragment(current))
                current = []
            current.append(char)
        if current:
            fragments.append(make_fragment(current))
        fragments = [f for f in fragments if f.text]
        if fragments:
            lines.append(fragments)
    return lines


def make_fragment(chars: list[dict]) -> Fragment:
    text = re.sub(r"\s+", " ", "".join(c["text"] for c in chars)).strip()
    visible = [c for c in chars if not c["text"].isspace()] or chars
    return Fragment(
        x0=visible[0]["x0"],
        x1=visible[-1]["x1"],
        top=min(c["top"] for c in chars),
        bottom=max(c["bottom"] for c in chars),
        text=text,
    )


def locate_header(lines: list[list[Fragment]], header: dict) -> tuple[list[Column], float] | None:
    """Find the header labels in left-to-right order; returns the columns with their header
    positions and the bottom of the header block, or None if this page does not carry the header."""
    columns = flatten_header(header)
    vocabulary = header_vocabulary(header)
    words = sorted(
        (
            (word, fragment.x0 + offset * (fragment.x1 - fragment.x0) / max(len(fragment.text), 1), fragment)
            for line in lines
            for fragment in line
            for word, offset in split_with_offsets(fragment.text)
        ),
        key=lambda item: item[1],
    )
    previous_x0 = float("-inf")
    header_fragments = []
    for column in columns:
        match = next((w for w in words if w[0] == column.header_word and w[1] > previous_x0), None)
        if match is None:
            return None
        column.header_x0 = previous_x0 = match[1]
        header_fragments.append(match[2])
    # the header block may continue below the label words ("location", "exam"): extend it over
    # the following lines as long as they consist mostly of header vocabulary
    header_bottom = max(f.bottom for f in header_fragments)
    for line in lines:
        if line[0].top <= header_bottom:
            continue
        line_words = " ".join(f.text for f in line).lower().split()
        if line[0].top - header_bottom < HEADER_LINE_GAP and sum(w in vocabulary for w in line_words) * 2 >= len(line_words):
            header_bottom = max(f.bottom for f in line)
        else:
            break
    return columns, header_bottom


def split_with_offsets(text: str) -> list[tuple[str, int]]:
    return [(m.group(0), m.start()) for m in re.finditer(r"\S+", text)]


def assign_column_starts(columns: list[Column], body: list[list[Fragment]]) -> None:
    """Column left edges = the most frequent fragment left edge between the previous column's
    start and the column's (centred) header label. A column empty on this page falls back to
    its header position, which is harmless since no fragment needs to land in it."""
    clusters: list[list[float]] = []
    for x0 in sorted(f.x0 for line in body for f in line):
        if clusters and x0 - clusters[-1][0] <= CLUSTER_TOLERANCE:
            clusters[-1].append(x0)
        else:
            clusters.append([x0])
    previous_start = float("-inf")
    for column in columns:
        candidates = [
            c for c in clusters if previous_start + CLUSTER_TOLERANCE < c[0] <= column.header_x0 + CLUSTER_TOLERANCE
        ]
        if candidates:
            best = max(candidates, key=lambda c: (len(c), c[0]))
            column.start = best[0]
        else:
            column.start = max(column.header_x0, previous_start + CLUSTER_TOLERANCE)
        previous_start = column.start


def column_index(columns: list[Column], x0: float) -> int:
    index = 0
    for i, column in enumerate(columns):
        if x0 >= column.start - CLUSTER_TOLERANCE:
            index = i
    return index


def join_cell(parts: list[str]) -> str:
    text = ""
    for part in parts:
        if not text:
            text = part
        elif re.search(r"[A-Za-z]-$", text) or (re.search(r"[ACGT]$", text) and re.match(r"[ACGT]{2,}", part)):
            # "severe-" + "profound", or a nucleotide sequence wrapped onto the next line
            text += part
        else:
            text += " " + part
    return text


def page_rows(page, header: dict) -> list[list[str]] | None:
    """Rows of the table on this page, or None if the page does not carry the table header."""
    lines = page_lines(page)
    located = locate_header(lines, header)
    if located is None:
        return None
    columns, header_bottom = located
    second_column_x0 = columns[1].header_x0

    # body: below the header, up to the first line whose first-column text is not a row number
    # (footnotes, page footer)
    body: list[list[Fragment]] = []
    for line in lines:
        if line[0].top <= header_bottom:
            continue
        first = line[0]
        if first.x0 < second_column_x0 and not ROW_NUMBER_PATTERN.match(first.text):
            break
        body.append(line)
    assign_column_starts(columns, body)

    anchors = [
        (line[0].top + line[0].bottom) / 2
        for line in body
        if column_index(columns, line[0].x0) == 0 and ROW_NUMBER_PATTERN.match(line[0].text)
    ]
    if not anchors:
        warn(f"page {page.page_number}: table header found but no numbered rows")
        return []

    cells = [[[] for _ in columns] for _ in anchors]
    for line in body:
        middle = (line[0].top + line[0].bottom) / 2
        row = min(range(len(anchors)), key=lambda i: abs(anchors[i] - middle))
        for fragment in line:
            cells[row][column_index(columns, fragment.x0)].append(fragment.text)
    return [[join_cell(parts) for parts in row] for row in cells]


def collect_table_rows(pdf, table_id: str, header: dict) -> list[list[str]]:
    own_caption = caption_regex(table_id)
    pages = pdf.pages
    if not pages:
        raise InputError("the PDF has no pages")

    start_index = None
    for index, page in enumerate(pages):
        if own_caption.search(page_text(page)) and page_rows(page, header) is not None:
            start_index = index
            break
    if start_index is None:
        raise InputError(f"no page with caption 'Table {table_id}' followed by the expected header was found")

    rows: list[list[str]] = [[column.name for column in flatten_header(header)]]
    last_index = start_index
    for index in range(start_index, len(pages)):
        page_data = page_rows(pages[index], header)
        if page_data is None:
            break
        rows += page_data
        last_index = index
    check_row_numbers(rows[1:])
    print(f"Table {table_id}: pages {start_index + 1}-{last_index + 1}, {len(rows) - 1} data rows", file=sys.stderr)
    return rows


def check_row_numbers(rows: list[list[str]]) -> None:
    numbers = [int(row[0]) for row in rows]
    for previous, current in zip(numbers, numbers[1:]):
        if current != previous + 1:
            warn(f"row numbers jump from {previous} to {current} - rows may be missing or merged")


def select_gene_rows(rows: list[list[str]], gene: str, gene_column: str) -> list[list[str]]:
    gene_regex = re.compile(
        rf"(?<![{GENE_TOKEN_CHARS}]){re.escape(gene)}(?![{GENE_TOKEN_CHARS}])",
        flags=re.IGNORECASE,
    )
    header, data = rows[0], rows[1:]
    index = header.index(gene_column)
    return [header] + [row for row in data if gene_regex.search(row[index])]


#######################################
# output
def write_tsv(rows: list[list[str]], out_path: Path) -> None:
    # write to a temp file in the same directory, then rename: no half-written output on failure
    fd, tmp_name = tempfile.mkstemp(dir=out_path.parent, prefix=f".{out_path.stem}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
            writer.writerows(rows)
        os.replace(tmp_name, out_path)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract rows mentioning a gene from a numbered table in a PDF into a .tsv file."
    )
    parser.add_argument("pdf", help="path to the input PDF")
    parser.add_argument("output_dir", help="existing directory to write <pdf stem>.tsv into")
    parser.add_argument("gene", help="gene name to look for, e.g. MYO7A")
    parser.add_argument("table_id", help="table id as in the caption, e.g. 3 or S3")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite an existing output file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        pdf_path = check_pdf_path(args.pdf)
        out_dir = check_output_dir(args.output_dir)
        gene = check_gene_name(args.gene)
        table_id = check_table_id(args.table_id)
        header = check_table_layout(table_id)
    except InputError as exc:
        die(str(exc))

    out_path = out_dir / f"{pdf_path.stem}.tsv"
    if out_path.exists() and not args.force:
        die(f"output file '{out_path}' already exists (use --force to overwrite)")
    if out_path.exists() and not out_path.is_file():
        die(f"output path '{out_path}' exists and is not a regular file")

    try:
        import pdfplumber
        from pdfminer.pdfparser import PDFSyntaxError
    except ImportError:
        die("pdfplumber is not installed; run with 'uv run' or 'pip install pdfplumber'")

    try:
        with pdfplumber.open(pdf_path) as pdf:
            rows = collect_table_rows(pdf, table_id, header)
    except InputError as exc:
        die(str(exc))
    except PDFSyntaxError as exc:
        die(f"'{pdf_path}' is not a valid PDF: {exc}")
    except Exception as exc:  # encrypted, truncated, or otherwise unreadable PDFs
        die(f"could not read '{pdf_path}': {type(exc).__name__}: {exc}")

    selected = select_gene_rows(rows, gene, GENE_COLUMN)
    if len(selected) == 1:
        die(f"no rows mentioning '{gene}' found in Table {table_id}", exit_code=2)

    try:
        write_tsv(selected, out_path)
    except OSError as exc:
        die(f"could not write '{out_path}': {exc}")
    print(f"Wrote {len(selected) - 1} rows to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
