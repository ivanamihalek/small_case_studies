#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = ["openpyxl"]
# ///
"""
Permissive version of 03_allele_pair_visualization.py: draw the protein-level allele pairs
of compound heterozygous / homozygous patients, collected from any number of .tsv and .xlsx
tables, as arcs over a protein sequence line.

Usage:
    uv run 04_allele_pair_visualization_permissive.py <table> [<table> ...] [-o OUT.svg]
        [--max-age AGE] [--mark POS1 POS2 SIDE] [--force]

Example:
    uv run 04_allele_pair_visualization_permissive.py data/sloan-heggen.suppl1.tsv \
        data/ush1B_lit_review_genotypes.xlsx -o allele_pairs.svg --max-age 10

Each table (every sheet, for .xlsx) must have a header row - not necessarily the first row -
with a pubmed, a patient id, an allele 1, an allele 2, an age and a phenotype column. Column
names are matched case-insensitively against the equivalence classes defined below
(*_COLUMNS). Age-of-onset columns are not age columns and are ignored. Optional columns:
allele 1 zygosity, other variants, gene.

The age, if given (a number such as 34, or a range such as "51-60" or "5-6 yrs"), is stored
as the decade of life it falls in ("0-10", "11-20", "21-30", ...); for a range, the upper
end decides. A patient's arc is drawn faded (reduced opacity) when the age is unknown or above
--max-age (default 10).

The phenotype must fall into one of two classes (see USH1B_PHENOTYPES and
NSHL_PHENOTYPE_PATTERNS): Usher syndrome 1B, drawn below the protein line, or non-syndromic
hearing loss, drawn above it. A homozygous variant is drawn as a small loop at its position.

The protein position of a variant is taken from its p. description. If there is none, and the
variant is a single-nucleotide substitution within the coding sequence (c.1234A>G; no intronic
+/- offsets, no dup, del, ins), the position is approximated as the codon it falls in and
labelled "~p.<codon>".

Rows are skipped (with a note on stderr) when
  - an "other variants" column exists and is not empty,
  - an allele cell holds more than one variant (complex allele),
  - an allele cell is not empty but has neither a p. variant nor a coding substitution,
  - there is no second allele (allele 2 empty and allele 1 not marked homozygous),
  - the phenotype is in neither class,
  - the same pubmed id + patient id was already seen.

The figure goes to -o (default, for a single input: the input with its extension replaced by
.svg). The patients in the figure are listed on stdout.
"""

import argparse
import csv
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from html import escape
from pathlib import Path

#######################################
# equivalence classes

# column names, compared after lowercasing and collapsing whitespace; when a table has
# more than one column of a class, the one listed first is used
PUBMED_COLUMNS = ("pubmed", "pmid", "pubmed id", "pubmed_id")
PATIENT_COLUMNS = ("patient id", "patient_id", "patient", "patient #")
ALLELE1_COLUMNS = ("allele 1", "allele #1", "allele_1", "allele #1: hgvs variant")
ALLELE2_COLUMNS = ("allele 2", "allele #2", "allele_2", "allele #2: hgvs variant")
AGE_COLUMNS = ("age", "age at test", "age range", "age (range)")
PHENOTYPE_COLUMNS = ("phenotype short", "phenotype", "provided diagnosis", "diagnosis")
# optional
ZYGOSITY1_COLUMNS = ("allele 1 zygosity", "allele #1: zygosity", "zygosity 1")
OTHER_COLUMNS = ("other variants", "other variants: hgvs variant")
GENE_COLUMNS = ("gene",)

REQUIRED_COLUMNS = {
    "pubmed": PUBMED_COLUMNS,
    "patient": PATIENT_COLUMNS,
    "allele1": ALLELE1_COLUMNS,
    "allele2": ALLELE2_COLUMNS,
    "age": AGE_COLUMNS,
    "phenotype": PHENOTYPE_COLUMNS,
}
OPTIONAL_COLUMNS = {
    "zygosity1": ZYGOSITY1_COLUMNS,
    "other": OTHER_COLUMNS,
    "gene": GENE_COLUMNS,
}

# phenotypes: Usher syndrome 1B is matched exactly (case-insensitively, whitespace collapsed),
# non-syndromic hearing loss by any of the patterns occurring anywhere in the cell
USH1B_PHENOTYPES = {"ush1b", "usher syndrome 1b", "usher syndrome type 1b"}
NSHL_PHENOTYPE_PATTERNS = [re.compile(r"\bnon-?syndromic hearing loss\b", re.I), re.compile(r"\bnshl\b", re.I)]
ABOVE_LABEL = "non-syndromic hearing loss"
BELOW_LABEL = "Usher syndrome 1B"

HOMOZYGOUS = {"hom", "homozygous"}
EMPTY_CELL = {"", "_", "-", "n/a", "na", "none"}
DEFAULT_MAX_AGE = 10  # arcs of patients older than this (or of unknown age) are faded

# p.Asp1266Gly, p.Cys31Stop, p.(Arg212His), p.R212H, p.*2216Gln
PROTEIN_VARIANT = re.compile(r"\bp\.\(?(?:[A-Z][a-z]{2}|[A-Z*])(\d+)")
# every variant in a cell carries its own coding-level description
CODING_VARIANT = re.compile(r"(?<![A-Za-z])c\.")
# c.1234A>G only: no intronic offset (c.1234+1G>A, c.-5C>T, c.*10A>G), no dup/del/ins
CODING_SUBSTITUTION = re.compile(r"(?<![A-Za-z])c\.(\d+)[ACGT]>[ACGT](?![\w+\-])")
NUMBER = re.compile(r"\d+(?:\.\d+)?")

HEADER_SEARCH_ROWS = 10  # how far down a table the header row is looked for

TICK_STEP = 100

# layout, in SVG user units (px)
WIDTH = 1200
MARGIN_LEFT = 60
MARGIN_RIGHT = 40
PANEL_HEIGHT = 260  # room for the tallest arc on either side of the protein line
TITLE_HEIGHT = 50
LEGEND_HEIGHT = 30
AXIS_LABEL_SPACE = 40  # tick labels and axis title below the protein line
LOOP_HEIGHT = 22  # how far a homozygous loop reaches above/below the protein line
# A homozygous pair is drawn as a cubic Bezier loop that leaves and returns to the same
# point on the protein line. Its two control points sit this far left and right of that
# point; the loop itself bulges out only ~0.29 of this distance to each side.
LOOP_CONTROL_POINT_OFFSET = 28

# line thicknesses, in SVG user units (px)
ARC_STROKE_WIDTH = 5  # arcs between two different positions (and their legend swatches)
LOOP_STROKE_WIDTH = 3  # homozygous loops
PROTEIN_LINE_STROKE_WIDTH = 2
TICK_STROKE_WIDTH = 1
MARK_STROKE_WIDTH = 5  # the --mark arc / loop
# arc opacity; patients older than --max-age (or of unknown age) are drawn faded
ARC_OPACITY = 0.75
FADED_ARC_OPACITY = 0.2

# font sizes, in SVG user units (px)
TITLE_FONT_SIZE = 22
SOURCE_FONT_SIZE = 16  # "source: <file>" line under the title
LEGEND_FONT_SIZE = 12
TICK_LABEL_FONT_SIZE = 16
AXIS_TITLE_FONT_SIZE = 12

# palette: categorical slots 1 and 2 of the reference data-viz palette (light mode)
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
AXIS_COLOR = "#52514e"
ABOVE_COLOR = "#785ef0"
BELOW_COLOR = "#dc267f"
MARK_COLOR = "#ffb000"
FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"


class InputError(Exception):
    """Raised when an input does not satisfy the script's assumptions."""


def die(message: str, exit_code: int = 1) -> None:
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(exit_code)


def note(message: str) -> None:
    print(f"Note: {message}", file=sys.stderr)


@dataclass
class AllelePair:
    pubmed: str
    patient: str
    first: int
    second: int
    above: bool
    age_range: str  # decade of life, "" if unknown
    faded: bool  # age unknown or above --max-age
    label: str  # tooltip text


@dataclass
class Table:
    source: str  # file name, plus sheet name for .xlsx
    rows: list[dict[str, str]]  # keyed by the class names in REQUIRED_COLUMNS / OPTIONAL_COLUMNS


#######################################
# input
def check_input_path(raw_path: str) -> Path:
    path = Path(raw_path).expanduser()
    if not path.exists():
        raise InputError(f"input file '{path}' does not exist")
    if not path.is_file():
        raise InputError(f"input path '{path}' is not a regular file")
    if path.suffix.lower() not in (".tsv", ".xlsx"):
        raise InputError(f"input file '{path}' has neither a .tsv nor an .xlsx extension")
    if not os.access(path, os.R_OK):
        raise InputError(f"input file '{path}' is not readable")
    if path.stat().st_size == 0:
        raise InputError(f"input file '{path}' is empty")
    return path


def normalize_header(cell: str) -> str:
    return " ".join(cell.split()).lower()


def cell_text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return " ".join(str(value).split())


def column_indices(header: list[str]) -> dict[str, int] | None:
    """Column index per class, or None if a required class is missing from this row."""
    names = [normalize_header(cell) for cell in header]
    indices = {}
    for key, aliases in {**REQUIRED_COLUMNS, **OPTIONAL_COLUMNS}.items():
        for alias in aliases:
            if alias in names:
                indices[key] = names.index(alias)
                break
    if not all(key in indices for key in REQUIRED_COLUMNS):
        return None
    return indices


def make_table(source: str, raw_rows: list[list[str]]) -> Table:
    for header_number, header in enumerate(raw_rows[:HEADER_SEARCH_ROWS]):
        indices = column_indices(header)
        if indices is not None:
            break
    else:
        expected = "; ".join(f"{key}: {' / '.join(aliases)}" for key, aliases in REQUIRED_COLUMNS.items())
        raise InputError(f"'{source}': no header row with all required columns in the first {HEADER_SEARCH_ROWS} rows ({expected})")
    rows = []
    for raw in raw_rows[header_number + 1 :]:
        rows.append({key: raw[index] if index < len(raw) else "" for key, index in indices.items()})
    return Table(source, rows)


def read_tsv(path: Path) -> list[Table]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            raw_rows = [[cell_text(cell) for cell in row] for row in csv.reader(handle, delimiter="\t")]
    except UnicodeDecodeError as exc:
        raise InputError(f"'{path}' is not UTF-8 text: {exc}") from exc
    except csv.Error as exc:
        raise InputError(f"'{path}' is not a valid TSV file: {exc}") from exc
    return [make_table(path.name, raw_rows)]


def read_xlsx(path: Path) -> list[Table]:
    import openpyxl  # only needed for .xlsx input
    from zipfile import BadZipFile

    try:
        workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    except (BadZipFile, KeyError, ValueError) as exc:
        raise InputError(f"'{path}' is not a valid .xlsx file: {exc}") from exc
    tables = []
    try:
        for sheet in workbook.worksheets:
            source = f"{path.name} [{sheet.title}]" if len(workbook.worksheets) > 1 else path.name
            raw_rows = [[cell_text(cell) for cell in row] for row in sheet.iter_rows(values_only=True)]
            try:
                tables.append(make_table(source, raw_rows))
            except InputError as exc:
                if len(workbook.worksheets) == 1:
                    raise
                note(f"sheet skipped - {exc}")
    finally:
        workbook.close()
    if not tables:
        raise InputError(f"'{path}': no sheet has all required columns")
    return tables


def read_tables(path: Path) -> list[Table]:
    return read_xlsx(path) if path.suffix.lower() == ".xlsx" else read_tsv(path)


def is_empty(cell: str) -> bool:
    return cell.lower() in EMPTY_CELL


def allele_position(cell: str) -> tuple[int | None, str | None]:
    """Protein position of the single variant in an allele cell, or (None, reason to skip)."""
    if len(CODING_VARIANT.findall(cell)) > 1 or len(PROTEIN_VARIANT.findall(cell)) > 1:
        return None, f"more than one variant in allele '{cell}'"
    match = PROTEIN_VARIANT.search(cell)
    if match is not None:
        return int(match.group(1)), None
    match = CODING_SUBSTITUTION.search(cell)
    if match is not None:
        return codon(int(match.group(1))), None
    return None, f"neither a p. variant nor a coding substitution in allele '{cell}'"


def codon(coding_position: int) -> int:
    return (coding_position + 2) // 3


def protein_change(cell: str) -> str:
    match = re.search(r"\bp\.\S+", cell)
    if match is not None:
        return match.group(0).rstrip(",;")
    match = CODING_SUBSTITUTION.search(cell)
    if match is not None:
        return f"{match.group(0)} (~p.{codon(int(match.group(1)))})"
    return cell


def phenotype_side(phenotype: str) -> bool | None:
    """True for non-syndromic hearing loss (above), False for USH1B (below), else None."""
    if normalize_header(phenotype) in USH1B_PHENOTYPES:
        return False
    if any(pattern.search(phenotype) for pattern in NSHL_PHENOTYPE_PATTERNS):
        return True
    return None


def age_upper_bound(cell: str) -> float | None:
    """The age, or the upper end of an age range, in years; None if no number is given."""
    numbers = [float(n) for n in NUMBER.findall(cell)]
    return max(numbers) if numbers else None


def decade_of_life(age: float) -> str:
    """'0-10' for ages up to 10, then '11-20', '21-30', ..."""
    if age <= 10:
        return "0-10"
    low = (math.ceil(age) - 1) // 10 * 10 + 1
    return f"{low}-{low + 9}"


def allele_pairs(tables: list[Table], max_age: float) -> tuple[list[AllelePair], int]:
    """The usable allele pairs, and the number of rows skipped."""
    pairs = []
    seen = set()
    skipped = 0
    for table in tables:
        for row in table.rows:
            if all(is_empty(value) for value in row.values()):
                continue
            pubmed, patient = row["pubmed"] or "?", row["patient"] or "?"

            def skip(reason: str) -> None:
                nonlocal skipped
                skipped += 1
                note(f"{table.source}: PMID {pubmed}, patient {patient}: skipped - {reason}")

            if not is_empty(row.get("other", "")):
                skip("other variants present")
                continue

            allele1, allele2 = row["allele1"], row["allele2"]
            if row.get("zygosity1", "").lower() in HOMOZYGOUS:
                if not is_empty(allele2) and allele2 != allele1:
                    skip("allele 1 is homozygous, but a different allele 2 is also given")
                    continue
                allele2 = allele1
            if is_empty(allele1) or is_empty(allele2):
                skip("only one allele given" if not (is_empty(allele1) and is_empty(allele2)) else "no alleles given")
                continue

            first, reason = allele_position(allele1)
            if first is None:
                skip(reason)
                continue
            second, reason = allele_position(allele2)
            if second is None:
                skip(reason)
                continue

            above = phenotype_side(row["phenotype"])
            if above is None:
                skip(f"phenotype '{row['phenotype']}'")
                continue

            key = (pubmed, patient)
            if key in seen:
                skip("pubmed id and patient id already seen")
                continue
            seen.add(key)

            age = age_upper_bound(row["age"])
            age_range = decade_of_life(age) if age is not None else ""
            label = (
                f"PMID {pubmed}, patient {patient}: {protein_change(allele1)} / {protein_change(allele2)} "
                f"({row['phenotype']}, age {row['age'] or 'unknown'})"
            )
            faded = age is None or age > max_age
            pairs.append(AllelePair(pubmed, patient, first, second, above, age_range, faded, label))
    return pairs, skipped


#######################################
# drawing
def pair_path(first: int, second: int, above: bool, x_of, axis_y: float, height_per_px: float) -> str:
    """SVG path of an arc between two positions, or of a loop if they are the same."""
    direction = -1 if above else 1
    x1, x2 = sorted((x_of(first), x_of(second)))
    if first == second:
        tip = axis_y + direction * LOOP_HEIGHT * 4 / 3  # a cubic Bezier reaches 3/4 of its control height
        return (
            f"M {x1:.1f} {axis_y} C {x1 - LOOP_CONTROL_POINT_OFFSET:.1f} {tip:.1f} "
            f"{x1 + LOOP_CONTROL_POINT_OFFSET:.1f} {tip:.1f} {x1:.1f} {axis_y}"
        )
    rx = (x2 - x1) / 2
    ry = rx * height_per_px
    sweep = 1 if above else 0
    return f"M {x1:.1f} {axis_y} A {rx:.1f} {ry:.1f} 0 0 {sweep} {x2:.1f} {axis_y}"


def format_age(age: float) -> str:
    return f"{age:g}"


def svg_document(
    pairs: list[AllelePair], gene: str, source_name: str, marks: list[int], mark_above: bool, max_age: float
) -> str:
    max_position = max([max(p.first, p.second) for p in pairs] + marks)
    axis_end = math.ceil(max_position / TICK_STEP) * TICK_STEP
    plot_width = WIDTH - MARGIN_LEFT - MARGIN_RIGHT
    axis_y = TITLE_HEIGHT + LEGEND_HEIGHT + PANEL_HEIGHT
    height = axis_y + PANEL_HEIGHT + AXIS_LABEL_SPACE

    def x_of(position: int) -> float:
        return MARGIN_LEFT + (position - 1) / (axis_end - 1) * plot_width

    # an arc spanning the whole axis reaches the edge of its panel
    height_per_px = (PANEL_HEIGHT - 10) / (plot_width / 2)

    n_above = sum(p.above for p in pairs)
    n_below = len(pairs) - n_above
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" '
        f'viewBox="0 0 {WIDTH} {height}" font-family="{FONT}" role="img" '
        f'aria-label="{escape(gene)} allele pairs by phenotype">',
        f'<rect width="100%" height="100%" fill="{SURFACE}"/>',
        f'<text x="{MARGIN_LEFT}" y="28" font-size="{TITLE_FONT_SIZE}" font-weight="600" fill="{TEXT_PRIMARY}">'
        f"{escape(gene)}: protein positions of the two alleles per patient</text>",
        f'<text x="{MARGIN_LEFT}" y="46" font-size="{SOURCE_FONT_SIZE}" fill="{TEXT_SECONDARY}">'
        f"Source: {escape(source_name)}. Hover an arc for the patient and variants.</text>",
    ]

    # legend
    legend_y = TITLE_HEIGHT + 14
    legend_x = MARGIN_LEFT
    legend_entries = [
        (ABOVE_COLOR, f"Above: {ABOVE_LABEL} (n={n_above})"),
        (BELOW_COLOR, f"Below: {BELOW_LABEL} (n={n_below})"),
    ]
    if marks:
        legend_entries.append((MARK_COLOR, f"Query pair: residues {marks[0]} / {marks[1]}"))
    if any(p.faded for p in pairs):
        legend_entries.append((AXIS_COLOR, f"Faded: age above {format_age(max_age)} or unknown", FADED_ARC_OPACITY))
    for color, text, *opacity in legend_entries:
        swatch_length = 20
        opacity_attribute = f' stroke-opacity="{opacity[0]}"' if opacity else ""
        out.append(
            f'<line x1="{legend_x}" y1="{legend_y - 4}" x2="{legend_x + swatch_length}" y2="{legend_y - 4}" '
            f'stroke="{color}" stroke-width="{ARC_STROKE_WIDTH}" stroke-linecap="round"{opacity_attribute}/>'
        )
        text_x = legend_x + swatch_length + 6
        out.append(f'<text x="{text_x}" y="{legend_y}" font-size="{LEGEND_FONT_SIZE}" fill="{TEXT_PRIMARY}">{escape(text)}</text>')
        legend_x = text_x + 5.6 * len(text) + 30  # ~5.6 px per character at LEGEND_FONT_SIZE 12

    # arcs and loops
    out.append(f'<g fill="none" stroke-width="{ARC_STROKE_WIDTH}" stroke-linecap="round">')
    for pair in sorted(pairs, key=lambda p: -abs(p.second - p.first)):  # short arcs drawn on top
        color = ABOVE_COLOR if pair.above else BELOW_COLOR
        path = pair_path(pair.first, pair.second, pair.above, x_of, axis_y, height_per_px)
        width = f' stroke-width="{LOOP_STROKE_WIDTH}"' if pair.first == pair.second else ""
        opacity = FADED_ARC_OPACITY if pair.faded else ARC_OPACITY
        out.append(f'<path d="{path}" stroke="{color}"{width} stroke-opacity="{opacity}"><title>{escape(pair.label)}</title></path>')
    out.append("</g>")

    # the --mark pair: arc (or loop) on the requested side, drawn on top of the patient arcs
    if marks:
        label = f"Query pair: residues {marks[0]} / {marks[1]}"
        path = pair_path(marks[0], marks[1], mark_above, x_of, axis_y, height_per_px)
        out.append(
            f'<path d="{path}" fill="none" stroke="{MARK_COLOR}" stroke-width="{MARK_STROKE_WIDTH}" '
            f'stroke-linecap="round"><title>{escape(label)}</title></path>'
        )

    # variant positions
    positions = {(p.first, ABOVE_COLOR if p.above else BELOW_COLOR) for p in pairs}
    positions |= {(p.second, ABOVE_COLOR if p.above else BELOW_COLOR) for p in pairs}
    positions |= {(m, MARK_COLOR) for m in marks}
    for position, color in sorted(positions):
        out.append(
            f'<circle cx="{x_of(position):.1f}" cy="{axis_y}" r="3" fill="{color}" stroke="{SURFACE}" stroke-width="1.5">'
            f"<title>residue {position}</title></circle>"
        )

    # protein line and ticks
    out.append(
        f'<line x1="{MARGIN_LEFT}" y1="{axis_y}" x2="{MARGIN_LEFT + plot_width}" y2="{axis_y}" '
        f'stroke="{AXIS_COLOR}" stroke-width="{PROTEIN_LINE_STROKE_WIDTH}"/>'
    )
    for position in [1] + list(range(TICK_STEP, axis_end + 1, TICK_STEP)):
        x = x_of(position)
        major = position == 1 or position % 500 == 0
        out.append(
            f'<line x1="{x:.1f}" y1="{axis_y - 4}" x2="{x:.1f}" y2="{axis_y + (8 if major else 5)}" '
            f'stroke="{AXIS_COLOR}" stroke-width="{TICK_STROKE_WIDTH}"/>'
        )
        if major:
            out.append(
                f'<text x="{x:.1f}" y="{axis_y + 22}" font-size="{TICK_LABEL_FONT_SIZE}" text-anchor="middle" '
                f'fill="{TEXT_SECONDARY}">{position}</text>'
            )
    out.append(
        f'<text x="{MARGIN_LEFT + plot_width / 2:.1f}" y="{height - 6}" font-size="{AXIS_TITLE_FONT_SIZE}" text-anchor="middle" '
        f'fill="{TEXT_SECONDARY}">{escape(gene)} residue (tick every {TICK_STEP})</text>'
    )
    out.append("</svg>")
    return "\n".join(out) + "\n"


#######################################
# output
def write_text(text: str, out_path: Path) -> None:
    # write to a temp file in the same directory, then rename: no half-written output on failure
    fd, tmp_name = tempfile.mkstemp(dir=out_path.parent, prefix=f".{out_path.stem}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(tmp_name, out_path)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw allele pairs from any number of .tsv / .xlsx tables as arcs over the protein, into an .svg."
    )
    parser.add_argument("tables", nargs="+", help="input .tsv or .xlsx files")
    parser.add_argument(
        "-o", "--output",
        help="output .svg file (default for a single input: the input with the extension replaced by .svg; "
        "required for several inputs)",
    )
    parser.add_argument(
        "-a", "--max-age", type=float, default=DEFAULT_MAX_AGE,
        help=f"arcs of patients older than this, or of unknown age, are faded (default {DEFAULT_MAX_AGE})",
    )
    parser.add_argument(
        "-m", "--mark", nargs=3, metavar=("POS1", "POS2", "SIDE"), default=[],
        help="a pair of protein positions to draw as a highlighted arc (a loop if equal), on the auditory side "
        "(SIDE 'a', above) or the USH1B side (SIDE 'u', below); the axis is extended to include them",
    )
    parser.add_argument("-f", "--force", action="store_true", help="overwrite an existing output file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_age < 0:
        die(f"--max-age must not be negative, got {format_age(args.max_age)}")

    paths: list[Path] = []
    tables: list[Table] = []
    for raw_path in args.tables:
        try:
            path = check_input_path(raw_path)
            tables.extend(read_tables(path))
        except InputError as exc:
            die(str(exc))
        except OSError as exc:
            die(f"could not read '{raw_path}': {exc}")
        paths.append(path)

    marks: list[int] = []
    mark_above = False
    if args.mark:
        pos1, pos2, side = args.mark
        try:
            marks = [int(pos1), int(pos2)]
        except ValueError:
            die(f"--mark positions must be integers, got '{pos1}' '{pos2}'")
        if any(position < 1 for position in marks):
            die(f"--mark positions must be at least 1, got {pos1} {pos2}")
        if side not in ("a", "u"):
            die(f"--mark side must be 'a' (auditory, above) or 'u' (USH1B, below), got '{side}'")
        mark_above = side == "a"

    if args.output:
        out_path = Path(args.output).expanduser()
    elif len(paths) == 1:
        out_path = paths[0].with_suffix(".svg")
    else:
        die("with several input files, the output file must be given with -o")
    if out_path.suffix.lower() != ".svg":
        die(f"output file '{out_path}' does not have an .svg extension")
    if out_path.exists() and not args.force:
        die(f"output file '{out_path}' already exists (use --force to overwrite)")
    if out_path.exists() and not out_path.is_file():
        die(f"output path '{out_path}' exists and is not a regular file")
    if not out_path.parent.is_dir():
        die(f"output directory '{out_path.parent}' does not exist")

    pairs, skipped = allele_pairs(tables, args.max_age)
    if not pairs:
        die("no row has a usable pair of variants", exit_code=2)

    genes = sorted({row.get("gene", "") for table in tables for row in table.rows} - {""})
    gene = "/".join(genes) if genes else "protein"
    source_name = ", ".join(path.name for path in paths)
    try:
        write_text(svg_document(pairs, gene, source_name, marks, mark_above, args.max_age), out_path)
    except OSError as exc:
        die(f"could not write '{out_path}': {exc}")
    print(f"Wrote {len(pairs)} allele pairs ({skipped} rows skipped) to {out_path}", file=sys.stderr)

    # the patients that made it into the figure, on stdout
    print("PMID\tPatient ID\tAge (decade)")
    for pair in pairs:
        print(f"{pair.pubmed}\t{pair.patient}\t{pair.age_range or 'unknown'}")


if __name__ == "__main__":
    main()
