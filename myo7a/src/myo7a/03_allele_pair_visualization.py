#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = []
# ///
"""
Draw the protein-level allele pairs of compound heterozygous / homozygous patients
as arcs over a protein sequence line.

Usage:
    uv run "03_allele_pair visualization.py" <table.tsv> [--mark POS1 POS2 SIDE] [--force]

Example:
    uv run "03_allele_pair visualization.py" data/sloan-heggen.suppl1.tsv --mark 100 1500 u

The input is the tab-separated output of 01_pdf_table_parser.py. The output is an SVG
next to the input, with the .tsv extension replaced by .svg.

The protein is drawn as a horizontal line with a tick at every 100th residue, ending at
the largest variant position rounded up to the next 100. Each patient's pair of variant
positions is connected by an arc: above the line for autosomal recessive non-syndromic
hearing loss, below it for Usher syndrome 1B. A homozygous variant is drawn as a small
loop at its position. The optional --mark pair is drawn the same way, in green, on the side
given by SIDE: 'a' above (auditory, non-syndromic hearing loss) or 'u' below (Usher 1B);
it is a loop if the two positions are equal, and the axis is extended if it lies beyond the
largest variant position. A patient's arc or loop is drawn with long dashes when the
"Age (range)" column is anything other than "0-10" (including empty).

Rows are skipped (with a note on stderr) when
  - the "Other variants: HGVS variant" column is not empty,
  - an allele cell holds more than one variant (complex allele),
  - an allele cell is not empty but has no protein-level (p.) variant,
  - there is no second allele (Allele #2 empty and Allele #1 not homozygous),
  - the diagnosis is neither of the two above.
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

ALLELE1_COLUMN = "Allele #1: HGVS variant"
ALLELE2_COLUMN = "Allele #2: HGVS variant"
ZYGOSITY1_COLUMN = "Allele #1: zygosity"
OTHER_COLUMN = "Other variants: HGVS variant"
DIAGNOSIS_COLUMN = "Provided Diagnosis"
PATIENT_COLUMN = "Patient #"
GENE_COLUMN = "Gene"
AGE_COLUMN = "Age (range)"
REQUIRED_COLUMNS = [ALLELE1_COLUMN, ALLELE2_COLUMN, ZYGOSITY1_COLUMN, OTHER_COLUMN, DIAGNOSIS_COLUMN, AGE_COLUMN]

ABOVE_DIAGNOSIS = "autosomal recessive non-syndromic hearing loss"
BELOW_DIAGNOSIS = "usher syndrome 1b"
SOLID_AGE = "0-10"  # arcs of patients in any other age range are dashed

# p.Asp1266Gly, p.Cys31Stop, p.(Arg212His), p.R212H, p.*2216Gln
PROTEIN_VARIANT = re.compile(r"\bp\.\(?(?:[A-Z][a-z]{2}|[A-Z*])(\d+)")
# every variant in a cell carries its own coding-level description
CODING_VARIANT = re.compile(r"(?<![A-Za-z])c\.")
EMPTY_CELL = {"", "_", "-"}

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
MARK_STROKE_WIDTH = 5  # the green --mark arc / loop
# long dash for patients older than SOLID_AGE: "dash gap" lengths; the round line caps
# eat into the gap by one stroke width, so the visible gap is shorter than given here
LONG_DASH = "18 14"

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
    patient: str
    first: int
    second: int
    above: bool
    age: str  # "Age (range)" as given, possibly empty
    dashed: bool  # age range other than SOLID_AGE
    label: str  # tooltip text


#######################################
# input
def check_tsv_path(raw_path: str) -> Path:
    path = Path(raw_path).expanduser()
    if not path.exists():
        raise InputError(f"input file '{path}' does not exist")
    if not path.is_file():
        raise InputError(f"input path '{path}' is not a regular file")
    if path.suffix.lower() != ".tsv":
        raise InputError(f"input file '{path}' does not have a .tsv extension")
    if not os.access(path, os.R_OK):
        raise InputError(f"input file '{path}' is not readable")
    if path.stat().st_size == 0:
        raise InputError(f"input file '{path}' is empty")
    return path


def read_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
            if missing:
                raise InputError(f"'{path}' lacks the column(s): {', '.join(missing)}")
            return [{key: (value or "").strip() for key, value in row.items() if key} for row in reader]
    except UnicodeDecodeError as exc:
        raise InputError(f"'{path}' is not UTF-8 text: {exc}") from exc
    except csv.Error as exc:
        raise InputError(f"'{path}' is not a valid TSV file: {exc}") from exc


def is_empty(cell: str) -> bool:
    return cell in EMPTY_CELL


def allele_position(cell: str) -> tuple[int | None, str | None]:
    """Protein position of the single variant in an allele cell, or (None, reason to skip)."""
    if len(CODING_VARIANT.findall(cell)) > 1 or len(PROTEIN_VARIANT.findall(cell)) > 1:
        return None, f"more than one variant in allele '{cell}'"
    match = PROTEIN_VARIANT.search(cell)
    if match is None:
        return None, f"no protein-level variant in allele '{cell}'"
    return int(match.group(1)), None


def protein_change(cell: str) -> str:
    match = re.search(r"\bp\.\S+", cell)
    return match.group(0).rstrip(",;") if match else cell


def allele_pairs(rows: list[dict[str, str]]) -> list[AllelePair]:
    pairs = []
    for row in rows:
        patient = row.get(PATIENT_COLUMN, "?")

        def skip(reason: str) -> None:
            note(f"patient {patient}: skipped - {reason}")

        if not is_empty(row[OTHER_COLUMN]):
            skip("other variants present")
            continue

        allele1, allele2 = row[ALLELE1_COLUMN], row[ALLELE2_COLUMN]
        if row[ZYGOSITY1_COLUMN].lower() == "hom":
            if not is_empty(allele2):
                skip("allele #1 is homozygous, but allele #2 is also given")
                continue
            allele2 = allele1
        if is_empty(allele1) or is_empty(allele2):
            skip("only one allele given")
            continue

        first, reason = allele_position(allele1)
        if first is None:
            skip(reason)
            continue
        second, reason = allele_position(allele2)
        if second is None:
            skip(reason)
            continue

        diagnosis = row[DIAGNOSIS_COLUMN].lower()
        if diagnosis not in (ABOVE_DIAGNOSIS, BELOW_DIAGNOSIS):
            skip(f"diagnosis '{row[DIAGNOSIS_COLUMN]}'")
            continue

        age = row[AGE_COLUMN]
        label = (
            f"Patient {patient}: {protein_change(allele1)} / {protein_change(allele2)} "
            f"({row[DIAGNOSIS_COLUMN]}, age {age or 'unknown'})"
        )
        pairs.append(AllelePair(patient, first, second, diagnosis == ABOVE_DIAGNOSIS, age, age != SOLID_AGE, label))
    return pairs


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


def svg_document(pairs: list[AllelePair], gene: str, source_name: str, marks: list[int], mark_above: bool) -> str:
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
        f'aria-label="{escape(gene)} allele pairs by diagnosis">',
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
        (ABOVE_COLOR, f"Above: autosomal recessive non-syndromic hearing loss (n={n_above})"),
        (BELOW_COLOR, f"Below: Usher syndrome 1B (n={n_below})"),
    ]
    if marks:
        legend_entries.append((MARK_COLOR, f"Query pair: residues {marks[0]} / {marks[1]}"))
    if any(p.dashed for p in pairs):
        legend_entries.append((AXIS_COLOR, f"Dashed: age other than {SOLID_AGE}", LONG_DASH))
    for color, text, *dash in legend_entries:
        swatch_length = 40 if dash else 20
        dash_attribute = f' stroke-dasharray="{dash[0]}"' if dash else ""
        out.append(
            f'<line x1="{legend_x}" y1="{legend_y - 4}" x2="{legend_x + swatch_length}" y2="{legend_y - 4}" '
            f'stroke="{color}" stroke-width="{ARC_STROKE_WIDTH}" stroke-linecap="round"{dash_attribute}/>'
        )
        text_x = legend_x + swatch_length + 6
        out.append(f'<text x="{text_x}" y="{legend_y}" font-size="{LEGEND_FONT_SIZE}" fill="{TEXT_PRIMARY}">{escape(text)}</text>')
        legend_x = text_x + 5.6 * len(text) + 30  # ~5.6 px per character at LEGEND_FONT_SIZE 12

    # arcs and loops
    out.append(f'<g fill="none" stroke-width="{ARC_STROKE_WIDTH}" stroke-linecap="round" stroke-opacity="0.75">')
    for pair in sorted(pairs, key=lambda p: -abs(p.second - p.first)):  # short arcs drawn on top
        color = ABOVE_COLOR if pair.above else BELOW_COLOR
        path = pair_path(pair.first, pair.second, pair.above, x_of, axis_y, height_per_px)
        width = f' stroke-width="{LOOP_STROKE_WIDTH}"' if pair.first == pair.second else ""
        dash = f' stroke-dasharray="{LONG_DASH}"' if pair.dashed else ""
        out.append(f'<path d="{path}" stroke="{color}"{width}{dash}><title>{escape(pair.label)}</title></path>')
    out.append("</g>")

    # the --mark pair: green arc (or loop) on the requested side, drawn on top of the patient arcs
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
        description="Draw allele pairs from a 01_pdf_table_parser.py .tsv file as arcs over the protein, into an .svg."
    )
    parser.add_argument("tsv", help="path to the input .tsv file")
    parser.add_argument(
        "-m", "--mark", nargs=3, metavar=("POS1", "POS2", "SIDE"), default=[],
        help="a pair of protein positions to draw as a green arc (a loop if equal), on the auditory side "
        "(SIDE 'a', above) or the USH1B side (SIDE 'u', below); the axis is extended to include them",
    )
    parser.add_argument("-f", "--force", action="store_true", help="overwrite an existing output file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        tsv_path = check_tsv_path(args.tsv)
        rows = read_rows(tsv_path)
    except InputError as exc:
        die(str(exc))
    except OSError as exc:
        die(f"could not read '{args.tsv}': {exc}")

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

    out_path = tsv_path.with_suffix(".svg")
    if out_path.exists() and not args.force:
        die(f"output file '{out_path}' already exists (use --force to overwrite)")
    if out_path.exists() and not out_path.is_file():
        die(f"output path '{out_path}' exists and is not a regular file")

    pairs = allele_pairs(rows)
    if not pairs:
        die("no row has a usable pair of protein-level variants", exit_code=2)

    genes = sorted({row.get(GENE_COLUMN, "") for row in rows} - {""})
    gene = "/".join(genes) if genes else "protein"
    try:
        write_text(svg_document(pairs, gene, tsv_path.name, marks, mark_above), out_path)
    except OSError as exc:
        die(f"could not write '{out_path}': {exc}")
    print(f"Wrote {len(pairs)} allele pairs ({len(rows) - len(pairs)} rows skipped) to {out_path}", file=sys.stderr)

    # the patients that made it into the figure, on stdout
    print(f"{PATIENT_COLUMN}\t{AGE_COLUMN}")
    for pair in pairs:
        print(f"{pair.patient}\t{pair.age or 'unknown'}")


if __name__ == "__main__":
    main()
