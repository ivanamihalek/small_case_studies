"""
SVG drawing shared by 03_allele_pair_visualization.py, 04_truncating_allele_count.py and
05_combined_visualization.py.

The protein is drawn by a ProteinTrack: by default (LINE_TRACK) a horizontal line, or a strip of domain
boxes (protein_domain_svg.py). Below it, a tick at every 100th residue, ending at the largest drawn position
(or the protein length, if the track knows it) rounded up to the next 100. The arcs and loops start at the
upper or lower edge of the track; on it, any combination of
  - patient allele pairs: arcs between the two variant positions, above the line for autosomal
    recessive non-syndromic hearing loss, below it for Usher syndrome 1B; a homozygous variant is a
    small loop; dotted when the patient's age range is other than SOLID_AGE,
  - patient stop / missense allele pairs: a lollipop at the missense position, on the same side
    as the arcs for the diagnosis, dotted by the same age rule,
  - gnomAD missense homozygotes: green loops above the line (the carriers are not known to be
    affected), sized by the number of homozygotes at the residue,
  - a query pair (--mark): an amber arc or loop on the requested side.
Only the standard library is used; the SVG is written as text.
"""

import math
import os
import tempfile
from dataclasses import dataclass, field
from html import escape
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

SOLID_AGE = "0-10"  # arcs of patients in any other age range are dotted
SCALES = ("linear", "log")

TICK_STEP = 100

# layout, in SVG user units (px)
WIDTH = 1200
MARGIN_LEFT = 60
MARGIN_RIGHT = 40
TITLE_HEIGHT = 50
LEGEND_HEIGHT = 30  # one legend row; every further row adds LEGEND_ROW_HEIGHT
LEGEND_ROW_HEIGHT = 18
LEGEND_CHAR_WIDTH = 5.6  # rough width of a character at LEGEND_FONT_SIZE, for laying out the legend
PANEL_HEIGHT = 260  # room for the tallest arc or loop on either side of the protein line
AXIS_LABEL_SPACE = 40  # tick labels and axis title below the protein line
# A loop leaves and returns to the same point on the protein line as a cubic Bezier. Its two
# control points sit left and right of that point; the loop itself bulges out only ~0.29 of
# this distance to each side.
PAIR_LOOP_HEIGHT = 22  # how far a homozygous patient loop reaches above/below the protein line
PAIR_LOOP_CONTROL_POINT_OFFSET = 28
SAME_PAIR_STEP = 6  # px added to the arc / loop height for each further patient with the same pair of positions
MAX_GNOMAD_LOOP_HEIGHT = PANEL_HEIGHT - 40  # the gnomAD residue with most homozygotes; room left for its label
GNOMAD_LOOP_WIDTH_RATIO = 0.8  # control point offset of a gnomAD loop, relative to its height
LABEL_FRACTION = 0.1  # gnomAD residues with at least this fraction of the largest count get a text label

# line thicknesses, in SVG user units (px)
ARC_STROKE_WIDTH = 5  # patient arcs between two different positions (and the legend swatches)
PAIR_LOOP_STROKE_WIDTH = 3  # homozygous patient loops
GNOMAD_LOOP_STROKE_WIDTH = 3
PROTEIN_LINE_STROKE_WIDTH = 2
TICK_STROKE_WIDTH = 1
MARK_STROKE_WIDTH = 5  # the --mark arc / loop
LOLLIPOP_STROKE_WIDTH = 3  # stop / missense patient lollipop stems
LOLLIPOP_HEIGHT = 50  # how far a lollipop head sits above/below the protein line
LOLLIPOP_HEAD_RADIUS = 6
# fine dots for patients older than SOLID_AGE: a zero-length dash with round line caps is drawn as a
# dot one stroke width across; the gap is measured between dot centres
FINE_DOTS = "0 7"

# font sizes, in SVG user units (px)
TITLE_FONT_SIZE = 22
SOURCE_FONT_SIZE = 16  # source line under the title
LEGEND_FONT_SIZE = 12
LOOP_LABEL_FONT_SIZE = 12
TICK_LABEL_FONT_SIZE = 16
AXIS_TITLE_FONT_SIZE = 12

# smaller gnomAD loops are drawn at the height of the axis tick labels, so that they stay visible (illustrative)
MIN_GNOMAD_LOOP_HEIGHT = TICK_LABEL_FONT_SIZE

# palette: categorical slots of the reference data-viz palette (light mode), plus a green for gnomAD
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
AXIS_COLOR = "#52514e"
ABOVE_COLOR = "#785ef0"
BELOW_COLOR = "#dc267f"
MARK_COLOR = "#ffb000"
GNOMAD_COLOR = "#1e9e55"
FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"


@dataclass
class AllelePair:
    patient: str
    first: int
    second: int
    above: bool
    age: str  # "Age (range)" as given, possibly empty
    dotted: bool  # age range other than SOLID_AGE
    label: str  # tooltip text


@dataclass
class Lollipop:
    patient: str
    position: int  # of the missense variant; the other allele is a stop codon
    above: bool
    age: str  # "Age (range)" as given, possibly empty
    dotted: bool  # age range other than SOLID_AGE
    label: str  # tooltip text


@dataclass
class MissenseVariant:
    variant_id: str
    protein_consequence: str
    allele_count: int
    homozygote_count: int


@dataclass
class ResidueLoop:
    position: int
    homozygotes: int
    variants: List[MissenseVariant]


@dataclass
class Figure:
    gene: str
    title: str
    subtitle: str
    pairs: List[AllelePair] = field(default_factory=list)  # patient arcs and loops
    lollipops: List[Lollipop] = field(default_factory=list)  # stop / missense patients
    loops: List[ResidueLoop] = field(default_factory=list)  # gnomAD homozygote loops
    scale: str = "linear"  # of the gnomAD loop size, one of SCALES
    marks: List[int] = field(default_factory=list)  # the query pair, if any
    mark_above: bool = False


# (x_of, axis_y, axis_end) -> SVG elements of the protein itself, centred on axis_y
TrackDrawer = Callable[[Callable[[int], float], float, int], List[str]]


@dataclass
class ProteinTrack:
    """How the protein itself is drawn; the arcs, loops and lollipops start at its upper or lower edge."""
    half_height: float  # px from the middle of the track to its upper / lower edge
    length: Optional[int]  # protein length in residues, if known; the axis reaches at least this far
    draw: TrackDrawer


def protein_line_svg(x_of: Callable[[int], float], axis_y: float, axis_end: int) -> List[str]:
    return [f'<line x1="{x_of(1):.2f}" y1="{axis_y}" x2="{x_of(axis_end):.2f}" y2="{axis_y}" '
            f'stroke="{AXIS_COLOR}" stroke-width="{PROTEIN_LINE_STROKE_WIDTH}"/>']


LINE_TRACK = ProteinTrack(0, None, protein_line_svg)


#######################################
# shapes
def pair_path(first: int, second: int, above: bool, x_of: Callable[[int], float], axis_y: float,
              height_per_px: float, extra_height: float = 0.0) -> str:
    """SVG path of an arc between two positions, or of a small loop if they are the same, extra_height px taller."""
    direction = -1 if above else 1
    x1, x2 = sorted((x_of(first), x_of(second)))
    if first == second:
        loop_height = PAIR_LOOP_HEIGHT + extra_height
        tip = axis_y + direction * loop_height * 4 / 3  # a cubic Bezier reaches 3/4 of its control height
        offset = PAIR_LOOP_CONTROL_POINT_OFFSET * loop_height / PAIR_LOOP_HEIGHT
        left, right = x1 - offset, x1 + offset
        return f"M {x1:.2f} {axis_y:.2f} C {left:.2f} {tip:.2f} {right:.2f} {tip:.2f} {x1:.2f} {axis_y:.2f}"
    rx = (x2 - x1) / 2
    ry = rx * height_per_px + extra_height
    sweep = 1 if above else 0
    return f"M {x1:.2f} {axis_y:.2f} A {rx:.2f} {ry:.2f} 0 0 {sweep} {x2:.2f} {axis_y:.2f}"


def gnomad_loop_height(homozygotes: int, largest: int, scale: str) -> float:
    if scale == "log":
        height = MAX_GNOMAD_LOOP_HEIGHT * math.log10(homozygotes + 1) / math.log10(largest + 1)
    else:
        height = MAX_GNOMAD_LOOP_HEIGHT * homozygotes / largest
    return max(height, MIN_GNOMAD_LOOP_HEIGHT)


def gnomad_loop_path(x: float, height: float, axis_y: float) -> str:
    """SVG path of a loop leaving and returning to (x, axis_y), reaching height above the line."""
    tip = axis_y - height * 4 / 3  # a cubic Bezier reaches 3/4 of its control height
    offset = height * GNOMAD_LOOP_WIDTH_RATIO
    return f"M {x:.2f} {axis_y:.2f} C {x - offset:.2f} {tip:.2f} {x + offset:.2f} {tip:.2f} {x:.2f} {axis_y:.2f}"


#######################################
# legend
@dataclass
class LegendEntry:
    color: str
    text: str
    dash: Optional[str] = None


def legend_entries(figure: Figure) -> List[LegendEntry]:
    entries: List[LegendEntry] = []
    patients = figure.pairs + figure.lollipops
    if patients:
        n_above = sum(p.above for p in patients)
        n_below = len(patients) - n_above
        entries.append(LegendEntry(ABOVE_COLOR, f"Above: autosomal recessive non-syndromic hearing loss (n={n_above})"))
        entries.append(LegendEntry(BELOW_COLOR, f"Below: Usher syndrome 1B (n={n_below})"))
    if figure.loops:
        largest = max(loop.homozygotes for loop in figure.loops)
        total = sum(loop.homozygotes for loop in figure.loops)
        measure = "homozygotes" if figure.scale == "linear" else "log10(homozygotes + 1)"
        text = (f"gnomAD missense homozygotes, loop size ~ {measure} (largest: {largest:,}; "
                f"{len(figure.loops)} residues, {total:,} homozygotes; small loops enlarged to stay visible)")
        entries.append(LegendEntry(GNOMAD_COLOR, text))
    if figure.lollipops:
        entries.append(LegendEntry(AXIS_COLOR, "Lollipop: stop codon / missense, at the missense residue"))
    if figure.marks:
        entries.append(LegendEntry(MARK_COLOR, f"Query pair: residues {figure.marks[0]} / {figure.marks[1]}"))
    if any(p.dotted for p in patients):
        entries.append(LegendEntry(AXIS_COLOR, f"Dotted: age other than {SOLID_AGE}", FINE_DOTS))
    return entries


def swatch_length(entry: LegendEntry) -> int:
    return 40 if entry.dash else 20


def legend_rows(entries: List[LegendEntry]) -> List[List[LegendEntry]]:
    """Lay the entries out left to right, starting a new row when the next one would not fit."""
    rows: List[List[LegendEntry]] = [[]]
    x = MARGIN_LEFT
    for entry in entries:
        entry_width = swatch_length(entry) + 6 + LEGEND_CHAR_WIDTH * len(entry.text)
        if rows[-1] and x + entry_width > WIDTH - MARGIN_RIGHT:
            rows.append([])
            x = MARGIN_LEFT
        rows[-1].append(entry)
        x += entry_width + 30
    return rows


def legend_svg(rows: List[List[LegendEntry]]) -> List[str]:
    out = []
    for row_index, row in enumerate(rows):
        y = TITLE_HEIGHT + 14 + row_index * LEGEND_ROW_HEIGHT
        x = MARGIN_LEFT
        for entry in row:
            length = swatch_length(entry)
            dash = f' stroke-dasharray="{entry.dash}"' if entry.dash else ""
            out.append(f'<line x1="{x}" y1="{y - 4}" x2="{x + length}" y2="{y - 4}" stroke="{entry.color}" '
                       f'stroke-width="{ARC_STROKE_WIDTH}" stroke-linecap="round"{dash}/>')
            text_x = x + length + 6
            out.append(f'<text x="{text_x}" y="{y}" font-size="{LEGEND_FONT_SIZE}" fill="{TEXT_PRIMARY}">'
                       f"{escape(entry.text)}</text>")
            x = text_x + LEGEND_CHAR_WIDTH * len(entry.text) + 30
    return out


#######################################
# figure
def svg_document(figure: Figure, track: ProteinTrack = LINE_TRACK) -> str:
    """The whole figure: patient arcs, gnomAD homozygote loops and the query pair, whichever are given."""
    positions = [p.first for p in figure.pairs] + [p.second for p in figure.pairs]
    positions += [p.position for p in figure.lollipops]
    positions += [loop.position for loop in figure.loops] + figure.marks
    positions += [track.length] if track.length else []
    axis_end = math.ceil(max(positions) / TICK_STEP) * TICK_STEP
    plot_width = WIDTH - MARGIN_LEFT - MARGIN_RIGHT
    rows = legend_rows(legend_entries(figure))
    legend_height = LEGEND_HEIGHT + (len(rows) - 1) * LEGEND_ROW_HEIGHT
    axis_y = TITLE_HEIGHT + legend_height + PANEL_HEIGHT + track.half_height
    top_y, bottom_y = axis_y - track.half_height, axis_y + track.half_height
    # the panel below the line holds the Usher 1B arcs (or a query pair placed there)
    has_lower_panel = bool(figure.pairs) or bool(figure.lollipops) or (bool(figure.marks) and not figure.mark_above)
    height = bottom_y + (PANEL_HEIGHT if has_lower_panel else 0) + AXIS_LABEL_SPACE

    def x_of(position: int) -> float:
        return MARGIN_LEFT + (position - 1) / (axis_end - 1) * plot_width

    def edge_y(above: bool) -> float:
        """Where the arcs, loops and lollipops of the given side leave the protein track."""
        return top_y if above else bottom_y

    # an arc spanning the whole axis reaches the edge of its panel
    height_per_px = (PANEL_HEIGHT - 10) / (plot_width / 2)

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" viewBox="0 0 {WIDTH} {height}" '
        f'font-family="{FONT}" role="img" aria-label="{escape(figure.title)}">',
        f'<rect width="100%" height="100%" fill="{SURFACE}"/>',
        f'<text x="{MARGIN_LEFT}" y="28" font-size="{TITLE_FONT_SIZE}" font-weight="600" fill="{TEXT_PRIMARY}">'
        f"{escape(figure.title)}</text>",
        f'<text x="{MARGIN_LEFT}" y="46" font-size="{SOURCE_FONT_SIZE}" fill="{TEXT_SECONDARY}">'
        f"{escape(figure.subtitle)}</text>",
    ]
    out.extend(legend_svg(rows))

    # patient arcs and loops, the short arcs drawn on top; each further patient with the same pair of positions
    # (on the same side) is drawn SAME_PAIR_STEP taller, so that the repeats are nested rather than hidden
    if figure.pairs:
        repeats: Dict[Tuple[int, int, bool], int] = {}
        nested: List[Tuple[AllelePair, int]] = []
        for pair in figure.pairs:
            key = (min(pair.first, pair.second), max(pair.first, pair.second), pair.above)
            nested.append((pair, repeats.get(key, 0)))
            repeats[key] = repeats.get(key, 0) + 1
        out.append(f'<g fill="none" stroke-width="{ARC_STROKE_WIDTH}" stroke-linecap="round" stroke-opacity="0.75">')
        for pair, repeat in sorted(nested, key=lambda item: (-abs(item[0].second - item[0].first), -item[1])):
            color = ABOVE_COLOR if pair.above else BELOW_COLOR
            extra_height = repeat * SAME_PAIR_STEP
            path = pair_path(pair.first, pair.second, pair.above, x_of, edge_y(pair.above), height_per_px, extra_height)
            width = f' stroke-width="{PAIR_LOOP_STROKE_WIDTH}"' if pair.first == pair.second else ""
            dash = f' stroke-dasharray="{FINE_DOTS}"' if pair.dotted else ""
            out.append(f'<path d="{path}" stroke="{color}"{width}{dash}><title>{escape(pair.label)}</title></path>')
        out.append("</g>")

    # stop / missense patients: a stem from the protein line to a round head at the missense position
    if figure.lollipops:
        out.append(f'<g stroke-width="{LOLLIPOP_STROKE_WIDTH}" stroke-linecap="round" stroke-opacity="0.75" '
                   f'fill-opacity="0.75">')
        for lollipop in figure.lollipops:
            color = ABOVE_COLOR if lollipop.above else BELOW_COLOR
            x = x_of(lollipop.position)
            stem_y = edge_y(lollipop.above)
            head_y = stem_y + (-1 if lollipop.above else 1) * LOLLIPOP_HEIGHT
            dash = f' stroke-dasharray="{FINE_DOTS}"' if lollipop.dotted else ""
            tooltip = f"<title>{escape(lollipop.label)}</title>"
            out.append(f'<line x1="{x:.2f}" y1="{stem_y}" x2="{x:.2f}" y2="{head_y:.2f}" stroke="{color}"{dash}>'
                       f"{tooltip}</line>")
            out.append(f'<circle cx="{x:.2f}" cy="{head_y:.2f}" r="{LOLLIPOP_HEAD_RADIUS}" fill="{color}">'
                       f"{tooltip}</circle>")
        out.append("</g>")

    # gnomAD homozygote loops, the largest drawn first so that the small ones stay visible on top
    if figure.loops:
        largest = max(loop.homozygotes for loop in figure.loops)
        out.append(f'<g fill="none" stroke="{GNOMAD_COLOR}" stroke-width="{GNOMAD_LOOP_STROKE_WIDTH}" '
                   f'stroke-linecap="round" stroke-opacity="0.75">')
        for loop in sorted(figure.loops, key=lambda loop: -loop.homozygotes):
            loop_height = gnomad_loop_height(loop.homozygotes, largest, figure.scale)
            path = gnomad_loop_path(x_of(loop.position), loop_height, top_y)
            variants = "; ".join(f"{v.protein_consequence} ({v.homozygote_count:,})" for v in loop.variants)
            tooltip = f"residue {loop.position}: {loop.homozygotes:,} gnomAD homozygotes - {variants}"
            out.append(f'<path d="{path}"><title>{escape(tooltip)}</title></path>')
        out.append("</g>")

        # labels on the loops of the residues with the most homozygotes, with a halo to read over the arcs
        for loop in figure.loops:
            if loop.homozygotes < LABEL_FRACTION * largest:
                continue
            label_y = top_y - gnomad_loop_height(loop.homozygotes, largest, figure.scale) - 8
            changes = "/".join(v.protein_consequence for v in loop.variants)
            label = f"{changes} ({loop.homozygotes:,})"
            out.append(f'<text x="{x_of(loop.position):.2f}" y="{label_y:.2f}" font-size="{LOOP_LABEL_FONT_SIZE}" '
                       f'text-anchor="middle" fill="{TEXT_PRIMARY}" stroke="{SURFACE}" stroke-width="3" '
                       f'paint-order="stroke">{escape(label)}</text>')

    # the query pair, drawn on top of the patient arcs
    if figure.marks:
        label = f"Query pair: residues {figure.marks[0]} / {figure.marks[1]}"
        path = pair_path(figure.marks[0], figure.marks[1], figure.mark_above, x_of, edge_y(figure.mark_above),
                         height_per_px)
        out.append(f'<path d="{path}" fill="none" stroke="{MARK_COLOR}" stroke-width="{MARK_STROKE_WIDTH}" '
                   f'stroke-linecap="round"><title>{escape(label)}</title></path>')

    # the protein itself, over the ends of the arcs and loops
    out.extend(track.draw(x_of, axis_y, axis_end))

    # variant positions, on the edge of the protein track where their arcs, loops or lollipops leave it
    dots = {(p.first, edge_y(p.above), ABOVE_COLOR if p.above else BELOW_COLOR) for p in figure.pairs}
    dots |= {(p.second, edge_y(p.above), ABOVE_COLOR if p.above else BELOW_COLOR) for p in figure.pairs}
    dots |= {(p.position, edge_y(p.above), ABOVE_COLOR if p.above else BELOW_COLOR) for p in figure.lollipops}
    dots |= {(loop.position, top_y, GNOMAD_COLOR) for loop in figure.loops}
    dots |= {(m, edge_y(figure.mark_above), MARK_COLOR) for m in figure.marks}
    for position, y, color in sorted(dots):
        out.append(f'<circle cx="{x_of(position):.2f}" cy="{y}" r="3" fill="{color}" stroke="{SURFACE}" '
                   f'stroke-width="1.5"><title>residue {position}</title></circle>')

    # ticks: across a bare protein line, below a protein track with some height
    tick_top = bottom_y if track.half_height else axis_y - 4
    for position in [1] + list(range(TICK_STEP, axis_end + 1, TICK_STEP)):
        x = x_of(position)
        major = position == 1 or position % 500 == 0
        out.append(f'<line x1="{x:.2f}" y1="{tick_top}" x2="{x:.2f}" y2="{bottom_y + (8 if major else 5)}" '
                   f'stroke="{AXIS_COLOR}" stroke-width="{TICK_STROKE_WIDTH}"/>')
        if major:
            out.append(f'<text x="{x:.2f}" y="{bottom_y + 22}" font-size="{TICK_LABEL_FONT_SIZE}" '
                       f'text-anchor="middle" fill="{TEXT_SECONDARY}">{position}</text>')
    axis_title = f"{figure.gene} residue (tick every {TICK_STEP})"
    out.append(f'<text x="{MARGIN_LEFT + plot_width / 2:.2f}" y="{height - 6}" font-size="{AXIS_TITLE_FONT_SIZE}" '
               f'text-anchor="middle" fill="{TEXT_SECONDARY}">{escape(axis_title)}</text>')
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
