"""
The MYO7A protein as a strip of domain boxes, an alternative to the bare protein line of protein_svg.py.

A thin backbone runs from residue 1 to MYO7A_LENGTH; on it, every domain scoring E <= MAX_E_VALUE is a
black and white box with rounded corners and a vertical gradient (for a slightly 3D look), labelled with
the domain name. The box spans the domain's residues on the same axis as the arcs and loops, so the
backbone between the boxes shows where the linkers lie in the sequence. A label too long for its box is
written just below the strip. The domains are SMART annotations; the unscored ones (coiled coil, low
complexity) are not drawn.

PROTEIN_TRACKS holds both protein representations, keyed by the name used on the command line of
03_allele_pair_visualization.py and 05_combined_visualization.py.
"""

from dataclasses import dataclass
from html import escape
from typing import Callable, Dict, List, Optional

from protein_svg import LINE_TRACK, SURFACE, TEXT_PRIMARY, ProteinTrack

MYO7A_LENGTH = 2215  # residues, UniProt Q13402
MAX_E_VALUE = 1e-3  # domains scoring worse than this (or not scored at all) are not drawn

# layout, in SVG user units (px)
DOMAIN_HALF_HEIGHT = 13  # half the height of a domain box
DOMAIN_CORNER_RADIUS = 6
DOMAIN_STROKE_WIDTH = 1.2
BACKBONE_STROKE_WIDTH = 3
DOMAIN_FONT_SIZE = 13
MIN_DOMAIN_FONT_SIZE = 9  # below this, the label goes under the strip instead
CHAR_WIDTH_RATIO = 0.62  # rough width of a (bold) character relative to the font size
LABEL_PADDING = 4  # px left free on either side of a label inside its box

# black and white, shaded from a light top to a darker bottom edge
BACKBONE_COLOR = "#0b0b0b"
GRADIENT_ID = "domain-shade"
GRADIENT_STOPS = [(0, "#ffffff"), (45, "#e8e8e8"), (100, "#9c9c9c")]  # (offset in %, color)


@dataclass(frozen=True)
class ProteinDomain:
    name: str
    start: int  # 1-based, inclusive
    end: int  # 1-based, inclusive
    e_value: Optional[float]  # None where the source table reports N/A


MYO7A_DOMAINS: List[ProteinDomain] = [
    ProteinDomain("MYSc", 59, 742, 0.0),
    ProteinDomain("IQ", 743, 765, 2.99e+00),
    ProteinDomain("IQ", 766, 788, 4.56e-05),
    ProteinDomain("IQ", 812, 834, 3.43e+00),
    ProteinDomain("IQ", 835, 857, 1.03e+01),
    ProteinDomain("coiled coil", 866, 935, None),
    ProteinDomain("low complexity", 979, 996, None),
    ProteinDomain("MyTH4", 1017, 1253, 1.73e-71),
    ProteinDomain("B41", 1254, 1469, 1.61e-44),
    ProteinDomain("SH3", 1606, 1671, 1.84e-06),
    ProteinDomain("MyTH4", 1747, 1896, 2.06e-60),
    ProteinDomain("B41", 1898, 2115, 1.57e-5),
]


def well_scoring(domains: List[ProteinDomain], max_e_value: float = MAX_E_VALUE) -> List[ProteinDomain]:
    return [d for d in domains if d.e_value is not None and d.e_value <= max_e_value]


def gradient_svg() -> List[str]:
    stops = [f'<stop offset="{offset}%" stop-color="{color}"/>' for offset, color in GRADIENT_STOPS]
    return [f'<defs><linearGradient id="{GRADIENT_ID}" x1="0" y1="0" x2="0" y2="1">', *stops,
            "</linearGradient></defs>"]


def domain_track_svg(x_of: Callable[[int], float], axis_y: float, axis_end: int) -> List[str]:
    """The backbone from residue 1 to MYO7A_LENGTH, with the well scoring domains as boxes on it."""
    top_y, bottom_y = axis_y - DOMAIN_HALF_HEIGHT, axis_y + DOMAIN_HALF_HEIGHT
    out = gradient_svg()
    out.append(f'<line x1="{x_of(1):.2f}" y1="{axis_y}" x2="{x_of(MYO7A_LENGTH):.2f}" y2="{axis_y}" '
               f'stroke="{BACKBONE_COLOR}" stroke-width="{BACKBONE_STROKE_WIDTH}" stroke-linecap="round"/>')
    for domain in well_scoring(MYO7A_DOMAINS):
        x1, x2 = x_of(domain.start), x_of(domain.end)
        width = x2 - x1
        radius = min(DOMAIN_CORNER_RADIUS, width / 2)
        description = f"{domain.name}: residues {domain.start}-{domain.end}, E = {domain.e_value:.2e}"
        tooltip = f"<title>{escape(description)}</title>"
        out.append(f'<rect x="{x1:.2f}" y="{top_y}" width="{width:.2f}" height="{2 * DOMAIN_HALF_HEIGHT}" '
                   f'rx="{radius:.2f}" fill="url(#{GRADIENT_ID})" stroke="{BACKBONE_COLOR}" '
                   f'stroke-width="{DOMAIN_STROKE_WIDTH}">{tooltip}</rect>')

        # the largest font (up to DOMAIN_FONT_SIZE) that fits the box, else a label under the strip
        font_size = min(DOMAIN_FONT_SIZE, (width - 2 * LABEL_PADDING) / (CHAR_WIDTH_RATIO * len(domain.name)))
        if font_size >= MIN_DOMAIN_FONT_SIZE:
            label_y = axis_y + 0.35 * font_size  # roughly centres the capitals vertically
            halo = ""
        else:
            font_size = DOMAIN_FONT_SIZE
            label_y = bottom_y + font_size + 2
            halo = f' stroke="{SURFACE}" stroke-width="3" paint-order="stroke"'
        out.append(f'<text x="{(x1 + x2) / 2:.2f}" y="{label_y:.2f}" font-size="{font_size:.2f}" font-weight="600" '
                   f'text-anchor="middle" fill="{TEXT_PRIMARY}"{halo}>{escape(domain.name)}{tooltip}</text>')
    return out


DOMAIN_TRACK = ProteinTrack(DOMAIN_HALF_HEIGHT, MYO7A_LENGTH, domain_track_svg)
PROTEIN_TRACKS: Dict[str, ProteinTrack] = {"domains": DOMAIN_TRACK, "line": LINE_TRACK}
DEFAULT_TRACK = "domains"
