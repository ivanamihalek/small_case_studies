#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = []
# ///
"""
Count truncating alleles and missense homozygotes in a gnomAD variant table, and draw the
missense homozygotes as loops over the protein sequence line.

Usage:
    uv run 04_truncating_allele_count.py <gnomad.csv> [--min-homozygotes N] [--max-homozygotes M]
        [--scale linear|log] [--gene NAME] [--force]

Example:
    uv run 04_truncating_allele_count.py data/gnomad_sample.csv

The input is a CSV exported from the gnomAD browser variant table for a gene/transcript.
Two tables are printed:
  1. Truncating variants: a variant is counted when its "VEP Annotation" column holds
     stop_gained or frameshift_variant (alone, or as one of several '&'- or ','-separated
     terms). The number of such variants and the sum of their "Allele Count" values are
     reported per consequence and in total, along with the largest "Allele Number" among
     them as a rough denominator.
  2. Missense homozygotes: every missense_variant with at least --min-homozygotes (default 3)
     and at most --max-homozygotes (default 1000) in its "Homozygote Count", sorted by that
     count, followed by the total over all missense variants. The upper bound leaves out the
     common polymorphisms, whose counts would dwarf everything else in the figure.

The variants of table 2 are also drawn into an SVG next to the input (.csv replaced by .svg),
by protein_svg.py, as in 03_allele_pair_visualization.py: the protein is a horizontal line with a
tick at every 100th residue, and each residue carrying such a variant gets a green loop. All loops
are above the line - the gnomAD homozygotes are not known to be affected, so there is no disease side.
The loop size is proportional to the number of homozygotes at the residue (summed over the
variants there); with --scale log it is proportional to log10(homozygotes + 1) instead.
"""

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from protein_svg import SCALES, Figure, MissenseVariant, ResidueLoop, svg_document, write_text

CONSEQUENCE_COLUMN = "VEP Annotation"
ALLELE_COUNT_COLUMN = "Allele Count"
ALLELE_NUMBER_COLUMN = "Allele Number"
HOMOZYGOTE_COUNT_COLUMN = "Homozygote Count"
VARIANT_ID_COLUMN = "gnomAD ID"
PROTEIN_CONSEQUENCE_COLUMN = "Protein Consequence"
TRUNCATING_CONSEQUENCES = ("stop_gained", "frameshift_variant")
MISSENSE = "missense_variant"
MIN_HOMOZYGOTES = 3
MAX_HOMOZYGOTES = 1000
DEFAULT_GENE = "MYO7A"

# p.Ser1666Cys, p.S1666C
PROTEIN_POSITION = re.compile(r"\bp\.\(?(?:[A-Z][a-z]{2}|[A-Z*])(\d+)")


@dataclass
class Counts:
    truncating_variants: Dict[str, int] = field(default_factory=lambda: dict.fromkeys(TRUNCATING_CONSEQUENCES, 0))
    truncating_alleles: Dict[str, int] = field(default_factory=lambda: dict.fromkeys(TRUNCATING_CONSEQUENCES, 0))
    max_allele_number: int = 0
    missense_variant_count: int = 0
    missense_allele_count: int = 0
    missense_homozygote_total: int = 0
    missense_listed: List[MissenseVariant] = field(default_factory=list)  # those with min to max homozygotes


def parse_args() -> argparse.Namespace:
    description = "Count truncating alleles and missense homozygotes in a gnomAD table, draw the latter as an SVG."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("table", type=Path, help="CSV exported from the gnomAD browser")
    min_help = f"list and draw missense variants with at least this many homozygotes (default {MIN_HOMOZYGOTES})"
    parser.add_argument("-n", "--min-homozygotes", type=int, default=MIN_HOMOZYGOTES, help=min_help)
    max_help = f"list and draw missense variants with at most this many homozygotes (default {MAX_HOMOZYGOTES})"
    parser.add_argument("-m", "--max-homozygotes", type=int, default=MAX_HOMOZYGOTES, help=max_help)
    scale_help = "loop size proportional to the homozygote count (linear, default) or to log10(count + 1)"
    parser.add_argument("-s", "--scale", choices=SCALES, default="linear", help=scale_help)
    gene_help = f"gene name for the figure title (default {DEFAULT_GENE})"
    parser.add_argument("-g", "--gene", default=DEFAULT_GENE, help=gene_help)
    parser.add_argument("-f", "--force", action="store_true", help="overwrite an existing output SVG")
    return parser.parse_args()


#######################################
# counting
def parse_int(row: Dict[str, str], column: str, location: str) -> Optional[int]:
    """Return the integer in row[column]; an empty cell counts as 0, an unparsable one is reported and gives None."""
    value = row.get(column, "").strip()
    if not value:
        return 0
    try:
        return int(value)
    except ValueError:
        warning = f"{location}: bad {column.lower()} {value!r}, row skipped"
        print(warning, file=sys.stderr)
        return None


def count_variants(table: Path, min_homozygotes: int, max_homozygotes: int) -> Counts:
    counts = Counts()
    with table.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {CONSEQUENCE_COLUMN, ALLELE_COUNT_COLUMN, HOMOZYGOTE_COUNT_COLUMN}
        missing = required - set(reader.fieldnames or [])
        if missing:
            error = f"{table}: missing column(s): {', '.join(sorted(missing))}"
            sys.exit(error)
        for line_number, row in enumerate(reader, start=2):
            location = f"{table}:{line_number}"
            terms = set(re.split(r"[&,]", row[CONSEQUENCE_COLUMN].strip()))
            truncating = [c for c in TRUNCATING_CONSEQUENCES if c in terms]
            if not truncating and MISSENSE not in terms:
                continue
            allele_count = parse_int(row, ALLELE_COUNT_COLUMN, location)
            if allele_count is None:
                continue
            if truncating:
                # a variant annotated with both terms is credited to the first (stop_gained)
                consequence = truncating[0]
                counts.truncating_variants[consequence] += 1
                counts.truncating_alleles[consequence] += allele_count
                allele_number = row.get(ALLELE_NUMBER_COLUMN, "")
                if allele_number.isdigit():
                    counts.max_allele_number = max(counts.max_allele_number, int(allele_number))
                continue
            homozygote_count = parse_int(row, HOMOZYGOTE_COUNT_COLUMN, location)
            if homozygote_count is None:
                continue
            counts.missense_variant_count += 1
            counts.missense_allele_count += allele_count
            counts.missense_homozygote_total += homozygote_count
            if min_homozygotes <= homozygote_count <= max_homozygotes:
                variant_id = row.get(VARIANT_ID_COLUMN, "")
                protein_consequence = row.get(PROTEIN_CONSEQUENCE_COLUMN, "")
                variant = MissenseVariant(variant_id, protein_consequence, allele_count, homozygote_count)
                counts.missense_listed.append(variant)
    counts.missense_listed.sort(key=lambda v: (-v.homozygote_count, v.variant_id))
    return counts


#######################################
# tables
def format_table(title: str, header: Sequence[str], rows: Sequence[Sequence[str]], footer: Sequence[str],
                 text_columns: int = 1) -> str:
    """Box-drawn table: the first text_columns columns left-aligned, the rest (numbers) right-aligned."""
    all_rows = [header, *rows, footer]
    widths = [max(len(r[i]) for r in all_rows) for i in range(len(header))]

    def rule(left: str, middle: str, right: str) -> str:
        return left + middle.join("─" * (w + 2) for w in widths) + right

    def cells(row: Sequence[str]) -> str:
        padded = [f" {c:<{w}} " if i < text_columns else f" {c:>{w}} " for i, (c, w) in enumerate(zip(row, widths))]
        return "│" + "│".join(padded) + "│"

    lines = [title, rule("┌", "┬", "┐"), cells(header), rule("├", "┼", "┤")]
    lines.extend(cells(r) for r in rows)
    if rows:
        lines.append(rule("├", "┼", "┤"))
    lines.extend([cells(footer), rule("└", "┴", "┘")])
    return "\n".join(lines)


def truncating_table(counts: Counts) -> str:
    header = ("consequence", "variants", "alleles")
    rows: List[Tuple[str, ...]] = []
    for c in TRUNCATING_CONSEQUENCES:
        rows.append((c, f"{counts.truncating_variants[c]:,}", f"{counts.truncating_alleles[c]:,}"))
    total_variants = sum(counts.truncating_variants.values())
    total_alleles = sum(counts.truncating_alleles.values())
    footer = ("total", f"{total_variants:,}", f"{total_alleles:,}")
    return format_table("Truncating variants", header, rows, footer)


def missense_homozygote_table(counts: Counts, min_homozygotes: int, max_homozygotes: int) -> str:
    header = ("gnomAD ID", "protein change", "alleles", "homozygotes")
    rows: List[Tuple[str, ...]] = []
    for v in counts.missense_listed:
        rows.append((v.variant_id, v.protein_consequence, f"{v.allele_count:,}", f"{v.homozygote_count:,}"))
    footer_label = f"total, all {counts.missense_variant_count:,} missense"
    footer = (footer_label, "", f"{counts.missense_allele_count:,}", f"{counts.missense_homozygote_total:,}")
    listed = len(counts.missense_listed)
    bounds = f"{min_homozygotes:,} to {max_homozygotes:,}"
    title = f"Missense variants with {bounds} homozygotes ({listed:,} of {counts.missense_variant_count:,})"
    return format_table(title, header, rows, footer, text_columns=2)


#######################################
# figure
def residue_loops(variants: List[MissenseVariant]) -> List[ResidueLoop]:
    """Group the variants by protein position; variants without a readable position are reported and dropped."""
    by_position: Dict[int, ResidueLoop] = {}
    for variant in variants:
        match = PROTEIN_POSITION.search(variant.protein_consequence)
        if match is None:
            warning = f"{variant.variant_id}: no protein position in {variant.protein_consequence!r}, not drawn"
            print(warning, file=sys.stderr)
            continue
        position = int(match.group(1))
        loop = by_position.setdefault(position, ResidueLoop(position, 0, []))
        loop.homozygotes += variant.homozygote_count
        loop.variants.append(variant)
    return sorted(by_position.values(), key=lambda loop: loop.position)


def main() -> None:
    args = parse_args()
    if not args.table.is_file():
        error = f"{args.table}: no such file"
        sys.exit(error)
    if args.min_homozygotes < 1:
        error = f"--min-homozygotes must be at least 1, got {args.min_homozygotes}"
        sys.exit(error)
    if args.max_homozygotes < args.min_homozygotes:
        error = f"--max-homozygotes ({args.max_homozygotes}) is below --min-homozygotes ({args.min_homozygotes})"
        sys.exit(error)
    out_path = args.table.with_suffix(".svg")
    if out_path.exists() and not args.force:
        error = f"{out_path}: output file already exists (use --force to overwrite)"
        sys.exit(error)

    print(f"parsing {args.table}\n")
    counts = count_variants(args.table, args.min_homozygotes, args.max_homozygotes)
    print(truncating_table(counts))
    if counts.max_allele_number:
        print(f"max allele number: {counts.max_allele_number:,}")
    print()
    print(missense_homozygote_table(counts, args.min_homozygotes, args.max_homozygotes))

    loops = residue_loops(counts.missense_listed)
    range_text = f"{args.min_homozygotes:,} to {args.max_homozygotes:,}"
    if not loops:
        print(f"\nno missense variant with {range_text} homozygotes: no figure written")
        return
    title = f"{args.gene}: missense variants with {range_text} homozygotes in gnomAD"
    subtitle = f"Source: {args.table.name}. Homozygotes not known to be affected. Hover a loop for the variants."
    figure = Figure(args.gene, title, subtitle, loops=loops, scale=args.scale)
    try:
        write_text(svg_document(figure), out_path)
    except OSError as exc:
        error = f"could not write {out_path}: {exc}"
        sys.exit(error)
    print(f"\nwrote {len(loops)} residue loops to {out_path}")


if __name__ == "__main__":
    main()
