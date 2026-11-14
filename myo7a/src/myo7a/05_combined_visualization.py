#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = ["openpyxl"]
# ///
"""
Draw the patient allele pairs of 03_allele_pair_visualization.py and the gnomAD missense
homozygotes of 04_truncating_allele_count.py on the same protein line.

Usage:
    uv run 05_combined_visualization.py <table.tsv|table.xlsx|table.csv> [<table> ...] <gnomad.csv> [--min-homozygotes N]
        [--max-homozygotes M] [--scale linear|log] [--mark POS1 POS2 SIDE] [--protein domains|line]
        [--output OUT.svg] [--force]

Example:
    uv run 05_combined_visualization.py data/sloan-heggen.suppl1.tsv data/moorfield.xlsx data/ambrosio_myo7a.csv \
        data/gnomAD_v4.1.2_ENST00000409709.csv

The patient pairs are arcs, above the line for autosomal recessive non-syndromic hearing loss and
below it for Usher syndrome 1B, with the same patient filters as 03_allele_pair_visualization.py;
the patients of all tables (read by patient_tables.py) are drawn together, each patient id prefixed
with the first three characters of its table's file name ("amb 1"). The gnomAD table is
always the last positional argument, so a patient table may also have the .csv extension. The
gnomAD missense variants with --min-homozygotes to --max-homozygotes homozygotes are green loops
above the line (the carriers are not known to be affected), sized by the number of homozygotes
at the residue as in 04_truncating_allele_count.py. The optional --mark pair is an amber arc or
loop, as in 03_allele_pair_visualization.py. The protein is a strip of its well scoring domains
(--protein domains, the default; protein_domain_svg.py) or a bare line (--protein line).

The output defaults to the directory of the first table, named after the prefixes of all the patient
tables, in the order given, with the extension .gnomad.svg (data/slo_moo_amb.gnomad.svg for the example).
"""

import argparse
import importlib
import sys
from pathlib import Path
from typing import List, Tuple

from patient_record import InputError
from patient_tables import output_path, read_patient_tables
from protein_domain_svg import DEFAULT_TRACK, PROTEIN_TRACKS
from protein_svg import SCALES, Figure, svg_document, write_text

# the two scripts have names that are not valid identifiers, hence importlib
allele_pairs_script = importlib.import_module("03_allele_pair_visualization")
gnomad_script = importlib.import_module("04_truncating_allele_count")

OUTPUT_SUFFIX = ".gnomad.svg"


def parse_args() -> argparse.Namespace:
    description = "Draw patient allele pairs (arcs) and gnomAD missense homozygotes (green loops) into one SVG."
    parser = argparse.ArgumentParser(description=description)
    tables_help = "patient tables: the .tsv output of 01_pdf_table_parser.py, the Moorfields .xlsx, the Ambrosio .csv"
    parser.add_argument("tables", nargs="+", help=tables_help)
    parser.add_argument("gnomad", type=Path, help="CSV exported from the gnomAD browser")
    min_help = f"draw missense variants with at least this many homozygotes (default {gnomad_script.MIN_HOMOZYGOTES})"
    parser.add_argument("-n", "--min-homozygotes", type=int, default=gnomad_script.MIN_HOMOZYGOTES, help=min_help)
    max_help = f"draw missense variants with at most this many homozygotes (default {gnomad_script.MAX_HOMOZYGOTES})"
    parser.add_argument("-m", "--max-homozygotes", type=int, default=gnomad_script.MAX_HOMOZYGOTES, help=max_help)
    scale_help = "gnomAD loop size proportional to the homozygote count (linear, default) or to log10(count + 1)"
    parser.add_argument("-s", "--scale", choices=SCALES, default="linear", help=scale_help)
    mark_help = ("a pair of protein positions to draw as an amber arc (a loop if equal), on the auditory side "
                 "(SIDE 'a', above) or the USH1B side (SIDE 'u', below)")
    parser.add_argument("--mark", nargs=3, metavar=("POS1", "POS2", "SIDE"), default=[], help=mark_help)
    protein_help = f"draw the protein as a strip of domain boxes or as a bare line (default {DEFAULT_TRACK})"
    parser.add_argument("-p", "--protein", choices=PROTEIN_TRACKS, default=DEFAULT_TRACK, help=protein_help)
    output_help = f"output SVG (default: the patient tables' name prefixes joined by '_', plus {OUTPUT_SUFFIX})"
    parser.add_argument("-o", "--output", type=Path, help=output_help)
    parser.add_argument("-f", "--force", action="store_true", help="overwrite an existing output file")
    return parser.parse_args()


def parse_marks(mark: List[str]) -> Tuple[List[int], bool]:
    if not mark:
        return [], False
    pos1, pos2, side = mark
    try:
        marks = [int(pos1), int(pos2)]
    except ValueError:
        error = f"--mark positions must be integers, got '{pos1}' '{pos2}'"
        sys.exit(error)
    if any(position < 1 for position in marks):
        error = f"--mark positions must be at least 1, got {pos1} {pos2}"
        sys.exit(error)
    if side not in ("a", "u"):
        error = f"--mark side must be 'a' (auditory, above) or 'u' (USH1B, below), got '{side}'"
        sys.exit(error)
    return marks, side == "a"


def main() -> None:
    args = parse_args()
    try:
        paths, records = read_patient_tables(args.tables, prefix_ids=True)
    except InputError as exc:
        sys.exit(str(exc))
    except OSError as exc:
        error = f"could not read the patient tables: {exc}"
        sys.exit(error)
    if not args.gnomad.is_file():
        error = f"{args.gnomad}: no such file"
        sys.exit(error)
    if args.min_homozygotes < 1:
        error = f"--min-homozygotes must be at least 1, got {args.min_homozygotes}"
        sys.exit(error)
    if args.max_homozygotes < args.min_homozygotes:
        error = f"--max-homozygotes ({args.max_homozygotes}) is below --min-homozygotes ({args.min_homozygotes})"
        sys.exit(error)
    marks, mark_above = parse_marks(args.mark)

    out_path = args.output or output_path(paths, OUTPUT_SUFFIX)
    if out_path.exists() and not args.force:
        error = f"{out_path}: output file already exists (use --force to overwrite)"
        sys.exit(error)

    pairs, lollipops = allele_pairs_script.allele_pairs(records)
    counts = gnomad_script.count_variants(args.gnomad, args.min_homozygotes, args.max_homozygotes)
    loops = gnomad_script.residue_loops(counts.missense_listed)
    if not pairs and not lollipops and not loops:
        error = "neither a usable patient allele pair nor a gnomAD missense homozygote: nothing to draw"
        sys.exit(error)

    genes = sorted({record.gene for record in records} - {""})
    gene = "/".join(genes) if genes else "protein"
    range_text = f"{args.min_homozygotes:,} to {args.max_homozygotes:,}"
    title = f"{gene}: patient allele pairs and gnomAD missense homozygotes ({range_text})"
    sources = ", ".join(path.name for path in paths + [args.gnomad])
    subtitle = f"Sources: {sources}. Hover an arc or loop for the details."
    figure = Figure(gene, title, subtitle, pairs=pairs, lollipops=lollipops, loops=loops, scale=args.scale,
                    marks=marks, mark_above=mark_above)
    try:
        write_text(svg_document(figure, PROTEIN_TRACKS[args.protein]), out_path)
    except OSError as exc:
        error = f"could not write {out_path}: {exc}"
        sys.exit(error)
    print(f"wrote {len(pairs)} allele pairs, {len(lollipops)} stop / missense lollipops and {len(loops)} gnomAD "
          f"residue loops to {out_path}")


if __name__ == "__main__":
    main()
