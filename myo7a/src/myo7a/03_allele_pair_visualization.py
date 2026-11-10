#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = ["openpyxl"]
# ///
"""
Draw the protein-level allele pairs of compound heterozygous / homozygous patients
as arcs over a protein sequence line.

Usage:
    uv run "03_allele_pair visualization.py" <table.tsv|table.xlsx|table.csv> [<table> ...] [--mark POS1 POS2 SIDE]
        [--protein domains|line] [--force]

Example:
    uv run "03_allele_pair visualization.py" data/sloan-heggen.suppl1.tsv data/moorfield.xlsx data/ambrosio_myo7a.csv \
        --mark 100 1500 u

The inputs are patient tables, read by patient_tables.py into input-format independent patient
records (patient_record.py) and drawn together: the tab-separated output of 01_pdf_table_parser.py
(sloan_heggen_tsv.py), the Moorfields .xlsx table (moorfields_xlsx.py) and the tab-separated
Ambrosio et al. .csv table (ambrosio_csv.py). Each patient id is prefixed with the first three
characters of its table's file name ("amb 1", "moo ...", "slo ..."). The output is an SVG in the
directory of the first input, named after the same prefixes of all inputs, in the order given
(data/slo_moo_amb.svg for the example); it is drawn by protein_svg.py.

The protein is drawn as a strip of its well scoring domains (--protein domains, the default;
protein_domain_svg.py) or as a bare horizontal line (--protein line), with a tick at every 100th
residue, ending at the largest variant position (or, for the domain strip, the protein length)
rounded up to the next 100. Each patient's pair of variant
positions is connected by an arc: above the line for autosomal recessive non-syndromic
hearing loss, below it for Usher syndrome 1B. A homozygous variant is drawn as a small
loop at its position. The optional --mark pair is drawn the same way, in amber, on the side
given by SIDE: 'a' above (auditory, non-syndromic hearing loss) or 'u' below (Usher 1B);
it is a loop if the two positions are equal, and the axis is extended if it lies beyond the
largest variant position. A patient with one truncating allele (a stop codon or a
frameshift) and one missense allele (a three-letter amino acid substitution, e.g. p.Asp1266Gly)
is drawn as a lollipop at the missense position, pointing up for non-syndromic hearing loss
and down for Usher syndrome 1B. A patient's arc, loop or lollipop is drawn with long dashes when the age range is anything
other than "0-10" (including unknown).

Patients are skipped (with a note on stderr) when
  - other variants are given for the patient,
  - an allele holds more than one variant (complex allele),
  - an allele is given but has no protein-level (p.) variant,
  - both alleles are truncating (stop codon: Ter, Stop or *; frameshift: fs), or one is
    truncating and the other is not a missense variant,
  - there is no second allele (allele #2 not given and allele #1 not homozygous),
    or allele #1 is homozygous and allele #2 is given as well,
  - the diagnosis is neither of the two above.
"""

import argparse
import re
import sys
from typing import List, Optional, Tuple

from patient_record import InputError, PatientRecord, Phenotype
from patient_tables import output_path, read_patient_tables
from protein_domain_svg import DEFAULT_TRACK, PROTEIN_TRACKS
from protein_svg import SOLID_AGE, AllelePair, Figure, Lollipop, svg_document, write_text

# p.Asp1266Gly, p.Cys31Stop, p.(Arg212His), p.R212H, p.*2216Gln
PROTEIN_VARIANT = re.compile(r"\bp\.\(?(?:[A-Z][a-z]{2}|[A-Z*])(\d+)")
# three-letter amino acid variant: reference, position, substitution (p.Asp1266Gly, p.Cys31Stop, p.(Arg212Ter))
CODING_VARIANT = re.compile(r"\bp\.\(?([A-Z][a-z]{2})(\d+)(Stop|[A-Z][a-z]{2}|\*)")
STOP_SUBSTITUTIONS = {"Ter", "Stop", "*"}
# p.(Asp2191ValfsTer4), p.Asp2191fs
FRAMESHIFT = re.compile(r"\bp\.\(?(?:[A-Z][a-z]{2}|[A-Z])\d+(?:[A-Z][a-z]{2}|[A-Z])?fs")


def die(message: str, exit_code: int = 1) -> None:
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(exit_code)


def note(message: str) -> None:
    print(f"Note: {message}", file=sys.stderr)


#######################################
# variants
def allele_position(cell: str) -> tuple[int | None, str | None]:
    """Protein position of the single variant in an allele cell, or (None, reason to skip)."""
    if len(CODING_VARIANT.findall(cell)) > 1 or len(PROTEIN_VARIANT.findall(cell)) > 1:
        return None, f"more than one variant in allele '{cell}'"
    if is_stop(cell):
        return None, f"truncating variant in allele '{cell}'"
    match = PROTEIN_VARIANT.search(cell)
    if match is None:
        return None, f"no protein-level variant in allele '{cell}'"
    return int(match.group(1)), None


def is_stop(cell: str) -> bool:
    """A stop codon or a frameshift: the protein is truncated."""
    coding = CODING_VARIANT.search(cell)
    return (coding is not None and coding.group(3) in STOP_SUBSTITUTIONS) or FRAMESHIFT.search(cell) is not None


def missense_position(cell: str) -> Tuple[Optional[int], Optional[str]]:
    """Protein position of the single missense variant in an allele cell, or (None, reason to skip)."""
    position, reason = allele_position(cell)
    if position is None:
        return None, reason
    if CODING_VARIANT.search(cell) is None:
        return None, f"allele '{cell}' is not a three-letter amino acid substitution"
    return position, None


def protein_change(cell: str) -> str:
    match = re.search(r"\bp\.\S+", cell)
    return match.group(0).rstrip(",;") if match else cell


def allele_pairs(records: List[PatientRecord]) -> Tuple[List[AllelePair], List[Lollipop]]:
    """The patients' allele pairs, and the stop / missense patients as lollipops at the missense position."""
    pairs: List[AllelePair] = []
    lollipops: List[Lollipop] = []
    for record in records:
        patient = record.patient

        def skip(reason: str) -> None:
            note(f"patient {patient}: skipped - {reason}")

        if record.other_variants:
            skip("other variants present")
            continue

        allele1, allele2 = record.allele1, record.allele2
        if record.homozygous:
            if allele2:
                skip("allele #1 is homozygous, but allele #2 is also given")
                continue
            allele2 = allele1
        if not allele1 or not allele2:
            skip("only one allele given")
            continue

        # one truncating allele and one missense allele: a lollipop at the missense position
        stop1, stop2 = is_stop(allele1), is_stop(allele2)
        if stop1 and stop2:
            skip("both alleles are truncating")
            continue
        if stop1 or stop2:
            stop_allele, missense_allele = (allele1, allele2) if stop1 else (allele2, allele1)
            if len(CODING_VARIANT.findall(stop_allele)) > 1 or len(PROTEIN_VARIANT.findall(stop_allele)) > 1:
                skip(f"more than one variant in allele '{stop_allele}'")
                continue
            first, reason = missense_position(missense_allele)
            second = first
        else:
            first, reason = allele_position(allele1)
            if first is not None:
                second, reason = allele_position(allele2)
        if first is None or second is None:
            skip(reason)
            continue

        if record.phenotype is None:
            skip(f"diagnosis '{record.diagnosis}'")
            continue

        age = record.age
        label = (f"Patient {patient}: {protein_change(allele1)} / {protein_change(allele2)} "
                 f"({record.diagnosis}, age {age or 'unknown'})")
        above = record.phenotype == Phenotype.NONSYNDROMIC
        if stop1 or stop2:
            lollipops.append(Lollipop(patient, first, above, age, age != SOLID_AGE, label))
        else:
            pairs.append(AllelePair(patient, first, second, above, age, age != SOLID_AGE, label))
    return pairs, lollipops


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw allele pairs from patient tables (.tsv, .xlsx, .csv) as arcs over the protein, into an .svg."
    )
    tables_help = "patient tables: 01_pdf_table_parser.py .tsv, Moorfields .xlsx, Ambrosio .csv"
    parser.add_argument("tables", nargs="+", help=tables_help)
    parser.add_argument(
        "-m", "--mark", nargs=3, metavar=("POS1", "POS2", "SIDE"), default=[],
        help="a pair of protein positions to draw as an amber arc (a loop if equal), on the auditory side "
        "(SIDE 'a', above) or the USH1B side (SIDE 'u', below); the axis is extended to include them",
    )
    protein_help = f"draw the protein as a strip of domain boxes or as a bare line (default {DEFAULT_TRACK})"
    parser.add_argument("-p", "--protein", choices=PROTEIN_TRACKS, default=DEFAULT_TRACK, help=protein_help)
    parser.add_argument("-f", "--force", action="store_true", help="overwrite an existing output file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        paths, records = read_patient_tables(args.tables, prefix_ids=True)
    except InputError as exc:
        die(str(exc))
    except OSError as exc:
        die(f"could not read the patient tables: {exc}")

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

    out_path = output_path(paths, ".svg")
    if out_path.exists() and not args.force:
        die(f"output file '{out_path}' already exists (use --force to overwrite)")
    if out_path.exists() and not out_path.is_file():
        die(f"output path '{out_path}' exists and is not a regular file")

    pairs, lollipops = allele_pairs(records)
    if not pairs and not lollipops:
        die("no patient has a usable pair of protein-level variants", exit_code=2)

    genes = sorted({record.gene for record in records} - {""})
    gene = "/".join(genes) if genes else "protein"
    title = f"{gene}: protein positions of the two alleles per patient"
    sources = ", ".join(path.name for path in paths)
    subtitle = f"Source: {sources}. Hover an arc for the patient and variants."
    figure = Figure(gene, title, subtitle, pairs=pairs, lollipops=lollipops, marks=marks, mark_above=mark_above)
    try:
        write_text(svg_document(figure, PROTEIN_TRACKS[args.protein]), out_path)
    except OSError as exc:
        die(f"could not write '{out_path}': {exc}")
    n_skipped = len(records) - len(pairs) - len(lollipops)
    summary = (f"Wrote {len(pairs)} allele pairs and {len(lollipops)} stop / missense lollipops "
               f"({n_skipped} patients skipped) to {out_path}")
    print(summary, file=sys.stderr)

    # the patients that made it into the figure, on stdout
    print("Patient #\tAge (range)")
    for patient in pairs + lollipops:
        print(f"{patient.patient}\t{patient.age or 'unknown'}")


if __name__ == "__main__":
    main()
