# Extracting Table 1 from `../literature/ambrosio.pdf`

Source: Ambrosio et al., Doc Ophthalmol (2021) 143:39–51, Table 1
("Genetic diagnosis and ophthalmic characteristics by study patients").

Output: `data/ambrosio_table1.tsv` (39 patients: 26 USH1, 12 USH2, 1 USH3).

The table is printed sideways and runs over journal pages 42–43 (PDF pages 4–5).
Rows 1–32 were parsed automatically from the PDF's text. Rows 33–39 were typed in by hand from the
page image, because the PDF's text for page 5 was scrambled.

Unless noted, commands are run from the project root.

## 1. Find which pages contain Table 1

Run from `data/`:

```bash
for p in $(seq 1 13); do echo "== $p"; pdftotext -f $p -l $p -layout ambrosio.pdf - | grep -n "Table 1" | head; done
```

The table starts on page 4 and continues on page 5.

## 2. Dump the text of those pages, keeping the column layout

```bash
pdftotext -f 4 -l 5 -layout ambrosio.pdf -
```

Page 4 comes out with clean columns. Page 5 is jumbled because its text is mixed with the body text that follows.

## 3. Render pages 4–5 as images and rotate them for reading

```bash
cd /tmp && pdftoppm -f 4 -l 5 -r 200 -png ~/academia/projects/small_case_studies/myo7a/data/ambrosio.pdf ambh
python3 -c "
from PIL import Image
for p in ['ambh-04.png','ambh-05.png']:
    im=Image.open(p).rotate(90, expand=True); im.save(p.replace('.png','_r.png')); print(im.size)
"
```

The rotated images `/tmp/ambh-04_r.png` and `/tmp/ambh-05_r.png` come out upside down but are still readable.
They were used to check rows 1–32, to find out what the footnote marks mean, and to transcribe rows 33–39.

## 4. Build the TSV

```bash
python3 - <<'EOF'
import re, subprocess
txt = subprocess.run(["pdftotext", "-f", "4", "-l", "4", "-layout", "data/ambrosio.pdf", "-"], capture_output=True, text=True).stdout
rows = []
for line in txt.splitlines():
    m = re.match(r"^\s+(\d+)\s{2,}(\S+)\s{2,}(.+)$", line)
    if not m or int(m.group(1)) > 39: continue
    rest = re.split(r"\s{2,}", m.group(3).strip())
    rows.append([m.group(1), m.group(2)] + rest)
# rows 33-39, transcribed by hand from the page-5 image (kept in the same garbled notation so one parser handles both)
page5 = """33|USH2A|c.6224G [ A; p.W2075*|c.12714 T [ G;p.Y4238*|8.9|Normal|Normal|Pigmentation|Altered pigment
34|USH2A|c.949C [ A;p.R317R|c.(?_5299)_(5572_?)del;Exon 27 deletion|10.1|Normal|Attenuated|Normal|Normal
35|USH2A|c.2299delG; p.E767fs|c.11231 ? 1G [ T§|10.7|Normal|Attenuated|Pigmentation|Altered pigment
36|USH2A|c.2299delG; p.E767fs|c.11231 ? 1G [ T§|13.8|Normal|Attenuated|Pigmentation|altered pigment
37|USH2A|c.11864G [ A; p.W3955*|c.3170A [ G;p.Q1057R§|17.5|Normal|Attenuated|Pigmentation|NORMAL
38|USH2A|c.12295-2A [ G§|c.1227G [ C;p.W409C§|23.1|Normal|Attenuated|Pigmentation|Normal
39|CLRN1|c.189C [ A;p.Y63*|c.189C [ A;p.Y63*|13.3|Mild pallor|Attenuated|Pigmentation|Altered pigment"""
rows += [l.split("|") for l in page5.splitlines()]
# footnote marks as they appear in pdftotext output: à = ‡, § = §, || = ‖; $ is a misread § (pt 29)
notes = {"à": "novel", "§": "clinvar_only", "||": "novel_prev_reported_PMID27631835", "$": "clinvar_only"}
def parse_mut(s):
    note = ""
    for k, v in notes.items():
        if s.endswith(k):
            note, s = v, s[: -len(k)]
    s = s.replace(" [ ", ">").replace(" ? ", "+").replace("c. ", "c.").replace(" T>", "T>")
    s = s.rstrip(")") if s.count(")") > s.count("(") else s
    parts = [p.strip() for p in s.split(";", 1)]
    cdna, prot = parts[0], (parts[1] if len(parts) > 1 else "")
    return cdna, prot, note
header = ["patient", "ush_type", "gene", "mutation1_cdna", "mutation1_protein", "mutation1_note", "mutation2_cdna",
          "mutation2_protein", "mutation2_note", "age_first_exam_years", "disc", "vessels", "fundus", "macula"]
out = ["\t".join(header)]
for r in rows:
    assert len(r) == 9, r
    pt = int(r[0])
    ush = "USH1" if pt <= 26 else ("USH2" if pt <= 38 else "USH3")
    m1, m2 = parse_mut(r[2]), parse_mut(r[3])
    out.append("\t".join([r[0], ush, r[1], *m1, *m2, *r[4:]]))
open("data/ambrosio_table1.tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
EOF
```

This is a lightly edited version of the script that was actually run:
- the two comments inside the script were added afterwards, to explain it;
- a no-op line, `prot = prot.replace("Q1758*", "Q1758*")`, was removed. The output is the same either way.

## Output columns

`patient`, `ush_type`, `gene`, then for each of the two mutations its DNA change (`_cdna`), protein change
(`_protein`) and footnote (`_note`), followed by `age_first_exam_years`, `disc`, `vessels`, `fundus`, `macula`.

## Changes from the printed table

- **Garbled symbols fixed:** the PDF's text has some symbols wrong: `[` is restored to `>`, `?` to `+`,
  and stray spaces are removed (e.g. `c. 999 T [ G` → `c.999T>G`).
- **Mutations split** at the semicolon into DNA and protein parts. Where the paper gives only a DNA change
  (splice variants like `c.2094+1G>A`), the protein column is empty.
- **Footnote marks turned into words** in the `_note` columns:
  - ‡ → `novel` (not reported in the literature or ClinVar)
  - § → `clinvar_only`
  - ‖ → `novel_prev_reported_PMID27631835` (new variant, but the patient was described in an earlier paper)

## Kept as printed

- Inconsistent capitalisation in the eye-exam columns ("altered pigment", "Altered Pigment", "NORMAL",
  "pigmentation").
- `p.Pro3804LeufsX13` (pt 30), which uses a different notation from the rest; only a stray `)` was dropped.
- `p.L196` (pt 17), which is missing the changed amino acid in the paper itself.

## Reproducibility

Rows 33–39 depend on a reading of the page image. If you re-run this, check those values by eye.
