# Arc diagram
 
## Is there a name for the type fo diagram produced in `03_allele_pair visualization.py`?
 
 Yes. It's called an arc diagram: points sit on a single line and arcs connect related pairs. Martin Wattenberg described it in "Arc Diagrams: Visualizing Structure in Strings" (IEEE InfoVis, 2002).                      
                                                                                                                                                                                                                             
  Related names you may see:                                                                                                                                                                                                 
                                                                                                                                                                                                                             
  - Linear arc plot / link plot: common in genomics. RNA secondary structure (base pairs drawn as arcs over the sequence) and Hi-C contacts are often shown this way.                                                        
  - Two-sided or mirrored arc diagram: arcs on both sides of the line, as in your script. Each side shows a different category of pair, here ARNSHL above and USH1B below.                                                   
  - Chord diagram: the circular version, with the line bent into a circle (as in Circos).                                                                                                                                    
                                                                                                                                                                                                                             
  The loops for homozygous pairs aren't a standard part of arc diagrams; they're a self-link convention borrowed from graph drawing. For a figure legend, something like "arc diagram of compound-heterozygous allele pairs  
  along MYO7A; loops denote homozygous variants" would describe it accurately.      
## Why does the subtitle say "Hover an arc for the patient and variants"?

The arcs have no text labels; with many overlapping arcs those would be unreadable. Instead each arc is an SVG `<path>` with a `<title>` child, which browsers show as a native tooltip on hover, e.g.

    Patient 207: p.Arg212His / p.Gly214Arg (ARNSHL, age 0-10)

(patient number, the two protein changes, diagnosis, and age range). The `--mark` arc ("Query pair: residues X / Y") and the dots on the protein line ("residue N") have tooltips too.

The tooltips only work when the SVG is viewed interactively (opened directly in a browser). Embedded as a static image in a paper or slide, or converted to PDF/PNG, they are gone, so for a print figure the hover sentence should be dropped from the subtitle.

## How can I share the interactive version with collaborators?

The SVG file itself is the interactive version; it only needs to be opened in a web browser.

- Simplest: send `data/sloan-heggen.suppl1.svg` and ask collaborators to download it and open it in Chrome, Firefox, Safari, or Edge (drag it into a browser window, or right-click → Open with). It works offline, nothing to install.
- For a link instead of an attachment, put the SVG on any static web host (GitHub Pages, university web space, Netlify Drop); a direct link to the `.svg` opens interactively.

What breaks the interactivity:

- Previews in email, Google Drive, Dropbox, or Slack usually render the SVG as a flat image; the file has to be downloaded and opened in a browser.
- GitHub's file view shows SVGs as images, and the "raw" link usually serves the source as text.
- Image viewers and Office apps (Preview, PowerPoint, Word, ...) don't show the tooltips.
- Phones and tablets have no hover; a long press may or may not show the tooltip.

Tips:

- Tooltips appear after holding the mouse still over an arc for about half a second; mention this, or people may think they don't work.
- For touch-device or print readers, send along the patient/age list the script prints on stdout (`... > included.tsv`).
