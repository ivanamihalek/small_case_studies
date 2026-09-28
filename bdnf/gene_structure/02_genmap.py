#!/usr/bin/python
from sys import argv, stderr
from gt import *

def main():

	width = 2400

	if len(argv) != 4:
		stderr.write("Usage: " + argv[0] +
					 " <style_file> <outfile> <gff3_file(input)>\n")
		stderr.write("Create png/svg representation of gff3 annotation file.\n")
		exit(1)

	outfile = argv[2]
	svg = outfile[-4:]==".svg"

	# load style file
	style = Style()
	style.load_file(argv[1])

	# create feature index
	feature_index = FeatureIndexMemory()

	# add GFF3 file to index
	feature_index.add_gff3file(argv[3])

	# create diagram for first sequence ID in feature index
	seqid = feature_index.get_first_seqid()
	range = feature_index.get_range_for_seqid(seqid)
	diagram = Diagram.from_index(feature_index, seqid, range, style)

	# create layout

	layout = Layout(diagram, width, style)
	height = layout.get_height()

	# create canvas
	canvas = CanvasCairoFileSVG(style, width, height) if svg else CanvasCairoFile(style, width, height)

	# sketch layout on canvas

	layout.sketch(canvas)

	# write canvas to file

	canvas.to_file(outfile)

	return

if __name__ == "__main__":
	main()
