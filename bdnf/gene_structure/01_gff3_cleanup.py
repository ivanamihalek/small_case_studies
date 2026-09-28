#!/usr/bin/python3

# this is made for consumption by http://genometools.org/cgi-bin/gff3validator.cgi
# and http://genometools.org/cgi-bin/annotationsketch_demo.cgi

# 'gold' tagget transctipts coincidde wiht manual annotation (HAVANA)
# (from https://m.ensembl.org/Help/View?id=143)

root = 'ENSG00000176697.19'
parents = [root]
havana = ['ENST00000439476.6', 'ENST00000525528.1', 'ENST00000395986.6',
		  'ENST00000356660.9', 'ENST00000418212.5', 'ENST00000438929.5',
		  'ENST00000532997.5', 'ENST00000395981.7', 'ENST00000314915.6']

exons = {}
with open('orig.bdnf.gff3') as inf:
	for line in inf:
		if line[0]=='#':
			print(line.strip())
			continue
		# field = dict([a.split("=") for a in line.strip().split()[0].split(";")])
		field = line.strip().split()
		if len(field) == 0 or field[0]!='11': continue
		description = dict([a.split("=") for a in field[-1].split(";")])
		if description['Name']== root:
			print(line.strip())
			continue
		if 'biotype' in description and description['biotype']!='protein_coding':
			continue
		if 'ENST' in description['Name'] and description['Name'] not in havana:
			continue
		if description.get('Parent', 'no parent') in parents:
			if field[2]=='CDS': field[7]='0' # some crap about  phase
			print('\t'.join(field))
			parents.append(description['Name'])
			# append exons at the bottom as children of the gene, so we can show them collapsed all on the gene
			if field[2]=='exon':
				exons[description['Name']]=field
			continue

	for exon, field in exons.items():
		description = dict([a.split("=") for a in field[-1].split(";")])
		description['Parent'] = root
		field[-1]= ";".join(["{}={}".format(k,v) for k,v in description.items()])
		print('\t'.join(field))