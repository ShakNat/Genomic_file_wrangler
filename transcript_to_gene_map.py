### Shakunthala Natarajan ###
### s64snata@uni-bonn.de ###
__version__ = "v0.1"

__usage__ = """
					Usage:
					python3 transcript_to_gene_map.py
					--gff <FULL_PATH_TO_GFF_FILE_INCLUDING_FILE_NAME>
					--transcripts <FULL_PATH_TO_TXT_FILE_WITH_TRANSCRIPT_NAMES_ONE_PER_LINE_INCLUDING_FILE_NAME>
					--genes <FULL_PATH_TO_OUTPUT_FILE_INCLUDING_FILE_NAME>
					"""


import re, os, sys, subprocess, gzip
import tempfile
import traceback
import numpy as np
from decimal import Decimal, ROUND_HALF_DOWN
from collections import defaultdict
try:
	import matplotlib.pyplot as plt
	from matplotlib.lines import Line2D
	import matplotlib.ticker as ticker
	from matplotlib.patches import FancyArrow
except ImportError:
	pass

# --- end of imports --- #


def load_gene_infos( gff_file, child_attribute, child_parent_linker, parent_attribute):
	"""! @brief load gene ID, position, and orientation from GFF3 file """
	
	gene_infos = {}
	mrna_infos = {}
	five_utr_infos={}
	cds_infos = {}
	genes_per_chromosome = {}
	transcripts_per_gene = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[2].upper() == "GENE" or parts[2].upper() == "TRANSPOSABLE_ELEMENT_GENE" or parts[2].upper() == "PSEUDOGENE" or parts[2].upper() == "PSEUDO_GENE":	#could be extended to other feature types
					ID = parts[-1].split(f'{parent_attribute}=')[-1]
					if ";" in ID:
						ID = ID.split(';')[0]
					gene_infos.update( { ID: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
					try:
						genes_per_chromosome[ parts[0] ].append( ID )
					except KeyError:
						genes_per_chromosome.update( { parts[0]: [ ID ] } )
				if parts[2].upper() == "MRNA":
					ID = parts[-1].split(f'{child_attribute}=')[-1]
					if ";" in ID:
						ID = ID.split(';')[0]
					Parent = parts[-1].split(f'{child_parent_linker}=')[-1]
					if ";" in Parent:
						Parent = Parent.split(';')[0]
					mrna_infos.update( { ID: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
					try:
						transcripts_per_gene[ Parent ].append( ID )
					except KeyError:
						transcripts_per_gene.update( { Parent: [ ID ] } )
				if parts[2].upper() == 'FIVE_PRIME_UTR':
					Parent = parts[-1].split(f'{child_parent_linker}=')[-1]#Parent of 5'UTR is transcript
					if ";" in Parent:
						Parent = Parent.split(';')[0]
					five_utr_infos.update({ Parent: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } })# key of this nested dictionary is the transcript name
				if parts[2].upper() == 'CDS':
					cds_parents = parts[-1].split(f'{child_parent_linker}=')[-1]
					if ";" in cds_parents:
						cds_parents = cds_parents.split(';')[0]
					for cds_parent in cds_parents.split(','):  # handle multiple parents
						cds_parent = cds_parent.strip()
						cds_tuple = (int(parts[3]), int(parts[4]))
						try:
							cds_infos[cds_parent].append(cds_tuple)
						except KeyError:
							cds_infos[cds_parent] = [cds_tuple]
			line = f.readline()
	for chromosome in genes_per_chromosome: #sort the genes in each contig/ chromosome in the ascending order of start positions
		genes_per_chromosome[chromosome].sort(key=lambda gene: (gene_infos[gene]['start'], gene))
	for gene in transcripts_per_gene:#sort the transcripts per gene in the ascending order of mRNA start positions for + strand and mRNA end positions for - strand
		if gene_infos[gene]['orientation']=='+':
			transcripts_per_gene[gene].sort(key=lambda transcript: (mrna_infos[transcript]['start'], transcript))
		elif gene_infos[gene]['orientation']=='-':
			transcripts_per_gene[gene].sort(key=lambda transcript: (mrna_infos[transcript]['end'], transcript))
	gene_atg_dic = {}
	for gene in transcripts_per_gene:
		if gene not in gene_infos:
			continue
		orientation = gene_infos[gene]['orientation']
		transcript_list = transcripts_per_gene[gene]

		# select most upstream transcript for +, most downstream for -
		if orientation == '+':
			selected_transcript = transcript_list[0]  # already sorted by start ascending
		else:
			selected_transcript = transcript_list[-1]  # already sorted by end ascending, last = highest end

		if selected_transcript not in cds_infos:
			continue  # no CDS annotated for this transcript

		cds_list = cds_infos[selected_transcript]

		if orientation == '+':
			cds_list.sort(key=lambda x: x[0])  # sort by start ascending
			gene_atg_dic[gene] = cds_list[0][0]  # start of most upstream CDS = ATG
		else:
			cds_list.sort(key=lambda x: x[1])  # sort by end ascending
			gene_atg_dic[gene] = cds_list[-1][1]  # end of most downstream CDS = ATG
	return gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos


def main( arguments ):
	"""! @brief run everything """
	
	#input: GFF file with gene positions (translation start sites)
	gff_file = arguments[ arguments.index('--gff')+1 ]

	#gff file config params
	if '--gff_config' in arguments:
		gff_config_file = arguments[arguments.index('--gff_config')+1]
		with open (gff_config_file, 'r') as f:
			for line in f:
				parts = line.strip().split()
				child_attribute = parts[0]
				child_parent_linker = parts[1]
				parent_attribute = parts[2]
	else:
		child_attribute = 'ID'
		child_parent_linker = 'Parent'
		parent_attribute = 'ID'

	#output file
	output_file = arguments[ arguments.index('--genes')+1 ]
	#file with transcript names one per line
	transcript_file = arguments[ arguments.index('--transcripts')+1 ]
	transcript_list = []
	with open (transcript_file, 'r') as f:
		for line in f:
			transcript = line.strip()
			transcript_list.append(transcript)
	gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos = load_gene_infos( gff_file, child_attribute, child_parent_linker, parent_attribute)
	transcript_gene_map = {}
	for gene in transcripts_per_gene:
		for transcript in transcript_list:
			if transcript in transcripts_per_gene[gene]:
				transcript_gene_map[transcript] = gene
	with open(output_file, 'w') as out:
		for transcript in transcript_gene_map:
			out.write(f'{transcript_gene_map[transcript]}\n')

if '--gff' in sys.argv and '--transcripts' in sys.argv and '--genes' in sys.argv:
	main(sys.argv)
else:
	sys.exit( __usage__ )
