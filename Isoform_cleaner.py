### Shakunthala Natarajan ###
### bug reports: s64snata@uni-bonn.de ###

__version__=0.1
__usage__="""
			python3 Isoform_cleaner.py
			--cds <Full path to CDS file>
			--gff <Full path to GFF file>
			--gff_config <Full path to GFF config file>
			--sample_name <Name of the sample you are analyzing>
			--out <Full path to output directory>
			"""

### --- start imports --- ###
import os,sys,glob, re
import gzip
### --- end of imports --- ###

def load_multiple_fasta_file(fasta_file):
	"""Load all sequences from a (possibly wrapped) FASTA file into a dict."""
	content = {}
	header = None
	seq_chunks = []
	#uncompressed FASTA file
	if fasta_file[-2:].lower() != 'gz':
		with open(fasta_file, "r") as f:
			for line in f:
				if not line.strip():
					continue
				if line.startswith(">"):
					if header is not None:
						content[header] = "".join(seq_chunks)
					header = line[1:].strip().split()[0]  # trim after first whitespace
					seq_chunks = []
				else:
					seq_chunks.append(line.strip())
			if header is not None:
				content[header] = "".join(seq_chunks)
	else:
		with gzip.open(fasta_file, "rt") as f:
			for line in f:
				if not line.strip():
					continue
				if line.startswith(">"):
					if header is not None:
						content[header] = "".join(seq_chunks)
					header = line[1:].strip().split()[0]  # trim after first whitespace
					seq_chunks = []
				else:
					seq_chunks.append(line.strip())
			if header is not None:
				content[header] = "".join(seq_chunks)
	return content

def translate(seq, genetic_code, unknown="X", internal_stop_to_x=False, keep_terminal_stop=True):
	"""Translate DNA to protein.
	   - Unknown/ambiguous codons -> 'X'
	   - Stop codons per table -> '*'
	   - Optionally convert internal '*' to 'X' (keeps final '*' if present).
	"""
	seq = seq.upper().replace("U", "T")
	pep = []
	n = len(seq) // 3
	for i in range(n):
		codon = seq[i*3:i*3+3]
		aa = genetic_code.get(codon)
		if aa is None:
			aa = unknown
		pep.append(aa)
	pep = "".join(pep)

	if internal_stop_to_x and "*" in pep:
		if keep_terminal_stop and pep.endswith("*"):
			pep = pep[:-1].replace("*", "X") + "*"
		else:
			pep = pep.replace("*", "X")
	return pep

def translate_file(in_fa, out_fa, genetic_code, internal_stop_to_x=False):
	seqs = load_multiple_fasta_file(in_fa)
	internal_stop_count = 0
	with open(out_fa, "w") as out:
		for header, nt in seqs.items():
			pep = translate(nt, genetic_code, unknown="X",
							internal_stop_to_x=internal_stop_to_x,
							keep_terminal_stop=True)
			if "*" in pep[:-1]:
				internal_stop_count += 1
			out.write(f">{header}\n{pep}\n")
	print(f" {os.path.basename(in_fa)} -> {os.path.basename(out_fa)} | " f"seqs: {len(seqs)} | internal-stops (pre-fix): {internal_stop_count}")

def gather_inputs(spec):
	"""Return list of input files. If directory, grab common FASTA extensions."""
	if os.path.isdir(spec):
		exts = ("*.fa", "*.fasta", "*.fna", "*.fas")
		files = []
		for ext in exts:
			files.extend(glob.glob(os.path.join(spec, ext)))
		files.sort()
		return files
	else:
		return [spec]

def make_output_path(out_spec, in_file):
	"""Decide output path: if out_spec is a dir or endswith '/', write there; else treat as file path."""
	if out_spec.endswith(os.sep) or os.path.isdir(out_spec):
		os.makedirs(out_spec, exist_ok=True)
		base = os.path.basename(in_file)
		root = os.path.splitext(base)[0]
		return os.path.join(out_spec, f"{root}.pep.fa")
	else:
		# Single input -> exact output file path
		# Multi-input with a single-file out_spec is not supported
		return out_spec

def load_IDs(filename):
	IDs = []
	if not filename.lower().endswith('.gz'):
		with open(filename, "r") as f:
			line = f.readline()
			while line:
				ID = line.strip()
				if len(ID) > 3:
					if "\t" in line:
						tmp = line.strip().split('\t')
						for each in tmp:
							if len(each) > 3:
								IDs.append(each)
					else:
						IDs.append(ID)
				line = f.readline()
	else:
		with gzip.open(filename, "rt") as f:
			line = f.readline()
			while line:
				ID = line.strip()
				if len(ID) > 3:
					if "\t" in line:
						tmp = line.strip().split('\t')
						for each in tmp:
							if len(each) > 3:
								IDs.append(each)
					else:
						IDs.append(ID)
				line = f.readline()
	return IDs

#functions to remove isoforms
def load_fasta(fasta_file):
	"""! @brief load FASTA alignment into dictionary	"""

	sequences = {}
	with open(fasta_file) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update({header: "".join(seq)})
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
				seq = []
			else:
				seq.append(line.strip())
			line = f.readline()
		sequences.update({header: "".join(seq)})
	return sequences

#function for removing alternate transcripts from the peptide FASTA file
def isoform_clean(gff3_input_file, cds_dict, no_trans_cds, child_attribute, child_parent_linker):
	no_gene_no_parent = False
	has_gene = False
	has_parent =False
	if gff3_input_file[-2:].lower() != 'gz':
		with open(gff3_input_file, "r") as f:
			gff_lines = f.readlines()
			# checking if gene feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)

	else:
		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			# checking if gene feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_mrna:
		coding_feature = 'MRNA'
	else:
		has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
		if has_transcript:
			coding_feature = 'TRANSCRIPT'
		else:
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			if has_cds:
				coding_feature = 'CDS'
			else:
				coding_feature = 'EXON'

	nogene_noparent_counter = 0
	if gff3_input_file[-2:].lower() != 'gz':  # uncompressed gff file
		transcripts_per_gene = {}
		with open(gff3_input_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					no_gene_no_parent = False  # Reset for each line
					parts = line.strip().split('\t')
					if len(parts) > 2:
						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker) + '=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								for each in partsnew:
									pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
									if re.match(pattern_par, each):
										partsnew1 = str(each).replace(str(child_parent_linker) + '=', "")
							for every in partsnew:
								pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
								if re.match(pattern_ID, every):
									partsnew0 = str(every).replace(str(child_attribute) + '=', "")
							try:
								transcripts_per_gene[partsnew1].append(partsnew0)
							except KeyError:
								transcripts_per_gene.update({partsnew1: [partsnew0]})
				line = f.readline()

	else:#compressed gff file
		transcripts_per_gene = {}
		with gzip.open(gff3_input_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					no_gene_no_parent = False  # Reset for each line
					parts = line.strip().split('\t')
					if len(parts) > 2:

						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker)+'=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								for each in partsnew:
									pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
									if re.match(pattern_par, each):
										partsnew1 = str(each).replace(str(child_parent_linker)+'=', "")
							for every in partsnew:
								pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
								if re.match(pattern_ID, every):
									partsnew0 = str(every).replace(str(child_attribute)+'=', "")
							try:
								transcripts_per_gene[partsnew1].append(partsnew0)
							except KeyError:
								transcripts_per_gene.update({partsnew1: [partsnew0]})
				line = f.readline()

	gene_names = list(transcripts_per_gene.keys())
	with open(no_trans_cds, "w") as out:
		for gene in gene_names:
			trans_length = []
			isoform_list = []
			for trans in transcripts_per_gene[gene]:
				if trans in cds_dict:
					if len(transcripts_per_gene[gene]) < 2:
						out.write('>' + str(trans) + '\n' + str(cds_dict[trans]) + '\n')
					else:
						trans_length.append((trans, cds_dict[trans]))
						isoform_list.append(trans)
			if trans_length:
				best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
				out.write('>' + str(best_trans) + "\n" + str(seq) + "\n")
				isoform_list.remove(best_trans)



def main(arguments):
	if '--sample_name' in arguments:
		orgname = arguments[arguments.index('--sample_name')+1]
	else:
		orgname = 'sample'

	cds_file = arguments[arguments.index('--cds') + 1]  # full path to CDS file

	if '--gff' in arguments:
		gff_file=arguments[arguments.index('--gff')+1]
	# gff file config params
	if '--gff_config' in arguments:
		gff_config_file = arguments[arguments.index('--gff_config') + 1]
		with open(gff_config_file, 'r') as f:
			for line in f:
				parts = line.strip().split()
				child_attribute = parts[0]
				child_parent_linker = parts[1]
				parent_attribute = parts[2]
	else:
		child_attribute = 'ID'
		child_parent_linker = 'Parent'
		parent_attribute = 'ID'

	outdir = arguments[arguments.index('--out') + 1]
	if outdir[-1] != "/":
		outdir += "/"
	if not os.path.exists(outdir):
		os.makedirs(outdir)


	#optional removal of isoforms
	isoform_reduced_cds_file = outdir + f"{orgname}_repr.cds.fasta"
	repr_output_spec = os.path.join(outdir, f'{orgname}_repr.pep.fasta')

	if os.path.exists(isoform_reduced_cds_file) and os.path.exists(repr_output_spec):
		pass
	else:
		print("Removing alternative isoforms")
		cds_dict = load_multiple_fasta_file(cds_file)
		isoform_clean(gff_file, cds_dict, isoform_reduced_cds_file, child_attribute, child_parent_linker)

		#code block to produce PEP file without isoforms
		repr_input_spec = isoform_reduced_cds_file
		internal_stop_to_x = ('--internal-stop-to-x' in arguments)

		# Standard code table (nuclear, 1)
		genetic_code = {
			'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R',
			'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I',
			'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q',
			'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P',
			'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
			'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S',
			'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V',
			'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
		}

		inputs = gather_inputs(repr_input_spec)

		# If multiple inputs but output_spec is a file path, abort to avoid overwriting
		if len(inputs) > 1 and (not repr_output_spec.endswith(os.sep) and not os.path.isdir(repr_output_spec)):
			sys.exit("[ERROR] For multiple input files, --out must be a directory or end with '/'.")
		if repr_output_spec.endswith(os.sep):
			os.makedirs(repr_output_spec, exist_ok=True)

		for in_file in inputs:
			pep_file = make_output_path(repr_output_spec, in_file)
			os.makedirs(os.path.dirname(pep_file), exist_ok=True)
			translate_file(in_file, pep_file, genetic_code, internal_stop_to_x=internal_stop_to_x)

if '--cds' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)
