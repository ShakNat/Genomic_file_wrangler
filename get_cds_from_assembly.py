### V0.1 ###
### Shakunthala Natarajan ###
"""
feature requests and bug reports to s64snata@uni-bonn.de
"""

__version__='0.1'
__usage__= """
			python3 get_cds_from_assembly.py
			--assembly <full path to file or folder of assembly files>
			--gff <full path to file or folder of gff files>
			--out <full path to output directory>
			
			NOTE_1: Assembly files should have extension of .fasta or .fasta.gz
			NOTE_2: GFF file should have extension of .gff or .gff.gz
			
			optional:
					--gff_config <full path to txt file containing the gff parameters to be used separated by spaces or tabs;
					Needs 4 columns in the order 
					(i) base file name - same as the base name you use for the gff file | all in case all the gff 
						have the same attribute pattern
					(ii) child_attribute: attribute field of the mRNA or transcript feature in the file like ID
					iii) child_parent_linker: attribute field of the mRNA or transcript, CDS, exon features that link them with their 
					respective parent feature like Parent - Note: base assumption by the tool is that all child levels
					have the same child-parent linker attribute fields. For eg., if Parent is the child-parent linker in the mRNA feature line,
					then Parent will be the child-parent linker for all other child-level feature lines in the GFF
					(iv) parent_attribute: attribute field of the gene feature like ID 
			"""

### --- begin imports --- ###
import os, sys, glob, re, subprocess, gzip
from operator import itemgetter


### --- end of imports --- ###

#Loading gene IDs from FASTA file
def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	if multiple_fasta_file[-2:].lower() != 'gz':#dealing with uncompressed fasta file
		sequences = {}
		with open(multiple_fasta_file) as f:
			header = f.readline()[1:].strip()
			if " " in header:
				header = header.split(' ')[0]
				if "\t" in header:
					header = header.split('\t')[0]
			elif "\t" in header:
				header = header.split('\t')[0]
			seq = []
			line = f.readline()
			while line:
				if line[0] == '>':
					sequences.update({header: "".join(seq).upper()})
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
						if "\t" in header:
							header = header.split('\t')[0]
					elif "\t" in header:
						header = header.split('\t')[0]
					seq = []
				else:
					seq.append(line.strip())
				line = f.readline()
			sequences.update({header: "".join(seq).upper()})
		return sequences

	else:#dealing with compressed FASTA file(s)
		sequences = {}
		with gzip.open (multiple_fasta_file, "rt") as f:
			header = f.readline()[1:].strip()
			if " " in header:
				header = header.split(' ')[0]
				if "\t" in header:
					header = header.split('\t')[0]
			elif "\t" in header:
				header = header.split('\t')[0]
			seq = []
			line = f.readline()
			while line:
				if line[0] == '>':
					sequences.update({header: "".join(seq).upper()})
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
						if "\t" in header:
							header = header.split('\t')[0]
					elif "\t" in header:
						header = header.split('\t')[0]
					seq = []
				else:
					seq.append(line.strip())
				line = f.readline()
			sequences.update({header: "".join(seq).upper()})
		return sequences

#function for loading transcript information from GFF3 file
def load_transcript_information_from_gff3( gff3_input_file,process_pseudos,child_attribute,child_parent_linker):
	"""! @brief load all transcript information from gff3 file """
	# --- load all data from file --- #
	message = []
	gff_pseudos = set()

	if gff3_input_file[-2:].lower() != 'gz':  # dealing with uncompressed gff file
		information = []
		mrna_dict = {}  # To store mRNA ID -> details mapping
		exon_dict = {}  # To store exon parent -> list of exons mapping
		#code to check for pseudogenes
		with open(gff3_input_file, "r") as f:
			gff_lines = f.readlines()
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
					  if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
					  if not line.startswith('#') and len(line.split('\t')) >= 3)
			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)
		if process_pseudos == 'no':
			if has_mrna:  # checking for pseudogenes when mrna feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith((child_attribute)+'='):
										gene_id = attr[(len(child_attribute)+1):]
										gff_pseudos_genes.add(gene_id)
										break
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'MRNA':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[(len(child_attribute)+1):]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[(len(child_parent_linker)+1):]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			elif has_transcript:  # checking for pseudogenes when transcript feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[3:]
										gff_pseudos_genes.add(gene_id)
										break

					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'TRANSCRIPT':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[(len(child_attribute)+1):]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[(len(child_parent_linker)+1):]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[
										8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			else:  # checking for pseudogenes when no mrna, transcript feature is present and only cds feature is present
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE' or 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[(len(child_attribute)+1):]
										gff_pseudos.add(gene_id)
										break

		with open(gff3_input_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) > 5:
						if has_cds:
							if parts[2].upper() == 'CDS':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None  # Changed from False to None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_parent_linker)+'=') in subp:
												parent = subp.replace((str(child_parent_linker)+'='), "")

										# Check if parent is pseudogene AFTER finding parent
										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue  # Skip this CDS

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if (str(child_parent_linker)+'=') in parts[-1]:
											parent = str(parts[-1]).replace((str(child_parent_linker)+'='), "")

										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
						else:
							# Handle MRNA/TRANSCRIPT and EXON features
							if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
								if len(parts[-1]) > len((str(child_attribute)+'=')):
									if ";" in parts[-1]:
										mrna_id = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_attribute)+'=') in subp:
												mrna_id = subp.replace((str(child_attribute)+'='), "")

										if mrna_id and mrna_id not in gff_pseudos:
											mrna_dict[mrna_id] = {
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
											}

							elif parts[2].upper() == 'EXON':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_parent_linker)+'=') in subp:
												parent = subp.replace((str(child_parent_linker)+'='), "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if (str(child_parent_linker)+'=') in parts[-1]:
											parent = str(parts[-1]).replace((str(child_parent_linker)+'='), "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
				line = f.readline()

		if has_cds:
			# --- sort data by parent --- #
			sorted_data = {}
			for each in information:
				try:
					sorted_data[each['parent']].append(each)
				except KeyError:
					sorted_data.update({each['parent']: [each]})

			final_data = []
			for key in sorted_data.keys():
				if sorted_data[key][0]['orientation'] == '+':
					final_data.append(sorted(sorted_data[key], key=itemgetter('start')))
				else:
					final_data.append(sorted(sorted_data[key], key=itemgetter('start'))[::-1])
		else:
			# Process exon/transcript based annotation
			final_data = []
			for mrna_id, mrna_info in mrna_dict.items():
				if mrna_id in exon_dict:
					exons = exon_dict[mrna_id]

					# Sort the exons based on orientation
					if mrna_info['orientation'] == '+':
						sorted_exons = sorted(exons, key=itemgetter('start'))
					else:
						sorted_exons = sorted(exons, key=itemgetter('start'), reverse=True)
					final_data.append(sorted_exons)
		return final_data, message

	else:#dealing with compressed GFF3 file(s)
		information = []
		mrna_dict = {}  # To store mRNA ID -> details mapping
		exon_dict = {}  # To store exon parent -> list of exons mapping


		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
				   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
						 if not line.startswith('#') and len(line.split('\t')) >= 3)
			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
				  if not line.startswith('#') and len(line.split('\t')) >= 3)

		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			# collecting pseudogenes to skip cds of pseudogenes
			for line in gff_lines:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
						for attr in parts[8].split(';'):
							if attr.startswith(str(child_attribute)+'='):
								gene_id = attr[len(str(child_attribute)+'=')+1:]
								gff_pseudos.add(gene_id)
								break

			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)

		if process_pseudos == 'no':
			if has_mrna:  # checking for pseudogenes when mrna feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[len(str(child_attribute)+'=')+1:]
										gff_pseudos_genes.add(gene_id)
										break
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'MRNA':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[len(str(child_attribute)+'=')+1:]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[len(str(child_parent_linker)+'=')+1:]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			elif has_transcript:  # checking for pseudogenes when transcript feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[len(str(child_attribute)+'=')+1:]
										gff_pseudos_genes.add(gene_id)
										break

					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'TRANSCRIPT':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[len(str(child_attribute)+'=')+1:]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[len(str(child_parent_linker)+'=')+1:]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[
										8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			else:  # checking for pseudogenes when no mrna, transcript feature is present and only cds feature is present
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE' or 'PSEUDO=TRUE' in parts[
								8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[len(str(child_attribute)+'=')+1:]
										gff_pseudos.add(gene_id)
										break

		with gzip.open(gff3_input_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) > 5:
						if has_cds:
							if parts[2].upper() == 'CDS':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None  # Changed from False to None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_parent_linker)+'=' in subp:
												parent = subp.replace(str(child_parent_linker)+'=', "")

										# Check if parent is pseudogene AFTER finding parent
										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue  # Skip this CDS

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if str(child_parent_linker)+'=' in parts[-1]:
											parent = str(parts[-1]).replace(str(child_parent_linker)+'=', "")

										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
						else:
							# Handle MRNA/TRANSCRIPT and EXON features
							if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
								if len(parts[-1]) > len(str(child_attribute)+'='):
									if ";" in parts[-1]:
										mrna_id = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_attribute)+'=' in subp:
												mrna_id = subp.replace(str(child_attribute)+'=', "")

										if mrna_id and mrna_id not in gff_pseudos:
											mrna_dict[mrna_id] = {
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
											}

							elif parts[2].upper() == 'EXON':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_parent_linker)+'=' in subp:
												parent = subp.replace(str(child_parent_linker)+'=', "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if str(child_parent_linker)+'=' in parts[-1]:
											parent = str(parts[-1]).replace(str(child_parent_linker)+'=', "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
				line = f.readline()

		if has_cds:
			# --- sort data by parent --- #
			sorted_data = {}
			for each in information:
				try:
					sorted_data[each['parent']].append(each)
				except KeyError:
					sorted_data.update({each['parent']: [each]})

			final_data = []
			for key in sorted_data.keys():
				if sorted_data[key][0]['orientation'] == '+':
					final_data.append(sorted(sorted_data[key], key=itemgetter('start')))
				else:
					final_data.append(sorted(sorted_data[key], key=itemgetter('start'))[::-1])
		else:
			# Process exon/transcript based annotation
			final_data = []
			for mrna_id, mrna_info in mrna_dict.items():
				if mrna_id in exon_dict:
					exons = exon_dict[mrna_id]
					# Sort the exons based on orientation
					if mrna_info['orientation'] == '+':
						sorted_exons = sorted(exons, key=itemgetter('start'))
					else:
						sorted_exons = sorted(exons, key=itemgetter('start'), reverse=True)
					final_data.append(sorted_exons)
		return final_data, message

#function of constructing the CDS FASTA file
def construct_CDS_file( transcript_info, CDS_file, assembly, child_parent_linker):
	"""! @brief construct file with all sequences for translation """
	with open( CDS_file, "w" ) as out:
		for transcript in transcript_info:
			seq = []
			revcomp_status = False
			if transcript[0]['orientation'] == '-':
				revcomp_status = True
			for part in transcript:
				if revcomp_status:
					seq.append( revcomp( assembly[ part['chr'] ][ part['start']-1:part['end'] ] ) )
				else:
					seq.append( assembly[ part['chr'] ][ part['start']-1:part['end'] ] )
				# Get parent ID, handling both formats of CDS feature or mRNA/ exon features' presence
				parent_id = transcript[0]['parent']
				if str(child_parent_linker)+"=" in parent_id:
					parent_id = parent_id.replace(str(child_parent_linker)+"=", "")
			out.write( '>' + str(parent_id) + '\n' + "".join( seq ) + '\n' )

#function for constructing the reverse complement
def revcomp( seq ):
	"""! @brief constructs revcomp """
	new_seq = []
	dictionary = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N','a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n' }
	for nt in seq:
		try:
			new_seq.append( dictionary[ nt ] )
		except KeyError:
			new_seq.append( "N")
	return ''.join( new_seq[::-1] )

def main (arguments):
	assembly_arg=arguments[arguments.index('--assembly')+1]#full path to file or folder of assembly files
	gff_arg=arguments[arguments.index('--gff')+1]#full path to file or folder of gff files
	outdir = arguments[arguments.index('--out') + 1]
	if outdir[-1] != "/":
		outdir += "/"
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	process_pseudos = 'yes'
	gff_config_parameters = {}
	if '--gff_config' in arguments:
		gff_config_file = arguments[arguments.index('--gff_config')+1]
		with open(gff_config_file,'r')as f:
			for line in f:
				parts=line.strip().split()
				if len(parts)==4:
					gff_config_parameters[parts[0]] = {'child_attribute': parts[1],'child_parent_linker': parts[2],'parent_attribute': parts[3]}
	else:
		gff_config_parameters['default'] = {
			'child_attribute': 'ID',
			'child_parent_linker': 'Parent',
			'parent_attribute': 'ID'
		}
	process_pseudos = 'yes'
	assembly_files=[]
	gff_files=[]
	if os.path.isfile(assembly_arg):
		assembly_files=[assembly_arg]
	elif os.path.isdir(assembly_arg):
		assembly_files = sorted(glob.glob(os.path.join(assembly_arg, "*")))
	if os.path.isfile(gff_arg):
		gff_files=[gff_arg]
	elif os.path.isdir(gff_arg):
		gff_files = sorted(glob.glob(os.path.join(gff_arg, "*")))
	for assembly in assembly_files:
		if not (assembly.lower()).endswith('gz'):
			assembly_name = os.path.basename(assembly).lower().replace('.fasta','').strip()
		else:
			assembly_name = os.path.basename(assembly).lower().replace('.fasta.gz', '').strip()
		for gff in gff_files:
			if not (gff.lower()).endswith('gz'):
				gff_name = os.path.basename(gff).lower().replace('.gff', '').strip()
			else:
				gff_name = os.path.basename(gff).lower().replace('.gff.gz', '').strip()
			if assembly_name == gff_name:
				cds_out = os.path.join(outdir,f"{gff_name}.cds.fasta")
				print(f"gff_name is {gff_name}")
				if gff_name in gff_config_parameters:
					child_attribute = gff_config_parameters[gff_name]['child_attribute']
					child_parent_linker = gff_config_parameters[gff_name]['child_parent_linker']
				else:
					child_attribute = gff_config_parameters['default']['child_attribute']
					child_parent_linker = gff_config_parameters['default']['child_parent_linker']
				seqs = load_sequences(str(assembly))
				transcript_information, message = load_transcript_information_from_gff3(str(gff), process_pseudos,child_attribute,child_parent_linker)
				construct_CDS_file(transcript_information, cds_out, seqs, child_parent_linker)

if '--assembly' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)