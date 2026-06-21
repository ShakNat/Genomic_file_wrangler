# Genomic_file_wrangler
Collection of scripts to process and wrangle various genomic files.

**(1) get_cds_from_assembly.py** - Python script to obtain CDS file from assembly and GFF files

Usage instructions of **get_cds_from_assembly.py**:

```
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

```
**(2) transcript_to_gene_map.py** - Python script to obtain list of genes corresponding to list of transcripts by mapping with the GFF file

Usage instructions of **transcript_to_gene_map.py**:

```
python3 transcript_to_gene_map.py
					--gff <full path to GFF file>
					--transcripts <full path to TXT file with transcripts (one per line)>
					--genes <full path to output file>

```
**(3) Isoform_cleaner.py** - Python script to clean laternate transcripts from CDS file and produce representative CDS, PEP files with only the primary transcripts

Usage instructions of **Isoform_cleaner.py**:

```
python3 Isoform_cleaner.py
                    --cds <Full path to CDS file>
			        --sample_name <Name of the sample you are analyzing>
			        --out <Full path to output directory>
```
