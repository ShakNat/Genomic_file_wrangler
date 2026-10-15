# Genomic_file_wrangler
Collection of scripts to process and wrangle various genomic files.

**(1) get_cds_from_assembly.py** - Python script to obtain CDS file from assembly and GFF files

Usage instructions of **get_cds_from_assembly.py**:

```
python3 get_cds_from_assembly.py
    --assembly <full path to file or folder of assembly files>
    --gff <full path to file or folder of gff files>
    --out <full path to output directory>

    NOTE_1: Assembly files should have extension of .fa, .fna or .fasta (optionally followed by .gz)
    NOTE_2: GFF files should have extension of .gff or .gff3 (optionally followed by .gz)
    NOTE_3: An assembly file and a GFF file are paired when their names are identical after removing the extension (case ignored),
            e.g. species1.fasta.gz and species1.gff3
    NOTE_4: GFF files without CDS features are skipped with a warning, because exon-based sequences would include UTRs
    NOTE_5: GFF seqids (column 1) must match the FASTA headers (first word only). Transcripts on contigs that are missing
            in the FASTA are skipped; if no seqid matches at all, the species is skipped with an error message

    optional:
        --gff_config <full path to txt file containing the gff parameters to be used separated by spaces or tabs>
        Needs 4 columns in the order
        (i) base file name - same as the base name you use for the gff file (without extension) | all in case all the gff
            have the same attribute pattern
        (ii) child_attribute: attribute field of the mRNA or transcript feature in the file like ID
        (iii) child_parent_linker: attribute field of the mRNA or transcript, CDS, exon features that link them with their
              respective parent feature like Parent - Note: base assumption by the tool is that all child levels
              have the same child-parent linker attribute fields. For eg., if Parent is the child-parent linker in the
              mRNA feature line, then Parent will be the child-parent linker for all other child-level feature lines in the GFF
        (iv) parent_attribute: attribute field of the gene feature like ID
        Default if no config is given (or no matching row): ID, Parent, ID

        --process_pseudos <yes|no>
        yes (default): CDS of pseudogenes are extracted like all others
        no: CDS of pseudogenes are skipped

        --transsplice <warn|skip|process>
        Controls transcripts that look trans-spliced (CDS segments on different strands or carrying part= attributes)
        warn (default): print a warning and extract using the standard single-strand rules (sequence might be wrong)
        skip: print a warning and leave these transcripts out
        process: print a warning and extract them, ordering segments by their part= numbers and handling the strand of
                 each segment separately; transcripts without part= numbers on all segments are skipped

        --contig_boundary_violation_process <yes|no>
        Controls CDS segments whose GFF coordinates run past the end of the contig. NCBI writes genes that cross the
        origin of circular genomes (mitochondria, plastids) this way, e.g. a CDS segment at 394625-395019 on a 367808 bp contig
        yes: such segments wrap around to the start of the contig (circular contig assumed) and a warning is printed for each
        no (default): coordinates are used as given, so such segments come out empty or cut short without a warning;
                      check the summary for sequences whose length is not a multiple of 3

    Transcripts with unstranded features (strand other than + or -) are always skipped with a warning.
    CDS lines with several comma-separated parents are assigned to each of the parents.

    Output files in the output directory:
        <name>.cds.fasta - CDS sequences, one per transcript
        <name>.cds_summary.txt - counts of sequences, length not a multiple of 3, missing ATG start,
                                 missing terminal stop, internal stops, non-ACGT characters, duplicate IDs
        <name>.cds_problems.tsv - IDs and lengths of the flagged sequences with the issues found
        cds_summary_all_species.tsv - one summary row per species

```

**(2) compare_cds_seq.py** - Compares CDS FASTA files produced by two different tools (e.g. this repository's `get_cds_from_assembly.py` 
versus a reference such as the NCBI) for the same annotation. The comparison is done with respect to the reference: sequences are matched 
by the hash that seqkit calculates for each sequence, so sequence IDs are not used for matching.

Instructions to download seqkit can be found here: https://bioinf.shenwei.me/seqkit/download/

Usage instructions of **compare_cds_seq.py**:

    python3 compare_cds_seq.py
    --reference <full path to folder of reference CDS files>
    --test <full path to folder of CDS files to be tested>
    --out <full path to output directory>

    NOTE_1: CDS files should have the extension .cds.fa, .cds.fna or .cds.fasta (optionally followed by .gz)
    NOTE_2: A reference file and a test file are paired when their names are identical after removing the extension (case ignored),
            e.g. species1.cds.fasta and Species1.cds.fa.gz. Files without a partner in the other folder are skipped with a warning
    NOTE_3: The comparison is done with respect to the reference. Sequences are compared by the hash calculated by seqkit
            (seqkit fx2tab -s); sequences occurring more than once must occur equally often in both files
    NOTE_4: seqkit has to be installed: conda install -c bioconda seqkit

    All three arguments are mandatory; there are no optional arguments.

    Output files in the output directory:
            <name>_differences.tsv - sequences found only in the reference (reference_only) or only in the test file (test_only);
                                     columns: status, id, sequence_hash (only written if there are differences)
            comparison_summary.tsv - one row per name; columns: name, reference_total, test_total, identical, reference_only, test_only


**(3) transcript_to_gene_map.py** - Python script to obtain list of genes corresponding to list of transcripts by mapping with the GFF file

Usage instructions of **transcript_to_gene_map.py**:

```
python3 transcript_to_gene_map.py
					--gff <full path to GFF file>
					--transcripts <full path to TXT file with transcripts (one per line)>
					--genes <full path to output file>

```
**(4) Isoform_cleaner.py** - Python script to clean laternate transcripts from CDS file and produce representative CDS, PEP files with only the primary transcripts

Usage instructions of **Isoform_cleaner.py**:

```
python3 Isoform_cleaner.py
                    --cds <Full path to CDS file>
                    --gff <Full path to GFF file>
					--sample_name <Name of the sample you are analyzing>
					--out <Full path to output directory>

					Optional:
                    --gff_config <Full path to GFF config TXT file>
					Needs 3 columns in the order

					(i) child_attribute: attribute field of the mRNA or transcript feature in the file like ID

					(ii) child_parent_linker: attribute field of the mRNA or transcript, CDS, exon features that link them with their 
					respective parent feature like Parent - Note: base assumption by the tool is that all child levels
					have the same child-parent linker attribute fields. For eg., if Parent is the child-parent linker in the mRNA 						feature line,then Parent will be the child-parent linker for all other child-level feature lines in the GFF

					(iii) parent_attribute: attribute field of the gene feature like ID			        
```

##### Kindly cite this repository if you use any scripts from here. Thank you!

### Third party tool references:
1. Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191. 
