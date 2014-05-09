
The scripts in this directory generate a multiple alignment across
all complete HIV-1 subtype B genomes in the LANL HIV genomic database
taken between the years of 1980 and 2014 in the United States.  A
complete alignment is established across all samples, about 1,100 at
the time of writing, and generates 100% consensus sequences for each
year to find conserved sections of the genome.

The file hiv-db-ALL.fasta contains the 1,100 source sequences from
which the scripts start work.

1.  The script 1_process.py processes the source file, creating a
    "sequences" subdirectory and partitioning the samples in the
    source file by year.

2.  The script 2_consensus.py determines the longest aligned 
    sequence in the collection and extends all sequences to that
    length by padding them on the left with dashes.  This makes all
    of the sequences in the multiple alignment the same length, so
    BioAlign can be used to generate a dumb consensus sequence.
    A dumb 100% consensus sequence is then generated for each year
    and written to a consensus sequence file.

3.  The script 3_import.py reads the expanded sequences and the
    consensus sequence file and generates a MySQL 5.6 script suitable
    for importing the sequence and consensus data into a MySQL database
    that has been prepared with the file schema/hiv_schema.sql.
    The output of this script is captured as schema/hiv_data.sql.
