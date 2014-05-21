
The scripts in this directory generate a multiple alignment across
all complete HIV-1 subtype B genomes in the LANL HIV genomic database
taken between the years of 1980 and 2014 in the United States.  A
complete alignment is established across all samples, about 1,100 at
the time of writing, and generates 100% consensus sequences for each
year to find conserved sections of the genome.

The file hiv-db-ALL.fasta contains the 1,100 source sequences from
which the scripts start work.

To prepare data automatically, just run ``make all'' on the script
directory.  Run ``make clean'' to clean up.  The makefile will call
the following scripts:

1.  The script process.py processes the source file, creating a
    "sequences" subdirectory and partitioning the samples in the
    source file by year.

2.  The script consensus.py determines the longest aligned 
    sequence in the collection and extends all sequences to that
    length by padding them on the left with dashes.  This makes all
    of the sequences in the multiple alignment the same length, so
    BioAlign can be used to generate a dumb consensus sequence.
    A set of 100%, 90%, 80%, and 70% consensus sequences are then
    generated for each year and written to a consensus sequence file.

3.  The script reduce.py examines all of the original sequences and
    calculates how many of those sequences contribute to the extended
    alignment at every position.  For every position in the extended
    alignment for which fewer than some number of sequences (default
    10) contribute to the extended alignment, that position of the 
    alignment will be deleted from all consensus sequences.  The reason
    for this is that, for the inital 1,100 sequences, the extended
    alignment was 12,526bp, which made it difficult to reconcile the
    consensus sequence with the protein map (at the time of writing,
    available at http://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html)

4.  The script import.py reads the reduced sequences and the reduced
    consensus sequence files and generates a MySQL 5.6 script suitable
    for importing the sequence and consensus data into a MySQL database
    that has been prepared with the file schema/hiv_schema.sql.
    The output of this script is captured as schema/hiv_data.sql.

5.  The script alignment.py reads the expanded sequences and generates
    a map of the alignment based on the number of reduced sequences that
    have a residue at every alignment index.  This is important because we
    need to "reduce" the alignment so that it falls into agreement with
    known protein maps of the HIV-1 subtype B genome.

Key to generated files:

    sequences/YYYY.fasta:
    Contains all sequences from year YYYY extracted from source HIV data.

    sequences/YYYY.fasta.extended:
    Contains all sequences in YYYY.fasta extended to a global alignment
    (all years, not just year YYYY).

    sequences/YYYY.fasta.extended.consensus:
    70%, 80%, 90%, 100% consensus sequences for all sequences in
    YYYY.fasta.extended.

    sequences/YYYY.fasta.extended.reduced: 
    Contains all sequences in YYYY.fasta.extended, but with poorly-covered
    alignment sequences removed.

    sequences/YYYY.fasta.extended.consensus.reduced:
    Contains all sequences in YYYY.fasta.extended.consensus, but with
    poorly-covered alignment segments removed.

    sequences/summary.consensus:
    Contains all sequences in YYYY.fasta.extended.consensus for all years.

    sequences/summary.consensus.reduced:
    Contains all sequences in summary.consensus, but with poorly-covered
    alignment segments removed.

    ./alignment_map.html:
    Contains a "heat map" of the number of reduced sequences (all sequences
    from YYYY.fasta.extended.reduced for all years) that participate at
    every index in the global alignment.  Originally, this was an attempt
    to visually identify the segments of the global alignment that were
    poorly covered by the original sequences; now, it's more of a check that
    we have not removed too much or too little of the global alignment 
    in the reduce.py script.
