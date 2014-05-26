
The scripts in this directory generate a multiple alignment across
all complete HIV-1 genomes in the LANL HIV genomic database
taken between the years of 1980 and 2014 in the United States.  A
complete alignment is established across all samples, about 1,110 at
the time of writing, and generates consensus sequences of varying
thresholds for each year to find conserved sections of the genome.
At present, the analysis is performed on subtypes B and C, but it
is easy to add support for other subtypes.

The source files, acquired from the LANL HIV genomic database, are:

    ./sequences_B.fasta
    All complete genomic sequences for HIV-1 subtype B.

    ./sequences_C.fasta
    All complete genomic sequences for HIV-1 subtype C.

    ./patients_B.txt
    All patient information corresponding to the sequences for
    subtype B. 

    ./patients_C.txt
    All patient information corresponding to the sequences for
    subtype C. 

To prepare data automatically, just run ``make all`` on the script
directory.  Run ``make clean`` to clean up.  The makefile will call
the following scripts:

1.  The script process.py processes the source file, creating a
    "sequences" subdirectory and partitioning the samples in the
    source file by year.

    This script also prevents more than one sequence from any single
    patient from being included in any year.  This is necessary
    because some patients have many sequences listed in a given year
    and their presence can bias the consensus sequence for that
    year.  Since the patient information is not generally deducible
    in the sequence header, we need to consult a patient_info.txt
    file that connects the accession ID of each sequence to the
    patient code and year the sample was sequenced.

    The script supports analysis of multiple subtypes, but at present
    only subtype B is analyzed.

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
    is 1, based on the number and structure of remaining sequences)
    contribute to the extended alignment, that position of the 
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

    ./sequences/SSSS:
    Contains all generated analysis for subtype SSSS.

    ./sequences/SSSS/YYYY.fasta:
    Contains all sequences from year YYYY and subtype SSSS extracted from
    source HIV data.

    ./sequences/SSSS/YYYY.fasta.extended:
    Contains all sequences in YYYY.fasta and subtype SSSS extended to a
    global alignment (all years for subtype SSSS, not just year YYYY).

    ./sequences/SSSS/YYYY.fasta.extended.consensus:
    70%, 80%, 90%, 100% consensus sequences for all sequences in
    YYYY.fasta.extended.

    ./sequences/SSSS/YYYY.fasta.extended.reduced: 
    Contains all sequences in YYYY.fasta.extended, but with poorly-covered
    alignment sequences removed.

    ./sequences/SSSS/YYYY.fasta.extended.consensus.reduced:
    Contains all sequences in YYYY.fasta.extended.consensus, but with
    poorly-covered alignment segments removed.

    ./sequences/SSSS/summary.consensus:
    Contains all sequences in YYYY.fasta.extended.consensus for all years.

    ./sequences/SSSS/summary.consensus.reduced:
    Contains all sequences in summary.consensus, but with poorly-covered
    alignment segments removed.

    ./alignment_SSSS.html:
    Contains a "heat map" of the number of reduced sequences (all sequences
    from YYYY.fasta.extended.reduced for all years) that participate at
    every index in the global alignment.  Originally, this was an attempt
    to visually identify the segments of the global alignment that were
    poorly covered by the original sequences; now, it's more of a check that
    we have not removed too much or too little of the global alignment 
    in the reduce.py script.
