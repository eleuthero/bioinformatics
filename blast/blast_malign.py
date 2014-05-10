from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from Bio import AlignIO

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# WARNING:
# ========
# I'm just writing this to show you that it can be done, not that you
# should do it.  It's rarely a good idea to roll your own algorithm when
# there are already excellent existing packages that do it better.  In the
# case of multiple alignment, such excellent packages include ClustalW,
# HMMer, and MUSCLE.

# Constant to specify the threshold for what's "good enough"
# to be considered a match.

E_VALUE_THRESHOLD = 0.04

# Number of initial nucleotides on which to perform multiple alignment.

SEQ_LENGTH = 400

# Open up the file that contains the BLAST results and
# have BioPython organize it into a collection of records
# for us.

handle = open("struthio_camelus.xml")
records = NCBIXML.parse(handle)

# An array of sequences we'll align.

align = MultipleSeqAlignment([])

# For each result ...

for record in records:
    for alignment in record.alignments:
        for hsp in alignment.hsps:

            if 0 == len(align):

                # Add query sequence to list of sequences to align.

                record = SeqRecord(
                    Seq(hsp.query[0:SEQ_LENGTH], generic_dna),
                    id="Source")

                align.append(record)

            # Only add to the alignment if the sequence's score
            # is within the threshold.

            if hsp.expect < E_VALUE_THRESHOLD:

                # Add subject sequence to the list of sequences to align.
         
                record = SeqRecord(
                    Seq(hsp.sbjct[0:SEQ_LENGTH], generic_dna),
                    id=str(len(align)) + " " + alignment.title)

                align.append(record)

# All of the candidate sequences have been gathered.
# Align them to an output file in PHYLIP format.

AlignIO.write(align, "struthio_camelus.phy", "phylip")
