#!/usr/bin/python

import re
from os            import listdir
from os.path       import isfile, isdir, join
from Bio           import SeqIO, AlignIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align     import AlignInfo
from Bio.Alphabet  import generic_dna

from lanl          import SEQUENCE_DIR

# When there are at least this many sequences for a year, generate
# 70% ~ 100% threshold dumb consensus sequences for the year.  If
# there are fewer this number, a single, simple majority consensus
# will be generated.

# =========
# Globals
# =========

THRESHOLD = 10

global regex

regex  = re.compile(r"^(\d{4})\.")
maxseq = 0

# =========
# Functions
# =========
 
# Returns the year from the sequence name.
def getYear(line):
    global regex

    match = re.search(regex, line) 

    if match:
        return match.group(1)
    else:
        return 0

def generateConsensusThreshold(subtype, year, summary, fsum, fout, threshold):
    consensus_seq = Seq(str(summary.dumb_consensus(threshold = threshold,
                                                   ambiguous = 'N',
                                                   require_multiple = 1)),
                            generic_dna)

    consensus_rec = SeqRecord(consensus_seq,
                              id="",
                              description=("%s.%s.%i%%.Consensus" % (subtype,
                                                                     year,
                                                                     100 * threshold)))
    SeqIO.write(consensus_rec, fout, "fasta")
    SeqIO.write(consensus_rec, fsum, "fasta")

def generateConsensus(subtype, year, summary, fsum, fout):
    for cvalue in range(7, 11):
        generateConsensusThreshold(subtype, year, summary, fsum, fout, float(cvalue) / 10)

def getAlignmentLength(subtype):

    # Determine length of longest aligned sequences per subtype.

    maxseq = 0

    path = join(SEQUENCE_DIR, subtype)
    for item in listdir(path):
        if isfile(join(path, item)) and item.endswith(".fasta"):
            for record in SeqIO.parse(join(path, item), "fasta"):
                if len(record.seq) > maxseq:
                    maxseq = len(record.seq)
    return maxseq

def extendAllSequences(subtype, maxseq):

    # Extend all sequences for this subtype to length of longest
    # aligned sequence by adding - to the ends of them.

    path = join(SEQUENCE_DIR, subtype)
    for item in listdir(path):
        if isfile(join(path, item)) and item.endswith(".fasta"):
            with open(join(path, item) + ".extended", "w") as fout:
               for record in SeqIO.parse(join(path, item), "fasta"):
                   record.seq += ("-" * (maxseq - len(record.seq)))
                   SeqIO.write(record, fout, "fasta")

def generateConsensusSequence(fsum, path, item):

    # Determine number of sequences for this year.

    seqcount = len( SeqIO.index(join(path, item), "fasta") )

    # Get alignment summary.  This will be required for us
    # to generate a dumb n% consensus.

    alignment = AlignIO.read(join(path, item), "fasta")
    summary = AlignInfo.SummaryInfo(alignment)

    with open(join(path, item) + ".consensus", "w") as fout:
        year = getYear(item)

        if seqcount >= THRESHOLD:
            generateConsensus(subtype, year, summary, fsum, fout)
        else:

            # Generate a majority and a 100% consensus; there aren't enough
            # sequences to justify the granularity of multiple consensus sequences.

            generateConsensusThreshold(subtype, year, summary, fsum, fout, 0.5)
            generateConsensusThreshold(subtype, year, summary, fsum, fout, 1.0)

def generateConsensusSequences(subtype):

    # Generate consensus sequences for each year and a
    # consensus summary file.

    path = join(SEQUENCE_DIR, subtype)
    with open(join(path, "summary.consensus"), "w") as fsum:
        for item in sorted(listdir(path)):
            if isfile(join(path, item)) and item.endswith(".fasta.extended"):
                generateConsensusSequence(fsum, path, item)

def generateConsensusBySubtype(subtype):

    print "Generating consensus sequences for subtype %s..." % subtype 

    maxseq = getAlignmentLength(subtype)

    print "Longest sequence length: %i." % maxseq 

    extendAllSequences(subtype, maxseq)

    print "Extended all sequences to %i residues." % maxseq

    generateConsensusSequences(subtype)

# =========
# Main
# =========

for subtype in listdir(SEQUENCE_DIR):
    if isdir(join(SEQUENCE_DIR, subtype)):
        generateConsensusBySubtype(subtype)
