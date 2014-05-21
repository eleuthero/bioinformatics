#!/usr/bin/python

import re
from os            import listdir
from os.path       import isfile, join
from Bio           import SeqIO, AlignIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align     import AlignInfo
from Bio.Alphabet  import generic_dna

FASTA_PATH = "./sequences/"

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

def generateConsensusThreshold(summary, fsum, fout, threshold):

    consensus_seq = Seq(str(summary.dumb_consensus(threshold = threshold,
                                                   ambiguous = 'N',
                                                   require_multiple = 1)),
                            generic_dna)

    consensus_rec = SeqRecord(consensus_seq,
                              id="",
                              description=(".%s.%i%%.Consensus" % (year, 100 * threshold)))

    SeqIO.write(consensus_rec, fout, "fasta")
    SeqIO.write(consensus_rec, fsum, "fasta")

def generateConsensus(summary, fsum, fout):
    for cvalue in range(7, 11):
        generateConsensusThreshold(summary, fsum, fout, float(cvalue) / 10)

# =========
# Main
# =========

# Determine length of longest aligned sequence

for item in listdir(FASTA_PATH):
    if isfile(join(FASTA_PATH, item)) and item.endswith(".fasta"):
        for record in SeqIO.parse(join(FASTA_PATH, item), "fasta"):
            # print "%s - %s" % (record.id, len(record.seq))
            if len(record.seq) > maxseq:
                maxseq = len(record.seq)

print "Longest sequence length: %i" % maxseq 

# Extend all sequences to length of longest aligned sequence
# by adding - to the ends of them.

for item in listdir(FASTA_PATH):
    if isfile(join(FASTA_PATH, item)) and item.endswith(".fasta"):

        fout = open(join(FASTA_PATH, item) + ".extended", "w")

        for record in SeqIO.parse(join(FASTA_PATH, item), "fasta"):
            record.seq += ("-" * (maxseq - len(record.seq)))
            SeqIO.write(record, fout, "fasta")

        fout.close()
            
print "Extended all sequences to %i residues." % maxseq

# Generate consensus sequences for each year and a
# consensus summary file.

fsum = open(join(FASTA_PATH, "summary.consensus"), "w")

for item in sorted(listdir(FASTA_PATH)):
    if isfile(join(FASTA_PATH, item)) and item.endswith(".fasta.extended"):
        alignment = AlignIO.read(join(FASTA_PATH, item), "fasta")
        summary = AlignInfo.SummaryInfo(alignment)

        fout = open(join(FASTA_PATH, item) + ".consensus", "w")
        year = getYear(item)

        generateConsensus(summary, fsum, fout)

        fout.close()

fsum.close()
