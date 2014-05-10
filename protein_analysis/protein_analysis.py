#!/usr/bin/python

from os           import listdir
from os.path      import isfile, join
from Bio          import SeqIO
from Bio.SeqUtils import ProtParam

PROTEINFILE = "./proteins.fasta"

# ====
# Main
# ====

print "Description\tpI\tMol. Weight\tAromaticity\tInstability Index"

# Perform a protein analysis on every protein sequence in the
# specified FASTA file. 

for record in SeqIO.parse(PROTEINFILE, "fasta"):
    ps = ProtParam.ProteinAnalysis(str(record.seq))

    print "%s\t%.5f\t%.5f\t%.5f\t%.5f" % (record.description[:40],
                                          ps.isoelectric_point(),
                                          ps.molecular_weight(),
                                          ps.aromaticity(),
                                          ps.instability_index())
