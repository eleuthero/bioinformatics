#!/usr/bin/python

import re
from os            import listdir
from os.path       import isfile, join
from Bio           import SeqIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

FASTA_PATH = "./sequences/"

global yregex
global tregex

yregex  = re.compile(r"\.(\d{4})\.")
tregex  = re.compile(r"\.(\d{2,3})%\.")

# =========
# Functions 
# =========

# Returns the year from the sequence name.
def getYear(line):
    global yregex

    match = re.search(yregex, line)

    if match:
        return match.group(1)
    else:
        return 0

# Returns the consensus threshold from the sequence name.
def getThreshold(line):
    global tregex

    match = re.search(tregex, line)

    if match:
        return float(match.group(1)) / 100
    else:
        return 0
    
# Insert a sequence into the database.
def insert(consensus, record):

    year = getYear(record.description)

    if consensus:
        print "INSERT INTO `sequence` (`year`, `consensus`, `threshold`, `description`, `sequence`) " \
              "VALUES ('%s', %s, %.2f, '%s', '%s');" % (getYear(record.description),
                                                        str(consensus),
                                                        getThreshold(record.description),
                                                        record.description,
                                                        record.seq)
    else:
        print "INSERT INTO `sequence` (`year`, `consensus`, `threshold`, `description`, `sequence`) " \
              "VALUES ('%s', %s, NULL, '%s', '%s');" % (getYear(record.description),
                                                        str(consensus),
                                                        record.description,
                                                        record.seq)

# =========
# Main 
# =========

# Iterate through individual sequences

for item in sorted(listdir(FASTA_PATH)):
    if isfile(join(FASTA_PATH, item)) and item.endswith(".fasta"):
        for record in SeqIO.parse(join(FASTA_PATH, item), "fasta"):
            insert(False, record)

# Iterate through consensus summaries

for record in SeqIO.parse(join(FASTA_PATH, "summary.consensus.reduced"), "fasta"):
    insert(True, record)
