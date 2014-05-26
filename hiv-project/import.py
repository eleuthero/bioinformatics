#!/usr/bin/python

import re
from os            import listdir
from os.path       import isfile, isdir, join
from Bio           import SeqIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

SEQUENCE_DIR = "./sequences/"

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
def insert(consensus, subtype, record):

    year = getYear(record.description)

    if consensus:
        print "INSERT INTO `sequence` (`year`, `consensus`, `threshold`, `subtype`, `description`, `sequence`) " \
              "VALUES ('%s', %s, %.2f, '%s', '%s', '%s');" % (getYear(record.description),
                                                              str(consensus),
                                                              getThreshold(record.description),
                                                              subtype,
                                                              record.description,
                                                              record.seq)
    else:
        print "INSERT INTO `sequence` (`year`, `consensus`, `threshold`, `subtype`, `description`, `sequence`) " \
              "VALUES ('%s', %s, NULL, '%s', '%s', '%s');" % (getYear(record.description),
                                                              str(consensus),
                                                              subtype,
                                                              record.description,
                                                              record.seq)

# =========
# Main 
# =========

# Iterate through individual sequences

for subtype in listdir(SEQUENCE_DIR):
    if isdir(join(SEQUENCE_DIR, subtype)):

        path = join(SEQUENCE_DIR, subtype)
        for item in sorted(listdir(path)):
            if isfile(join(path, item)) and item.endswith(".fasta"):
                for record in SeqIO.parse(join(path, item), "fasta"):
                    insert(False, subtype, record)

        # Iterate through consensus summaries

        for record in SeqIO.parse(join(path, "summary.consensus.reduced"), "fasta"):
            insert(True, subtype, record)
