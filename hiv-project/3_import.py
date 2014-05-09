#!/usr/bin/python

import re
from os            import listdir
from os.path       import isfile, join
from Bio           import SeqIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

FASTA_PATH = "./sequences/"

global regex

regex  = re.compile(r"\.(\d{4})\.")

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

# Insert a sequence into the database.
def insert(consensus, record):

    year = getYear(record.description)

    print "INSERT INTO `sequence` (`year`, `consensus`, `description`, `sequence`) " \
          "VALUES ('%s', %s, '%s', '%s');" % (getYear(record.description),
                                              str(consensus),
                                              record.description,
                                              record.seq)

# =========
# Main 
# =========

# Iterate through individual sequences

for item in listdir(FASTA_PATH):
    if isfile(join(FASTA_PATH, item)) and item.endswith(".fasta"):
        for record in SeqIO.parse(join(FASTA_PATH, item), "fasta"):
            insert(False, record)

# Iterate through consensus summaries

for item in listdir(FASTA_PATH):
    if isfile(join(FASTA_PATH, item)) and item.endswith(".consensus"):
        for record in SeqIO.parse(join(FASTA_PATH, item), "fasta"):
            insert(True, record)
