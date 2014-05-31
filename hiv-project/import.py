#!/usr/bin/python

import re
from os            import listdir
from os.path       import isfile, isdir, join
from Bio           import SeqIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from lanl          import *

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
def insertSequence(record, ptinfo):
    print "INSERT INTO `sequence` " \
          "(`year`, `consensus`, `threshold`, `subtype`, `country`, `description`, `sequence`) " \
          "VALUES ('%s', False, NULL, '%s', '%s', '%s', '%s');" % (str(ptinfo['year']),
                                                                   ptinfo['subtype'],
                                                                   ptinfo['country'],
                                                                   record.description,
                                                                   record.seq)

# Insert a consensus sequence into the database.
def insertConsensus(subtype, record):

    year = getYear(record.description)

    print "INSERT INTO `sequence` " \
          "(`year`, `consensus`, `threshold`, `subtype`, `country`, `description`, `sequence`) " \
          "VALUES ('%s', True, %.2f, '%s', NULL, '%s', '%s');" % (getYear(record.description),
                                                                  getThreshold(record.description),
                                                                  subtype,
                                                                  record.description,
                                                                  record.seq)

def importSource(source):
    global SEQUENCE_DIR

    subtype = source['subtype']
    ptinfos = getPatientInfo(source['patient'])

    if isdir(join(SEQUENCE_DIR, subtype)):
        path = join(SEQUENCE_DIR, subtype)

        # Iterate through sequences

        for item in sorted(listdir(path)):
            if isfile(join(path, item)) and item.endswith(".fasta"):
                for record in SeqIO.parse(join(path, item), "fasta"):
                    seqinfo = getSequenceInfo(record.description)
                    ptinfo = getPatientInfoBySequence(ptinfos, seqinfo)
                    insertSequence(record, ptinfo)

        # Iterate through consensus summaries

        for record in SeqIO.parse(join(path, "summary.consensus.reduced"), "fasta"):
            insertConsensus(subtype, record)

# =========
# Main 
# =========

for source in SOURCES:
    importSource(source)
