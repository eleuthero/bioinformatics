#!/usr/bin/python

import os
import re
from os import path

SOURCES = [ { 'subtype':   'B',
              'sequences': './sequences_B.fasta',
              'patient':   './patients_B.txt' } ]
            
SEQUENCE_DIR = "./sequences"

# =========
# Functions
# =========

def getPatientInfo(path):
    with open(path) as fin:
        infos = fin.readlines()    

    l = [ ]

    for info in infos:
        tokens = info.strip().split('\t')
        if 12 == len(tokens):
            if 'HIV-1' == tokens[11]:
                i = { 'patientid' : tokens[2],
                      'accession' : tokens[3],
                      'year'      : int(tokens[7]) if tokens[7] else None } 
                l.append(i)
    return l

def getPatientInfoBySequence(ptinfos, seqinfo):
    # return [ ptinfo for ptinfo in ptinfos if ptinfo['accession'] == seqinfo['accession'] ]
    for ptinfo in ptinfos:
        if ptinfo['accession'] == seqinfo['accession']:
            return ptinfo
    return None

def getSequenceInfo(line):
    tokens = line[1:].strip().split('.')
    return { 'subtype'   : tokens[0],
             'country'   : tokens[1],
             'year'      : int(tokens[2]) if tokens[2] else None,
             'name'      : tokens[3],
             'accession' : tokens[4] }

def getSequenceFileName(subtype, year):
    global SEQUENCE_DIR

    return "%s/%s/%i.fasta" % (SEQUENCE_DIR, subtype, year)
    
def partitionSequencesByYear(source):
    global SEQUENCE_DIR

    subtype = source['subtype']

    # Ensure sequences output directory exists.

    if not os.path.exists(SEQUENCE_DIR):
        os.mkdir(SEQUENCE_DIR)
    if not os.path.exists(SEQUENCE_DIR + '/' + subtype):
        os.mkdir(SEQUENCE_DIR + '/' + subtype)

    print "Processing sequences for subtype %s." % subtype

    # Open sequence file.

    with open(source['sequences']) as fin:
        lines = fin.readlines()

    fyear = 0
    fout = None
    seqs_loaded = 0
    seqs_written = 0
    seqs_skipped = 0

    # Get patient info.

    ptinfos = getPatientInfo(source['patient'])

    # Track patients per year; we do not want to include more than one
    # sequence from any one patient for a given year.

    patients_by_year = { }
    
    # Partition sequences into individual year files for alignment.

    for line in lines:
        if line.startswith(">"):
            seqs_loaded += 1

            # Get sequence information.

            seqinfo = getSequenceInfo(line)
            year = seqinfo['year']

            # Is there a well-formed year in the sequence ?

            if not year:
                seqs_skipped += 1

            else:

                # Get patient information from this accession id.

                ptinfo = getPatientInfoBySequence(ptinfos, seqinfo)

                if not ptinfo:
                    seqs_skipped += 1

                    print "Skipping sequence with corrupted description %s." % line.strip()
                    # print "Closing file %s" % getSequenceFileName(subtype, fyear)
                    fout.close()

                    # Restart loop to keep logic simple,

                    continue
                    
                # Have we seen this patient this year ?

                if year in patients_by_year:
                    if ptinfo['patientid'] in patients_by_year[year]:
                        seqs_skipped += 1

                        # Yep.  Skip this sequence by closing the file handle.

                        if fout:
                            if not fout.closed:
                                # print "Closing file %s" % getSequenceFileName(subtype, fyear)
                                fout.close()

                        # print "Skipping patient %s for year %i." % (ptinfo['patientid'],
                        #                                             year)

                        # Restart loop to keep logic simple,

                        continue

                # Nope, this is a new year and/or a new patient.
                # Add the year and patient to our tracking array.

                if not year in patients_by_year:
                    patients_by_year[year] = [ ]
                patients_by_year[year].append(ptinfo['patientid'])

                # print "Added patient %s to year %i." % (ptinfo['patientid'],
                #                                         year)
                # print "Seen %i patients in year %i." % (len(patients_by_year[year]),
                #                                         year)

                # Do we need to close the current output file handle ?

                if (year != fyear):
                    if fout:
                        if not fout.closed:
                            # print "Closing file %s" % getSequenceFileName(subtype, fyear)
                            fout.close()

                fyear = year

                # Do we need to open the current output file handle ?

                if (not fout) or fout.closed:

                    filename = getSequenceFileName(subtype, fyear)
                    # print "Opening file %s" % filename
                    fout = open(filename, "a")

                seqs_written += 1
                fout.write(line)
        else:
            if fout:
                if not fout.closed:
                    fout.write(line)

    if fout:
        if not fout.closed:
            # print "Closing file %s" % getSequenceFileName(subtype, fyear)
            fout.close()

    print "%i of %i sequences written (%i skipped)." % (seqs_written,
                                                        seqs_loaded,
                                                        seqs_skipped)

# ===
# Main
# ====

for source in SOURCES:
    partitionSequencesByYear(source)
