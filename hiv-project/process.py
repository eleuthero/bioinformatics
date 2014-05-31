#!/usr/bin/python

import os
import re
from os import path

SOURCES = [ { 'subtype':   'B',
              'sequences': './sequences_B.fasta',
              'patient':   './patients_B.txt' },
            { 'subtype':   'C',
              'sequences': './sequences_C.fasta',
              'patient':   './patients_C.txt' } ]
            
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

        if 'HIV-1' == tokens[-1]:
            if 13 == len(tokens):

                # Includes Georegion field as tokens[9].

                i = { 'patientid' : tokens[2],
                      'accession' : tokens[3],
                      'year'      : int(tokens[7]) if tokens[7] else None } 
                l.append(i)

            elif 12 == len(tokens):

                # Does not include Georegion field.

                i = { 'patientid' : tokens[2],
                      'accession' : tokens[3],
                      'year'      : int(tokens[7]) if tokens[7] else None } 
                l.append(i)

    print "Loaded %i patient information records from %s." % (len(l),
                                                              path)
    return l

def getPatientInfoBySequence(ptinfos, seqinfo):
    # return [ ptinfo for ptinfo in ptinfos if ptinfo['accession'] == seqinfo['accession'] ]
    for ptinfo in ptinfos:
        if ptinfo['accession'] == seqinfo['accession']:
            return ptinfo
    return None

def getSequenceInfo(line):

    # The sequence description contains lots of useful information joined with periods.

    tokens = line[1:].strip().split('.')

    # Check to make sure that the year token is numeric; sometimes there is no year
    # information (it's usually '-' in this case); we can't do anything without a year.

    if not tokens[2].isdigit():
        return None
 
    # However, sometimes the name of the sequence itself contains periods too.
    # We'll get around this by popping non-name tokens off of the sequence description
    # and whatever's left we'll call the sequence name.

    return { 'subtype'   : tokens.pop(0),
             'country'   : tokens.pop(0),
             'year'      : int( tokens.pop(0) ),
             'accession' : tokens.pop(-1),         # Pop accession off of the end of list.
             'name'      : '.'.join(tokens) }      # Whatever is left in the list is the name.

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

    print "Partitioning sequences for subtype %s..." % subtype

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

            # Check to see whether we successfully parsed the sequence information.

            if not seqinfo:
                seqs_skipped += 1

                print "Skipping sequence with unparseable description %s." % line.strip()

                if fout:
                    if not fout.closed:
                        # print "Closing file %s." % getSequenceFileName(subtype, fyear)
                        fout.close()

                # Restart loop to keep logic simple.

                continue
            
            year = seqinfo['year']

            # Get patient information from this accession id.

            ptinfo = getPatientInfoBySequence(ptinfos, seqinfo)

            if not ptinfo:
                seqs_skipped += 1

                print "Skipping sequence %s with no patient info." % line.strip()

                if fout:
                    if not fout.closed:
                        # print "Closing file %s." % getSequenceFileName(subtype, fyear)
                        fout.close()

                # Restart loop to keep logic simple.

                continue
                    
            # Have we seen this patient this year ?

            if year in patients_by_year:
                if ptinfo['patientid'] in patients_by_year[year]:
                    seqs_skipped += 1

                    # Yep.  Skip this sequence by closing the file handle.

                    if fout:
                        if not fout.closed:
                            # print "Closing file %s." % getSequenceFileName(subtype, fyear)
                            fout.close()

                    # print "Skipping patient %s for year %i." % (ptinfo['patientid'],
                    #                                             year)

                    # Restart loop to keep logic simple.

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
                        # print "Closing file %s." % getSequenceFileName(subtype, fyear)
                        fout.close()

            fyear = year

            # Do we need to open the current output file handle ?

            if (not fout) or fout.closed:
                filename = getSequenceFileName(subtype, fyear)
                # print "Opening file %s." % filename
                fout = open(filename, "a")

            seqs_written += 1
            fout.write(line)

        else:
            if fout:
                if not fout.closed:
                    fout.write(line)

    # All done with sequences for this subtype.
    # Is there still a file handle open ?

    if fout:
        if not fout.closed:
            # print "Closing file %s." % getSequenceFileName(subtype, fyear)
            fout.close()

    # Summary.

    print "%i of %i sequences written (%i skipped)." % (seqs_written,
                                                        seqs_loaded,
                                                        seqs_skipped)

# ===
# Main
# ====

for source in SOURCES:
    partitionSequencesByYear(source)
