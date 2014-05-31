#!/usr/bin/python

# Library of shared functions and globals.

# =========
# Globals
# =========

SEQUENCE_DIR = "./sequences"

SOURCES = [ { 'subtype':   'B',
              'sequences': './sequences_B.fasta',
              'patient':   './patients_B.txt' },
            { 'subtype':   'C',
              'sequences': './sequences_C.fasta',
              'patient':   './patients_C.txt' } ]

# =========
# Functions
# =========

def getPatientInfo(path):
    with open(path) as fin:
        infos = fin.readlines()

    l = list()

    for info in infos:
        tokens = info.strip().split('\t')

        if 'HIV-1' == tokens[-1]:
            if 13 == len(tokens):

                # Includes Georegion field as tokens[9].

                i = { 'patientid' : tokens[2],
                      'accession' : tokens[3],
                      'subtype'   : tokens[5],
                      'country'   : tokens[6],
                      'year'      : int(tokens[7]) if tokens[7] else None }
                l.append(i)

            elif 12 == len(tokens):

                # Does not include Georegion field.

                i = { 'patientid' : tokens[2],
                      'accession' : tokens[3],
                      'subtype'   : tokens[5],
                      'country'   : tokens[6],
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

    line = line.strip()

    # Strip off a leading > if present.

    if line.startswith('>'):
        line = line[1:]

    # The sequence description contains lots of useful information joined with periods.

    tokens = line.split('.')

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
