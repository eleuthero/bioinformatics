#!/usr/bin/python

# Remove all positions from alignment for which fewer than
# some threshold number of sequences contribute bases.

import sys
from os            import listdir
from os.path       import isfile, isdir, join
from Bio           import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq       import MutableSeq

SEQUENCE_DIR = "./sequences/"
THRESHOLD = int(sys.argv[1]) if (len(sys.argv) >= 2) else 1

# =========
# Functions
# =========

def reduceToFile(path, item):
    with open(join(path, item) + ".reduced", "w") as fout:
        for record in SeqIO.parse(join(path, item), "fasta"):

            # Create mutable sequence from record.

            mseq = record.seq.tomutable()

            # Remove each xrange from sequence.

            for xr in drop:
                mseq[min(xr) : max(xr) + 1] = ''

            # Update record sequence (preserves id, description, etc.)

            record.seq = mseq.toseq()

            # Write to .reduced file.

            SeqIO.write(record, fout, "fasta") 


# ====
# Main
# ====

for dir in listdir(SEQUENCE_DIR):
    if isdir(join(SEQUENCE_DIR, dir)):

        print "Reducing sequences for subtype %s..." % dir
        print "Dropping all alignment positions with fewer than %i contributors..." % (THRESHOLD)

        path = join(SEQUENCE_DIR, dir)
        records = list();

        # Read all sequence records into a single list.

        for item in sorted(listdir(path)):
            if isfile(join(path, item)) and item.endswith(".fasta.extended.consensus"):
                records.extend(SeqIO.parse(join(path, item), "fasta"))

        # Find maximum sequence length.

        maxlen = max( map(lambda item: len(item.seq), records) )

        # Find the number of sequences with a residue at the ith location.

        population = [0 for i in range(maxlen)]

        for record in records:
            for i in range(len(record.seq)):
                if record.seq[i] != 'N':
                    population[i] += 1

        # Determine alignment ranges to drop based on the supplied threshold.

        drop = list()
        end = None
        r = range(len(population))
        r.reverse();

        for i in r:
            if population[i] < THRESHOLD:
                if not end:
                    end = i + 1
            elif population[i] >= THRESHOLD:
                if end:
                    drop.append(xrange(i, end))
                    end = None
        if end:
            drop.append(xrange(0, end))

        # Determine length of dropped positions

        dlen = 0
        for xr in drop:
            dlen += (max(xr) - min(xr) + 1)

        print "Dropping %i of %i alignment positions." % (dlen, maxlen)
        print "New global alignment size: %i." % (maxlen - dlen)

        # Drop alignment ranges in extended FASTA files.

        for item in sorted(listdir(path)):
            if isfile(join(path, item)) and item.endswith(".extended"):
                reduceToFile(path, item)
            elif isfile(join(path, item)) and item.endswith(".consensus"):
                reduceToFile(path, item)

        # Drop alignment ranges in extended FASTA files.

        for item in sorted(listdir(path)):
            if isfile(join(path, item)) and item.endswith(".fasta.extended"):
                with open(join(path, item) + ".reduced", "w") as fout:
                    for record in SeqIO.parse(join(path, item), "fasta"):

                        # Create mutable sequence from record.

                        mseq = record.seq.tomutable()

                        # Remove each xrange from sequence.

                        for xr in drop:
                            mseq[min(xr) : max(xr) + 1] = ''

                        # Update record sequence (preserves id, description, etc.)

                        record.seq = mseq.toseq()

                        # Write to .reduced file.

                        SeqIO.write(record, fout, "fasta") 
