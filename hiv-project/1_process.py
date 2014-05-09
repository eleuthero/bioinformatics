#!/usr/bin/python

import re

with open("hiv-db-ALL.fasta") as fin:
    lines = fin.readlines()

regex = re.compile(r"\.(\d{4})\.")
currentyear = 0
fout = None

# Partition sequences into individual year files for alignment

for line in lines:
    if line.startswith(">"):
        match = re.search(regex, line)

        if match:
            if fout:
                if not fout.closed:
                    print "Closing file ", currentyear, ".fasta"
                    fout.close()

            currentyear = match.group(1)

            filename = "./sequences/" + currentyear + ".fasta"

            print "Opening file ", filename
            fout = open(filename, "a")
            fout.write(line)

    else:
        fout.write(line)

fout.close()
