#!/usr/bin/python

from os      import listdir
from os.path import isfile, isdir, join
from Bio     import SeqIO

SEQUENCE_DIR = "./sequences/"
ROWLEN = 50

# =========
# Functions
# =========

def createAlignmentHeatMap(fout, population, maxlen):
    fout.write(
        """
        <html>
          <head>
            <style>
              td.A     { width: 12px; height: 12px; font-size: 8px; color: black;
                         vertical-align: center; text-align: center; }
              td.A.OK  { background-color: green; }
              td.A._10 { background-color: yellow; }
              td.A._5  { font-size: 10px; background-color: orange; }
              td.A._4  { font-size: 10px; background-color: black; color: orange; }
              td.A._3  { font-size: 10px; background-color: black; color: orange; }
              td.A._2  { font-size: 10px; background-color: black; color: red; }
              td.A._1  { font-size: 10px; background-color: black; color: red; }
              td.A._0  { font-size: 10px; background-color: black; color: red; }
            </style>
          </head>
          <body>
            <h1>Alignment Map</h1>
            <h6>Numbers indicate number of sequences that have residues at the specified alignment offset.</h6>
            <table>""")

    for i in range((maxlen / ROWLEN) + 1):
        fout.write("<tr>")
        fout.write("<td>%i</td>" % (i * ROWLEN + 1))
        for j in range(ROWLEN):
            k = i * ROWLEN + j
            if k < maxlen:
                if population[k] > 99:
                    fout.write("<td class='A OK' title='position: %i, population: %i'>OK</td>" % (k, population[k]))
                else:
                    css = None
                    if population[k] >= 25:
                        css = "OK"
                    elif population[k] >= 10:
                        css = "_10"
                    elif population[k] >= 5:
                        css = "_5"
                    elif population[k] == 4:
                        css = "_4"
                    elif population[k] == 3:
                        css = "_3"
                    elif population[k] == 2:
                        css = "_2"
                    elif population[k] == 1:
                        css = "_1"
                    else:
                        css = "_0"

                    fout.write("<td class='A %s' title='%i'>%i</td>" % (css,
                                                                        k,
                                                                        population[k]))
            else:
                fout.write("<td />")
        fout.write("<td>%i</td>" % min([(i + 1) * ROWLEN, maxlen - 1]))
        fout.write("</tr>")

    fout.write(
        """
        </table>
     </body>
    </html>""")

# ====
# Main 
# ====

for subtype in listdir(SEQUENCE_DIR):
    if isdir(join(SEQUENCE_DIR, subtype)):

        path = join(SEQUENCE_DIR, subtype)
        records = list();

        # Read all sequence records into a single list.

        for item in sorted(listdir(path)):
            if isfile(join(path, item)) and item.endswith(".fasta.extended.reduced"):
                records.extend(SeqIO.parse(join(path, item), "fasta"))

        # Find maximum sequence length.

        maxlen = max( map(lambda item: len(item.seq), records) )

        # Find the number of sequences with a residue at the ith location.

        population = [0 for i in range(maxlen)]

        for record in records:
            for i in range(len(record.seq)):
                if record.seq[i] not in ['N', '-']:
                    population[i] += 1

        with open("alignment_%s.html" % subtype, "w") as fout:
            s = createAlignmentHeatMap(fout, population, maxlen)
