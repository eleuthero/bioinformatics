from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Constant to specify the threshold for what's "good enough".

E_VALUE_THRESHOLD = 0.04

# handle = open("opuntia.xml")
handle = open("struthio_camelus.xml")

# Read the entire file we just opened and organize it for us.

records = NCBIXML.parse(handle)

# Now we can, in a completely hierarchical and programmatic
# manner, just stroll through all of the information.

for record in records:
    for alignment in record.alignments:
        for hsp in alignment.hsps:

            # Only add to the alignment if the sequence's score
            # is within the threshold.

            if hsp.expect < E_VALUE_THRESHOLD:

                # Print out some information about the alignment.

                print "sequence: ", alignment.title
                print "score   : ", hsp.score
                print "length  : ", alignment.length
                print "e value : ", hsp.expect
                print hsp.query[0:75] + "..."
                print hsp.match[0:75] + "..."
                print hsp.sbjct[0:75] + "..."
                print
