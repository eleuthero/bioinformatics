from Bio.Blast import NCBIWWW
from Bio import SeqIO

# Read the FASTA file we've saved.

record = SeqIO.read("opuntia_1.fasta", format="fasta")
# record = SeqIO.read("struthio_camelus_titin.fasta", format="fasta")

# Run the FASTA file through Blast as a nucleotide (NT) sequence.

handle = NCBIWWW.qblast("blastn", "nt", record.seq)

# Save the BLASTed output in its default XML format locally.

# fout = open("struthio_camelus_titin_blast.xml", "w")
fout = open("opuntia.xml", "w")

fout.write( handle.read() )
fout.close()
handle.close()
