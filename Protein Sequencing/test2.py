from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
#for seq_record in SeqIO.parse("hiv.fasta", "fasta", IUPAC.unambiguous_dna):
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
template_dna= coding_dna.reverse_complement()
messenger_rna = coding_dna.transcribe()
print(messenger_rna)
print(messenger_rna.translate())
#print("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG".translate())


