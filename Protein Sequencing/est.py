from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import string
import re
for seq_record in SeqIO.parse("hiv.fasta", "fasta", IUPAC.unambiguous_dna):
        coding_dna = Seq(str(seq_record.seq), IUPAC.unambiguous_dna)
        template_dna= coding_dna.reverse_complement()
        messenger_rna = coding_dna.transcribe()
        
        mRNA_Split = re.sub('-+', '*', str(messenger_rna))
        final_split = mRNA_Split.split('*')
        
        fo = open("final_rna.txt", "wb")
        fo.write(str(len(coding_dna)))
        fo.close()

##        for subsequence in final_split:
##            #    print subsequence
##                fo = open("text.txt", "wb")
##                fo.write(subsequence)
##                fo.close()
##                mRNA = Seq(subsequence, IUPAC.unambiguous_rna)
##                print(mRNA.translate())
