This is a game-like simulation of the application of set covering to SNPs.
The point of the exercise is to find a (minimal) tag set for the given set of SNPs and haplotypes.

Click the checkbox to the right of the snp to add it to the tag set.  When few checkboxes are checked, nothing will appear to happen, but as you select more checkboxes, the snps you have selected will be sufficient to distinguish at least a couple of haplotypes.  When this happens, the haplotype identifier at the top of the page will light up in green.  When all of the haplotype identifiers at the top of the page are lit in green, it means that you have generated a tag set of snps that are capable of distinguishing all haplotypes.

Click the "Find the next cover candidate SNP" to invoke a greedy algorithm to make what is likely, but not guaranteed, to be the most effective guess on the next snp to generate a tag set (a minimal set of snps required to distinguish all haplotypes in the sample).
