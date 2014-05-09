#!/usr/bin/R

library("seqinr")
a <- read.fasta(file = "./hba_human.fasta")
b <- read.fasta(file = "./hba_mouse.fasta")
aseq <- translate(seq = a[[1]])
bseq <- translate(seq = b[[1]])
dotPlot(aseq, bseq, xlab="Human (Homo sapiens) hemoglobin alpha proteins", ylab="Mouse (mus musculus) hemoglobin alpha proteins")

