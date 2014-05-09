#!/usr/bin/R

library("seqinr")
a <- read.fasta(file = "./hba_human.fasta")
b <- read.fasta(file = "./hba_mouse.fasta")
aseq <- a[[1]]
bseq <- b[[1]]
dotPlot(aseq, bseq, xlab="Human (Homo sapiens) hemoglobin alpha", ylab="Mouse (mus musculus) hemoglobin alpha")

