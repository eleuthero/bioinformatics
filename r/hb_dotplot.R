#!/usr/bin/R

library("seqinr")
hba <- read.fasta(file = "./hba.fasta")
hbb <- read.fasta(file = "./hbb.fasta")
hbaseq <- hba[[1]]
hbbseq <- hbb[[1]]
dotPlot(hbaseq, hbbseq)

