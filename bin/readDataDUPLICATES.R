#!/usr/bin/env Rscript

##############################################
# Calcualte marginal posterior probability of alleles
# author:  Mareike Wendorff
# contact: m.wendorff@ikmb.uni-kiel.de
# March 2020
##############################################
args = commandArgs(T)

bim  = paste(args[1],".bim",sep="")

snp_file = read.csv(bim,stringsAsFactors = F,sep="",header = F)
exclude = snp_file$V2[duplicated(paste(snp_file$V1,snp_file$V4))|duplicated(paste(snp_file$V1,snp_file$V4),fromLast = T)]
write.table(exclude,file = "exclude_duplicates.txt",row.names = F,col.names = F,quote = F)

#check if duplicates / flipped
#multiallelic delete completely (check how many)
#plink -- --missing