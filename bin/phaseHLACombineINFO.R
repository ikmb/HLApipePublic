#!/usr/bin/env Rscript
##############################################
# HLA phasing INFO
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# April 2020
##############################################


################################
# SETTINGS
################################
args = commandArgs(trailingOnly=T)

source(args[3])
fam= read.table(paste0(args[1],".fam"), h=F, col.names = c("FID","IID","","","","PHENO"))
data= read.table(paste0(args[2]), h=T,sep="\t")

frq=make.frq.hwe(data[,-(1:8)], as.numeric(as.matrix(fam$PHENO)))
write.table(cbind(data[,1:8],frq), paste0("imputation_", args[1],".haplotypes.info"), quote=F, sep="\t", row.names=F)
