#!/usr/bin/env Rscript

##############################################
# Combine phasing results (PLINK)
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# April 2020
##############################################

################################
# Settings
################################
args = commandArgs(trailingOnly=T)

library("reshape2")
library("parallel")

print(args)
################################
# FUNCTIONS
################################
format_to_plink = function(fam, data){
  tr = apply(data, 1, function(x){
     translate(x[4], x[5], x[-(1:5)])
  })
  print(head(tr))
  ped = cbind(fam, tr[match(rownames(tr),fam$IID),])
  map = cbind(data[,c(1,2)], 0, data[,3])
  map[,2] = gsub("\\*|:","_", map[,2])
  return(list(ped=ped, map=map))
}

################################
# MAIN
################################

load(paste0("imputation_", args[1],".haplotypes.RData"))
fam= read.table(paste0(args[1],".fam"), h=F, col.names = c("FID","IID","","","","PHENO"))
# CONVERT TO PLINK OUTOUT
source(args[2])
data[,3]=1:nrow(data)

plink = format_to_plink(fam, data)

write.table(plink$map, paste("imputation_",args[1],".haplotypes.map", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
write.table(plink$ped, paste("imputation_",args[1],".haplotypes.ped", sep=""), quote=F, sep="\t", row.names=F, col.names=F)

cmd = paste("plink --file", paste0("imputation_",args[1],".haplotypes"),
            "--make-bed --out", paste0("imputation_",args[1],".haplotypes"))
system(cmd)
