#!/usr/bin/env Rscript
##############################################
# Combine imputation results (PLINK)
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
  ped = cbind(fam, tr[match(rownames(tr),fam$IID),])
  map = cbind(data[,c(1,2)], 0, data[,3])
  map[,2] = gsub("\\*|:","_", map[,2])
  return(list(ped=ped, map=map))
}

################################
# MAIN
################################

load(paste0("imputation_", args[1],".RData"))
source(args[2])
fam= read.table(paste0(args[1],".fam"), h=F, col.names = c("FID","IID","","","","PHENO"))
plink = format_to_plink(fam, data)

write.table(plink$map, paste("imputation_",args[1],".map", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
write.table(plink$ped, paste("imputation_",args[1],".ped", sep=""), quote=F, sep="\t", row.names=F, col.names=F)

cmd = paste("plink --file", paste0("imputation_",args[1]),
            "--make-bed --out", paste0("imputation_",args[1]))
system(cmd)
