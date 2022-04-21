#!/usr/bin/env Rscript

##############################################
# Combine phasing results (CSV)
# author: Frauke Degenhardt
# contact: f.degenahrdt@ikmb.uni-kiel.de
# February 2021
##############################################

args = commandArgs(T)

load(args[1])

tmp = data
id = tmp$id
data = tmp[,-(1:8)]

out = apply(data, 2, function(x){a =  id[x==1]; b= id[x==2]; c = c(as.matrix(a), rep(as.matrix(b), each=2)); return(paste(c, collapse=", "))})
out = cbind(names(out), out)
colnames(out)=c("IID", "haplotypes")
write.table(out,args[2], sep="\t", col.names=F, row.names=F)

