#!/usr/bin/env Rscript

##############################################
# Combine phasing results
# author: Frauke Degenhardt
# contact: f.degenahrdt@ikmb.uni-kiel.de
# February 2021
##############################################

args = commandArgs(T)

tmp = read.table(args[1], h=T,sep="\t")

id = tmp$id
data = tmp[,-(1:8)]

out = apply(data, 2, function(x){a =  id[x==1]; b= id[x==2]; c = c(as.matrix(a), rep(as.matrix(b), each=2)); return(paste(c, collapse=", "))})

write.table(out,args[2], sep="\t", col.names=F, row.names=F)

