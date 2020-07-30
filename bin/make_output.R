#!/usr/bin/env Rscript

##############################################
# Combine imputation results from .RData object
# author:  Mareike Wendorff
# contact: m.wendorff@ikmb.uni-kiel.de
# March 2020
##############################################

args = commandArgs(T)


load(args[1])
result = do.call(cbind,sapply(pred,"[",4))[,sort(c(1,(1:length(pred))*5-3,(1:length(pred))*5-2))]
colnames(result) = c("sample.id", paste(rep(names(pred),each=2),1:2, sep="."))
write.table(result,paste("imputation_", args[2],".csv",sep=""),row.names=F,quote=F)
