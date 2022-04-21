#!/usr/bin/env Rscript

##############################################
# Combine imputation results (CSV)
# author:  Mareike Wendorff
# contact: m.wendorff@ikmb.uni-kiel.de
# March 2020
##############################################

args = commandArgs(T)
files  = args[-length(args)]
output_name = args[length(args)]

files = sort(files)
output = list()
for (f in files){
  load(f)
  locus = gsub(".*_(.*).RData","\\1",f)
  output[[locus]] = pred
}
pred = output

result = do.call(cbind,sapply(pred,"[",4))[,sort(c(1,(1:length(pred))*5-3,(1:length(pred))*5-2))]
colnames(result) = c("sample.id", paste(rep(names(pred),each=2),1:2, sep="."))
write.table(result,paste("imputation_", output_name,".csv",sep=""),sep="\t", row.names=F,quote=F)
