#!/usr/bin/env Rscript

##############################################
# Combine imputation results into one .RData
# author:  Mareike Wendorff, Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################

################################
# Settings
################################


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

save(pred,file=output_name)
