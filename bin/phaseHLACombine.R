#!/usr/bin/env Rscript

##############################################
# Combine phasing results into one .txt
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

################################
# Main
################################
output = list()
for (f in files){
  pred = read.table(f, sep="\t", h=T)
  locus = gsub(".*_(.*).RData","\\1",f)
  output[[locus]] = pred
}

if(any(grep("DRB[345]", names(output)))){
  output = output[-grep("DRB[345]", names(output))]
 
}
phased = output

phased = do.call(rbind, phased)
library(ggplot2)
png("phased.png")
p = ggplot(data = phased, 
           aes(x=locus, y=phase_prob)) + geom_boxplot() + theme_bw() + xlim(c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1"))
print(p)
dev.off()

write.table(phased, output_name, quote=F, sep="\t", row.names=F)

