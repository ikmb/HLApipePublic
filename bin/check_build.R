#!/usr/bin/env Rscript

##############################################
# Compare study to reference
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# April 2020
##############################################

options(stringsAsFactors = F)
################################
# Settings
################################

args = commandArgs(T)

################################
# Main
################################
# CHECK IF BUILD IS hg19

reference = read.table(paste(args[1],".bim",sep=""), h=F)
study = read.table(paste(args[2],".bim",sep=""), h=F)


tab=table(study$V4%in%reference$V4)
tab = round(tab/sum(tab),4)*100

print(tab)

if(tab["TRUE"]< 0.8){
  stop("Check the genome build of your data. It should be hg19 (GRCh37). Or use --assembly flag to update from hg18 to hg19. Other builds are currently not supported.")
}