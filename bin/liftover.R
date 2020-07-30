#!/usr/bin/env Rscript

##############################################
# LiftOver
# author:  Mareike Wendorff
# contact: m.wendorff@ikmb.uni-kiel.de
# April 2020
##############################################
################################
# Settings
################################

args = commandArgs(T)
liftover.file = args[1]
unMapped.file = args[2]
prelift.files = args[3]
postlift.files = args[4]

################################
# Main
################################
bim.prelift.file = paste0(prelift.files,".bim")
bim.postlift.file = paste0(postlift.files,".bim")
#library(xtable)
#library(data.table)

liftover = read.csv(liftover.file, header = F, sep = ":")
liftover$V1 = gsub("chr","",liftover$V1)
liftover$V1 = as.numeric(sapply(liftover$V1,function(x) switch(x, X=23,Y=24,XY=25,MT=26,x)))
liftover$V2 = as.numeric(gsub("(.*)-.*","\\1",liftover$V2))

if (file.info(unMapped.file)$size != 0){
  unMapped = read.csv(unMapped.file, header = F, sep = ":")
  unMapped = unMapped[grep("chr",unMapped$V1),]
  unMapped$V1 = gsub("chr","",unMapped$V1)
  unMapped$V1 = as.numeric(sapply(unMapped$V1,function(x) switch(x, X=23,Y=24,XY=25,MT=26,x)))
  unMapped$V2 = as.numeric(gsub("(.*)-.*","\\1",unMapped$V2))
}

bim.prelift = read.csv(bim.prelift.file, header = F, sep = "\t") 

if (exists("unMapped")){
  snp.rm = bim.prelift$V2[paste(bim.prelift$V1,bim.prelift$V4) %in% paste(unMapped$V1,unMapped$V2)]
  bim.tmp = bim.prelift[!paste(bim.prelift$V1,bim.prelift$V4) %in% paste(unMapped$V1,unMapped$V2),]
} else {
  snp.rm = ""
  bim.tmp = bim.prelift
}
bim.tmp$V1 = liftover$V1
bim.tmp$V4 = liftover$V2

if (exists("unMapped")){
  unmappedCHR = as.data.frame(table(unMapped$V1))
  colnames(unmappedCHR) = c("CHR","unmappedSNPs")
  write.table(snp.rm, file = paste0("unmappedSNPs_",postlift.files,".txt"), quote = F, row.names = F, col.names = F)
  write.table(unmappedCHR, file = paste0("unmappedCHR_",postlift.files,".txt"), quote = F, row.names = F, col.names = T)
} else {
  write.table("", file = paste0("unmappedSNPs_",postlift.files,".txt"), quote = F, row.names = F, col.names = F)
  write.table("", file = paste0("unmappedCHR_",postlift.files,".txt"), quote = F, row.names = F, col.names = T)
}

cmd = paste0("plink --bfile ", prelift.files, " --exclude unmappedSNPs_",postlift.files,".txt --make-bed --out ", postlift.files)
system(cmd)

write.table(bim.tmp, file = bim.postlift.file, quote = F, row.names = F, col.names = F)

