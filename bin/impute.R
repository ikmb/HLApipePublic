#!/usr/bin/env Rscript

##############################################
# impute HLA alleles with HIBAG
# author:  Mareike Wendorff
# contact: m.wendorff@ikmb.uni-kiel.de
# March 2020
##############################################

################################
# Settings
################################

library(HIBAG)
args = commandArgs(T)

bed  = paste(args[1],".bed",sep="")
bim  = paste(args[1],".bim",sep="")
fam  = paste(args[1],".fam",sep="")

params.assembly = args[2]
params.model = args[3]
gen = args[4]


################################
# Functions
################################

loadObject = function(x){
  model = local({
    load(x)
    stopifnot(length(ls())==1)
    environment()[[ls()]]
    
  } 
  )
  return(model)}



################################
# Main
################################

## read PLINK data
geno_data <- hlaBED2Geno(bed.fn=bed, fam.fn=fam, bim.fn=bim,assembly=params.assembly)
summary(geno_data)


## import model
HLAModelList = loadObject(params.model)

## prediction
  hla.id  = which(names(HLAModelList) %in% gen) # format loci as "A,B,DRB1"
  model <- hlaModelFromObj(HLAModelList[[hla.id]])
  # HLA allele frequencies
  #	cbind(frequency = model$hla.freq)
  #	frequency
  
  
  # best-guess genotypes and all posterior probabilities
  pred <- predict(model, geno_data, type="response+prob", match.type="Position")
  summary(pred)
  pred$value
  pred$postprob

save(pred, file=paste0("imputation_",args[1],"_",gen,".RData"))
