#!/usr/bin/env Rscript
##############################################
# Combine phasing results (RDATA)
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# April 2020
##############################################


################################
# SETTINGS
################################
args = commandArgs(trailingOnly=T)
library(reshape2)
library(parallel)
library(data.table)

################################
# FUNCTIONS
################################

# convert to alleles to dosage
to.dosage = function(x,ind){
  if(length(ind)>1){
    x=x[,ind]
    x=apply(x,1,function(x){paste(x, collapse=" ")})
    # x[grep("NA",x)]=NA
    dose = as.data.frame.matrix(table(x,names(x)))
    dose[,colSums(dose)==0] = NA
    return(dose)
  }else{
    x=x[,ind]
    # x[grep("NA",x)]=NA
    dose = as.data.frame.matrix(table(x,names(x)))
    dose[,colSums(dose)==0] = NA
    return(dose)
  }
  
}

# prepare phased data for output
make_output = function(phase, locus, fam){
  # a) Get only individuals with phasing certainty > 0.8
  # b) Assign correct allele names
  locus = locus[locus%in%phase$loc] # get locus
  phase = phase[phase$loc%in%locus,] # get locus
  phase$loc = factor(phase$loc, levels=locus)
  id=paste(phase$id, phase$loc)[phase$phase_prob <0.8] # a)
  phase[paste(phase$id, phase$loc)%in%id ,c("X1","X2")]=NA
  phase = data.frame(phase)
  phase[c(grep("/", phase$X1), grep("/", phase$X2)),c("X1","X2")] = NA
  phase$chr1 = paste(phase$locus, phase$X1,sep="*") # b)
  phase$chr2 = paste(phase$locus, phase$X2,sep="*") # b)
  
  phase=phase[phase$id%in%fam$IID,]
  
  
  phase = phase[, c("id","locus", "chr1","chr2")]
  
  print(head(phase))
  data = rbind(do.call(rbind,tapply(phase$chr1, phase$id, invisible)),
               do.call(rbind,tapply(phase$chr2, phase$id, invisible)))
  
  comb=list(1:length(locus))
  print(comb)
  
  
  # Make dosage
  dose = mclapply(comb, function(x){
    x=as.numeric(unlist(x))
    print(x)
    return(to.dosage(data,x))
  },  mc.preschedule=F, mc.cores = detectCores())
  
  dose = do.call(rbind, dose)
  info = cbind(6,rownames(dose),"NA","A","P")
  
  colnames(info)[1:5]= c("chr","id","pos","REF","ALT")
  
  print(fam$IID)
  print(dim(dose))
  print(table(fam$IID%in%colnames(dose))) 
  print(head(colnames(dose)))
  dose = dose[,match(fam$IID, colnames(dose))]
  
  data= cbind(info, dose)
  data = data.frame(data)
  data$id=gsub("[0-9].*\\.","", data$id)
  id = (strsplit(as.matrix(data$id), " "))
  id = lapply(id, function(x){names(x) = gsub("\\*.*", "", x); x= x[match(locus, names(x))]; x=paste(x, collapse="-"); return(x)})
  data$id = unlist(id)
  return(data)
  
}
################################
# MAIN 
################################

# READ FILES
phase = read.table(paste0("imputation_",args[1],".META.PHASING.txt"), sep="\t", h=T)
fam = read.table(paste0(args[1],".fam"), h=F, col.names=c("FID","IID","","","",""))

# PREPARE OUTPUT PHASING (ACB, DP, DQ, DR-DQ-DP or whole haplotype)
list = c()
loci = unique(phase$locus)
if("C"%in%loci & "B"%in%loci & "A"%in%loci){
  out = make_output(phase, c("A","C","B"),fam)
  list[["ACB"]] = out
}

if("DPA1"%in%loci & "DPB1"%in%loci){
  out = make_output(phase, c("DPA1","DPB1"),fam)
  list[["DP"]] = out
}

if("DQA1"%in%loci & "DQB1"%in%loci){
  out = make_output(phase, c("DQA1","DQB1"),fam)
  list[["DQ"]] = out
}

if("DQA1"%in%loci & "DQB1"%in%loci & "DRB1"%in%loci){
  out = make_output(phase, c("DRB1", "DQA1","DQB1"),fam)
  list[["DRDP"]] = out
}

if("DQA1"%in%loci & "DQB1"%in%loci & "DRB1"%in%loci & "DPA1"%in%loci & "DPB1"%in%loci){
  out = make_output(phase, c("DRB1", "DQA1","DQB1","DPA1","DPB1"),fam)
  list[["DRDQDP"]] = out
}

loci = factor(loci, levels=c("A","E","C","B","DRB1","DQA1","DQB1","DPA1","DPB1"))
loci=sort(loci)
out = make_output(phase, loci,fam)

list[["all"]]= out

data = do.call(rbind, list)
if(any(-grep("NA", rownames(data)))){
  data = data[-grep("NA", rownames(data)),]
}

save(data, file=paste0("imputation_",args[1],".haplotypes.RData", sep=""))

