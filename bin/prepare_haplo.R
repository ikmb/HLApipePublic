#!/usr/bin/env Rscript
##############################################
# HLA phasing
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# April 2020
##############################################


library(reshape2)
library(parallel)
library(data.table)

################################
# Settings
################################
args = commandArgs(trailingOnly=T)

################################
# FUNCTIONS
################################
######################################################
# Make bootstrap samples to determine 95% confidence of
# r² (LD)
####################################################
bootstrap= function(x,ind){
  x=x[,ind]
  x[grep("NA",x)]=NA
  x = na.omit(x)
  boot=list()
  boot[[1]]=apply(x,1,function(x){paste(x, collapse=" ")})
  for(i in 2:(2001)){
    bsi = sample(1:nrow(x),nrow(x), replace=T)
    boot[[i]] = apply(x[bsi,],1,function(x){paste(x, collapse=" ")})
  }
  out= lapply(boot, linkage)
  x = out[[1]]
  out = do.call(rbind,out)
  alpha=0.05
  out = tapply(out[,1], rownames(out), function(x){return(quantile(x, probs=c(1-alpha/2,alpha/2), type=6))})
  out = do.call(rbind,out)
  #  out = cbind(2*x[,1] - out[,1], 2*x[,1]-out[,2])
  out = cbind(x,out[match(rownames(x), rownames(out)),])
  return(out)
}

#####################
# Calculate linkage
####################
linkage = function(x){
  ab = table(x)/length(x)
  x = do.call(rbind,strsplit(x, " "))
  a = table(x[,1])/(nrow(x))
  b = table(x[,2])/(nrow(x))
  ind_a = match(gsub(" .*","", names(ab)), names(a))
  ind_b = match(gsub(".* ","", names(ab)), names(b))
  a= a[ind_a]
  b = b[ind_b]
  
  x = cbind((ab -a*b)^2/(a*(1-a)*b*(1-b)),a,b,ab,ab/a, ab/b)
  
  return(x)
}

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

################################
# MAIN 
################################
################################
# MAIN:READ
################################

phase = read.table(paste0(args[1],".HLA.all.phased.txt"), sep="\t", h=T)
fam = read.table(paste0(args[1],".fam"), h=F, col.names=c("FID","IID","","","",""))

head(phase)
head(fam)
################################
# MAIN:CALCULATE R^2 LIKE MEASURES 
# only for larger sample sets (>100)
# a) Use only alleles that could 
# be phased with certainty >= 0.8
# at specificed loci
# b) Reformat HLA alleles to e.g. A*01:02
################################
phase.tmp = phase
phase = phase[!phase$loc%in%c("DRB3","DRB4","DRB5"),]
id=paste(phase$id, phase$loc)[phase$phase_prob <0.8] # a)
phase[paste(phase$id, phase$loc)%in%id ,c("X1","X2")]=NA
phase = data.frame(phase)
phase[c(grep("/", phase$X1), grep("/", phase$X2)),c("X1","X2")] = NA
phase$chr1 = paste(phase$locus, phase$X1,sep="*") # b)
phase$chr2 = paste(phase$locus, phase$X2,sep="*") # b)  


phase=phase[phase$id%in%fam$IID,]

phase = phase[, c("id","locus", "chr1","chr2")]

data = rbind(do.call(rbind,tapply(phase$chr1, phase$id, invisible)),
             do.call(rbind,tapply(phase$chr2, phase$id, invisible)))

write.table(data, paste0(args[1], ".phased_overview.txt"), col.names=F, row.names=F, quote=F)
if(nrow(fam)> 100){

  
  loci=sort(unique(phase$loc))
  colnames(data) = loci
  poss = c("A", "E", "C",  "B",  "DRB1",  "DQA1",  "DQB1",  "DPA1",  "DPB1") # possible LOCI
  poss = poss[poss%in%loci]
  data = data[,poss] # sort A C B DRB1 DQA1 DQB1 DPA1 DPB1
  if("DQA1"%in%poss & "DQB1"%in%poss){
    data = cbind(data,apply(data[,c("DQA1", "DQB1")],1,function(x){paste(x, collapse="-")}))
  }
  if("DPA1"%in%poss & "DPB1"%in%poss){
    data = cbind(data,apply(data[,c("DPA1", "DPB1")],1,function(x){paste(x, collapse="-")}))
  }
  ############################################
  # Make combination of all possible alleles
  # i.e. AB AC ADRB1 ADQA1-DQB1 etc.. to
  # calculate linkage between them
  ###########################################
  

  comb = as.data.frame.matrix(expand.grid(1:ncol(data),1:ncol(data)))
  out = mclapply(as.list(rownames(comb)), function(x){
    print(x)
    res = bootstrap(data,as.numeric(comb[x,]))
    return(res)
  }, mc.preschedule=F, mc.cores = detectCores())
  
  res = data.frame(do.call(rbind,out))
  
  colnames(res)=c("r.sq","a","b", "ab", "ab/a", "ab/b", "2.5","97.5")
  res$A=gsub(" .*","", rownames(res))
  res$B = gsub(".* ","", rownames(res))
  res = res[ res$A !=res$B,]
  
  write.table(res, paste0(args[1], ".whole.linkage"), quote=F, row.names=T, col.names=T, sep="\t")
  
}




################################
# MAIN: PREPARE FOR ASSOC
# a) Use only alleles that could 
# be phased with certainty >= 0.8
# at specificed loci
# b) Reformat HLA alleles to e.g. A*01:02
################################

phase = phase.tmp

make_output = function(phase, locus){
  locus = locus[locus%in%phase$loc]
  phase = phase[phase$loc%in%locus,]
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

  
  
  dose = mclapply(comb, function(x){
    x=as.numeric(unlist(x))
    print(x)
    return(to.dosage(data,x))
  },  mc.preschedule=F, mc.cores = detectCores())
  
  dose = do.call(rbind, dose)
  info = cbind(6,rownames(dose),"NA","A","P","NA",0,"HIBAG_HLA_hapl")
  
  colnames(info)[1:8]= c("chr","id","pos","REF","ALT","INFO","TYPE","source")
  
  
  dose = dose[,match(fam$IID, colnames(dose))]
  
  data= cbind(info, dose)
  data = data.frame(data)
  data$id=gsub("[0-9].*\\.","", data$id)
  id = (strsplit(as.matrix(data$id), " "))
  id = lapply(id, function(x){names(x) = gsub("\\*.*", "", x); x= x[match(locus, names(x))]; x=paste(x, collapse=" "); return(x)})
  data$id = unlist(id)
  return(data)
  
}
list = c()
loci = unique(phase$locus)
if("C"%in%loci & "B"%in%loci & "A"%in%loci){
  out = make_output(phase, c("A","C","B"))
  list[["ACB"]] = out
}

if("DPA1"%in%loci & "DPB1"%in%loci){
  out = make_output(phase, c("DPA1","DPB1"))
  list[["DP"]] = out
}

if("DQA1"%in%loci & "DQB1"%in%loci){
  out = make_output(phase, c("DQA1","DQB1"))
  list[["DQ"]] = out
}

if("DQA1"%in%loci & "DQB1"%in%loci & "DRB1"%in%loci){
  out = make_output(phase, c("DRB1", "DQA1","DQB1"))
  list[["DRDP"]] = out
}

if("DQA1"%in%loci & "DQB1"%in%loci & "DRB1"%in%loci & "DPA1"%in%loci & "DPB1"%in%loci){
  out = make_output(phase, c("DRB1", "DQA1","DQB1","DPA1","DPB1"))
  list[["DRDQDP"]] = out
}

loci = factor(loci, levels=c("A","E","C","B","DRB1","DQA1","DQB1","DPA1","DPB1"))
loci=sort(loci)
out = make_output(phase, loci)

list[["all"]]= out

data = do.call(rbind, list)
if(any(-grep("NA", rownames(data)))){
  data = data[-grep("NA", rownames(data)),]
}

write.table(data, paste0("imputation_",args[1],".hapl.data", sep=""), quote=F, sep="\t", row.names=F)
