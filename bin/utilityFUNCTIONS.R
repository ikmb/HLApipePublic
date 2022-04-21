#!/usr/bin/env Rscript

##############################################
# Prepare AF output
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# April 2020
##############################################

args = commandArgs(T)
################################
# Settings
################################
library(data.table)
################################
# Functions
################################
hwe.pval = function(x){
  x= factor(x, levels=c(0,1,2))
  x = table(x)
  
  p = (2*x[1]+ x[2])/(2*sum(x))
  q = (2*x[3]+ x[2])/(2*sum(x))

  exp = c(p^2*sum(x), 2*p*q*sum(x), q^2*sum(x))

  if(x[1]==sum(x) | x[3]==sum(x)){
    chi=0
  }else{                                                                                                                    
    chi = sum((x -exp)^2/exp)
  }                                                                                                                     
  p_chi = 1- pchisq(chi,1)
  return(c(P=p_chi))


}

make.frq.hwe =function(x, case_control, hwe.bol=T){

  case_control= as.numeric(case_control)
  #1=control, 2=case
  print(dim(x))
  x=t(x)
  print(dim(x))
  x= apply(x,1,as.numeric)
  print(dim(x))

  frq = data.frame(AF_ALL=apply(x,1,function(x){x=na.omit(x); return(round(sum(x)/(2*length(x)),6))}))
  frq$AF_CASE=NA
  frq$AF_CONTROL=NA
  frq$AF_UKN = NA
  no = data.frame(SAMPLES_ALL=apply(x,1,function(x){x=na.omit(x);return(length(x))}))
  no$SAMPLES_CASE=0
  no$SAMPLES_CONTROL=0
  no$SAMPLES_UKN=0
if(any(case_control==2) & sum(case_control==2)>1){
  frq$AF_CASE=apply(x[,case_control==2],1,function(x){x=na.omit(x); return(round(sum(x)/(2*length(x)),6))})
  no$SAMPLES_CASE = apply(x[,case_control==2],1,function(x){x=na.omit(x); return(length(x))})
}
 
if(any(case_control==1) & sum(case_control==1)>1){
  frq$AF_CONTROL=apply(x[,case_control==1],1,function(x){x=na.omit(x); return(round(sum(x)/(2*length(x)),6))})
  no$SAMPLES_CONTROL = apply(x[,case_control==1],1,function(x){x=na.omit(x); return(length(x))})
 }


if(any(case_control!=1 & case_control!=2) & (sum(case_control!=1) + sum(case_control!=2))>1){
  frq$AF_UKN=apply(x[,case_control!=1 & case_control!=2],1,function(x){x=na.omit(x); return(round(sum(x)/(2*length(x)),6))})
  no$SAMPLES_UKN = apply(x[,case_control!=1 & case_control!=2],1,function(x){x=na.omit(x); return(length(x))})
 }

  
 if(hwe.bol & any(case_control==1) & sum(case_control==1)>1){
 hwe = data.frame( P_HWE_CONTROL = apply(x[,case_control==1],1,hwe.pval))
print(head(hwe))
}else{hwe=NA}

return(cbind(frq,no, P_HWE_CONTROL=hwe))
}


translate= function(ALT, REF, dos){
  dos[dos==2] = paste0(REF, REF)
  dos[dos==1] = paste0(ALT, REF)
  dos[dos==0] = paste0(ALT, ALT)
  dos[is.na(dos)]=paste0("00")
  return(dos)
}


