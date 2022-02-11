#!/usr/bin/env Rscript

##############################################
# # Retreive additional information on the data 
# i.e. marginal posterior probablity
# i.e. allele frequencies
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################


##############################################
# Settings
##############################################
options(stringsAsFactors=F)
library(ggplot2)
library(data.table)
library(reshape)
args= commandArgs(T)

##############################################
# Functions
##############################################

########################################
# CONCATANATE DATA INTO ONE MATRIX WITH
# MARIGNAL PROBABILITY FOR ALL ALLELES
# AND SAMPLES
########################################
concat = function(x, digits=4){
  tmp = x$postprob
  # GET marginal probabilitites
  id = do.call(rbind, strsplit(rownames(tmp),"/"))
  di = id[,c(2,1)]
  ind = di[,1] != di[,2]
  # Double rows, so that each allele occurs the same number of times
  info = rbind(id, di[ind,])
  out = data.frame(rbind(tmp,tmp[ind,]))
  
  if(digits==2){ # DELETE DOUBLES & PREPARE FOR CALCULATION OF 2 digits
    a= gsub(":.*","$", info[,1])
    id= split(data.frame(info[,1:2]),a)
    id = lapply(id, function(x){x = apply(x, 1,sort);return(unique(t(x)))})
    for(i in names(id)){
      ind = grep(gsub("\\$", ":", i), id[[i]][,2])
      id[[i]][ind,1:2] =id[[i]][ind, 2:1] 
    }
    id = do.call(rbind, id)
    ind =paste(info[,1],info[,2])%in%paste(id[,1], id[,2])
    info= info[ind,]
    out = out[ind,]
    info[,1]= gsub(":.*","$", info[,1])
    info[,2]= gsub(":.*","$", info[,2])
  }
  out = lapply(split(out, info[,1]), colSums)

  # MARGINAL probability matrix
  info = do.call(rbind, out)
  colnames(info) = colnames(tmp)
  return(info)
}

#############################################
# CALCULATE MEAN PROB ACROSS ALL SAMPLES
#############################################
get_prop = function(val, info, digits=4){

  val = melt(val[,1:3], id=c("sample.id"))
  if(digits==2){
    val$value = gsub(":.*","$",val$value)
  }
  info = melt(info, id=row.names)
  info = info[match(paste(val$sample.id, val$value), paste(info$X2,info$X1)),]
  info = tapply(info[,3], as.matrix(info[,1]), mean)
  
  return(info)
  
}

get_freq =  function(val, digits=4){

  val = melt(val[,1:3], id=c("sample.id"))
  if(digits==2){
    val$value = gsub(":.*","$",val$value)
  }
  freq = table(val$value)/nrow(val)
  return(freq)
}

prep_prob = function(files){
  val = list()
  info = list()
  info.small = list()
  res = list()
  freq = list()

  ## PREPARE FOR ANALYSIS
  for(i in files){
    load(i)
    for(loc in names(pred)){
      
      val[[loc]] = rbind(val[[loc]],pred[[loc]]$value)
      info[[loc]] = cbind(info[[loc]], concat(pred[[loc]]))
      info.small[[loc]] = cbind(info.small[[loc]], concat(pred[[loc]], digits=2))
    }
  }
  ## GET MARGINAL PROPABILITIES
  for(loc in names(pred)){
    
   freq[[loc]] = c(get_freq(val[[loc]]), get_freq(val[[loc]], digits=2))                      
    
  }
  
  for(loc in names(pred)){
    res[[loc]] = c(get_prop(val[[loc]],info[[loc]]),
                    get_prop(val[[loc]], info.small[[loc]], digits=2))
  }                
                   
  res = cbind(locus=rep(names(res), times = unlist(lapply(res,length))), 
               id=  unlist(lapply(res, names)),
               prob = unlist(res),
               freq = unlist(freq))
  res = data.frame(res)
  
return(res)
}

##############################################
# Main
##############################################
#files = list.files(imputed.dir, pattern=pattern, full.names=T)
res = prep_prob(args[1])

res$prob = as.numeric(res$prob)

png("unsure.png")
res$digits=NA
res$digits[-grep(":", res$id)] ="2-digit"
res$digits[grep(":", res$id)] ="4-digit"

p = ggplot(res, aes(x=locus, y=prob)) + geom_boxplot() + facet_wrap(.~digits) + theme_bw() +  ggtitle("Marginal posterior probabilities") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0,1))
print(p)
dev.off()


res$id= gsub("\\$","", res$id)
write.table(res,paste("marginal_prob_", args[2],  ".txt",sep=""), sep="\t", row.names=F, quote=F, col.names=T)

# GET INFO ON ALLELE FREQUENCIES
source(args[3])
data = read.table(paste0("imputation_",args[2],".data"), h=T)
data = data.frame(data)
fam= read.table(paste0(args[2],".fam"), h=F, col.names = c("FID","IID","","","","PHENO"))

frq=make.frq.hwe(data[,-(1:8)], as.numeric(as.matrix(fam$PHENO)))
print(warnings())
write.table(cbind(data[,1:8],frq), paste0("imputation_", args[2],".info"), quote=F, sep="\t", row.names=F)
