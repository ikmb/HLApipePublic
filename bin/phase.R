#!/usr/bin/env Rscript

##############################################
# HLA phasing
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################


library(reshape2)
library(parallel)
library(data.table)

################################
# Settings
################################
args = commandArgs(trailingOnly=T)

options(stringsAsFactors=F)
name   = args[1]
pred   = args[2]
haps   = args[3]
sample = args[4]
certainty = args[5]
model = args[6]
gen = args[7]
out.dir="."
# 
save(args, file="args.RData")

################################
#  Setttings
################################
print(args)

print(paste("Reading phased data",haps))
phased = read.table(haps, h=F)
samples = read.table(sample,h=T)

print(paste("Reading certainty data",certainty))
certainty = read.table(certainty, h=T, sep="\t")


rownames(certainty) = gsub("chr6:", "",rownames(certainty))
colnames(certainty) = as.matrix(rep(samples$ID_2[-1],each=2))
colnames(phased) = c("chr","id","pos","al1","al2", as.matrix(rep(samples$ID_2[-1],each=2)))
rownames(phased) =phased$pos

# phased=phased.tmp
# haplotypes=snp.haplotypes.model
# predicted=pred.value
# certainty=certainty.tmp

################################
# Function
################################
assign = function(phased,haplotypes,  predicted, certainty){
  id = as.list(unique( predicted[,1]))
  chromosome=do.call(rbind, mclapply(id, function(x){
    sample.id=x
    # Get sample in sample.id
    
    
    phased = phased[,colnames(phased)==sample.id]
    certainty = certainty[,gsub("\\.", "-", gsub("X0:", "", gsub("\\.1","",colnames(certainty))))==sample.id] #TODO
    certainty = round(apply(certainty,2,mean),2)
    # Get alleles of this sample
    predicted =   predicted[predicted$sample.id==sample.id,"value"]
    
    if(predicted[1]!=predicted[2]){
      if(predicted[1]%in%haplotypes$allele & predicted[2]%in%haplotypes$allele){
        haplotypes= haplotypes[haplotypes$allele%in%predicted,]
        sum = c()
        # Compare to genotype of this sample
        for(i in 1:nrow(haplotypes)){
          sum = rbind(sum,colSums(abs( phased -as.numeric(as.matrix(t(haplotypes[i,-(1:2)]))))))
        }
        rownames(sum) = haplotypes$allele   
        
        sum.names = rownames(sum)
        out = apply(sum,2,function(x){min = min(x); x =paste(unique(sum.names[which(x==min)]), collapse="/"); return(c(x, min))})
        out = rbind(out, certainty)
        chromosome = c(as.vector(out),  nrow(phased), sample.id)
      }else{
        chromosome = c(rep(NA,7), sample.id)
      }
     
      

    }else{
      chromosome = c(predicted[1],-1,certainty[1],predicted[2],-1,certainty[2], nrow(phased), sample.id)
    }},
                                     mc.cores=8,mc.preschedule = F ))
  
  return(chromosome)
}

#########################################
#  Get called HLA types
#########################################

print(paste("Loading predicted data",pred))
load(pred)
pred= pred$value 
pred.value = melt(pred[,1:3], id="sample.id")

##########################################
# Get allele list and name with positions
##########################################
loadObject = function(x){
  model = local({
    load(x)
    stopifnot(length(ls())==1)
    environment()[[ls()]]
   
  } 
  )
  return(model)
  }


print(paste("Loading model",model))
model=loadObject(model)

#specify subset of classifier
classifier=sample(1:length(model[[gen]]$classifiers),10, replace=F)


# Get alleles and position of SNPs
model.snp.info= do.call(rbind,strsplit(model[[gen]]$snp.allele,"\\/"))
rownames(model.snp.info)=model[[gen]]$snp.position
vote = c()


print(paste("Phasing allele", gen, "..."))
###########################################
# Get Haplotypes and indeces of SNPs
###########################################
for(i in classifier){
  
  print(paste("Opening classifier",i))
  # EXTRACT SNP HAPLOTYPES
  snp.haplotypes.model = model[[gen]]$classifiers[[i]]$haplo
  snp.haplotypes.model = cbind(snp.haplotypes.model[,1:2], 
                               do.call(rbind,strsplit(snp.haplotypes.model[,3],"")))
  indx = model[[gen]]$classifiers[[i]]$snpidx
  colnames(snp.haplotypes.model) = c("freq","allele",rownames(model.snp.info)[indx])
  
  # EXTRACT SNPS ALSO PRESENT IN PHASED
  snp.haplotypes.model= cbind(snp.haplotypes.model[,1:2], 
                              snp.haplotypes.model[,colnames(snp.haplotypes.model)%in%phased$pos])
  
  phased.tmp = phased[match(colnames(snp.haplotypes.model[,-(1:2)]),phased$pos),]
  certainty.tmp= certainty[match(phased.tmp$pos, rownames(certainty)),]
  
  # GET Alleles SNP HLA
  model.snp.info.tmp = model.snp.info[match(phased.tmp$pos,rownames(model.snp.info)),]
  
  
  #############################################
  # Get SNPs to flip (not real flip, just means 0/1 and 1/0 coding of HIBAG and SNP2HLA)
  #############################################
  process.snps = cbind(same = model.snp.info.tmp[,1]==phased.tmp[,5] & model.snp.info.tmp[,2]==phased.tmp[,4], 
                      flip  = model.snp.info.tmp[,1]==phased.tmp[,4] & model.snp.info.tmp[,2]==phased.tmp[,5] )

  
  process.snps = data.frame(process.snps)
  process.snps$use =process.snps$same |process.snps $flip
 
  ########################
  # Print difference
  ########################
  
  phased.tmp[process.snps$flip,-c(1:5)] = abs(phased.tmp[process.snps$flip,-c(1:5)]-1)
  
  ########################
  # Print differences
  ########################
  print("Using SNPs")
  print(process.snps)
  phased.tmp = phased.tmp[process.snps$use ,]
  snp.haplotypes.model = snp.haplotypes.model[,c(T,T,process.snps$use)]
  certainty.tmp =  certainty.tmp[process.snps$use,]
  
  print(paste("Getting votes for classifier", i))
  out = cbind(assign(phased.tmp, snp.haplotypes.model, pred.value,  certainty.tmp),i)
  out[grep("Error",out)]=NA
  print(out)
  out[is.na(out[,8]),8] = unique(pred$sample.id[!pred$sample.id%in%out[,8]])
  vote = rbind(vote,out)

}



colnames(vote)=c("1","min_1", "shapeit_1", "2", "min_2","shapeit_1", "number.pos","id","classifier")
vote = data.frame(vote)
save(vote,file="vote.RData")

final_vote = tapply(paste(vote[,1], vote[,4]), vote$id, function(x){
 # Get most common combination 
 x=sort(table(x), decreasing=T); 
 al = unlist(strsplit(names(x[1]), " "))
 if(al[1]==al[2]  | length(grep("/", al[1])) + length(grep("/", al[2]))>0){ #If chromsomes equal not good
   return(c(al[1],al[2],0))
   }else{
     return(c(al[1],al[2], x[1]/10))}
})
final_vote = do.call(rbind,final_vote)
ind = paste(vote$id,vote$X1, vote$X2)%in%paste(rownames(final_vote),final_vote[,1], final_vote[,2])
range = apply(vote[,c(2,3,5:7)][ind,],
             2, function(x){
 out =tapply(x, vote$id[ind], function(x){x=as.numeric(as.matrix(x)); x=quantile(x, c(0,0.5,1),type=6, na.rm=T);
                                return(paste(round(x[2],2),";[",round(x[1],2),",", round(x[2],2),"]",sep=""))}); 
 return(out)})

# Set phasing certainty for homozygous to 1
homozygous = unique(vote[vote$min_1=="-1" &!is.na(vote$min_1),"id"])
final_vote[homozygous,3]=1

output = cbind(final_vote, range)
pred.prob=pred= pred[match(samples$ID_2[-1],pred$sample.id),"prob"]

output = data.frame(output[match(samples$ID_2[-1], rownames(output)),])

out=cbind(rownames(output),gen,output, pred.prob)
colnames(out)=c("id","locus","1","2","phase_prob", "min_diff_1", "shapeit_1","mind_diff_2","shapeit_2","pos_used","imp.prob")

write.table(out,file.path(out.dir,paste(name,".",gen,".HLA.phased.txt",sep="")), 
           col.names=T, row.names=F, quote=F,sep="\t")

write.table(vote,file.path(out.dir,paste(name,".",gen,".HLA.info.txt",sep="")), 
           col.names=T, row.names=F, quote=F,sep="\t")

