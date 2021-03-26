#!/usr/bin/env Rscript

##############################################
# HLA phasing
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# April 2020
##############################################

################################
# Settings
################################
args = commandArgs(trailingOnly=T)

library("reshape2")
library("parallel")

print(args)
################################
# FUNCTIONS
################################

to.dosage = function(x, maj){
  miss = x=="00";  
  x=gsub(maj,"",x); 
  x=sapply(x,nchar);
  x[miss]=NA
  
  return(x)
}

convert = function(data){
  data$id = gsub(":.*","", data$id)
  data[data$ALT=="A",-(1:8)] = 2-   data[data$ALT=="A",-(1:8)]
  data[data$ALT=="A",c("ALT","REF")] =  data[data$ALT=="A",c("REF","ALT")]
  info = unique(data[,1:8])
  no = table(data$id)
  #no = no[no > 1]
  
  data = split(data,data$id)
  data = lapply(data, function(x){colSums(x[,-(1:8)])})
  data = do.call(rbind,data)
  data = cbind(info[match(rownames(data),info$id),], data)
  data = data[data$id%in%names(no),]
  return(data)
}



assign.HLA=function(data, fam, position){
  
  
  data[data$postprob < 0,c("A.name","B.name")] =NA
  data = melt(data[,c("id","locus", "A.name","B.name")],id=c("id","locus"))
  
  dos=c()
  
  for(i in unique(data$locus)){
    tmp = data[data$locus==i,]
    
    tmp$value=factor(tmp$value)
    tmp= as.data.frame.matrix(table(tmp$value, tmp$id))
    tmp[,colSums(tmp)==0]=NA
    dos=rbind(dos,tmp)
  }
  dos = dos[,match(fam$IID,colnames(dos))]

  position$mid = position$start+ceiling((position$end -position$start+1)/2)
  
  out = cbind(chr=6, id=paste0("imputed_HLA_", rownames(dos)),
                   pos=position[match(gsub("\\*.*","",rownames(dos)),position$locus),"mid"],
                   REF="A",ALT="P",INFO=NA, TYPE=0, source="HIBAG_HLA",dos)
  
  
  return(out)
}
# GET SNP/AA from HLA data
assign.SNP.PROT_prepare=function(data, info, fam, type){   
  
  data[data$postprob < 0,c("A.name","B.name")] =NA
  id = unique(data$id)
  
  out = mclapply(as.list(id),function(x){
    all = data[data$id==x,c("locus","A","B")]
     out= c()
    for(i in 1:nrow(all)){
      locus=all[i,"locus"]
      if(locus%in%c("DRB1","DRB3","DRB4","DRB5")){
        lib_loc = "DRB"
      }else{
        lib_loc=locus
      }
      gen = paste(locus,"*",as.character(all[i,c("A","B")]),sep="")
      gen = gsub("G", "", gen)
      lib = info[[paste(lib_loc,"_",type,sep="")]]
      colnames(lib)=paste(locus,"_",colnames(lib),sep="")
      rownames(lib)=gsub("[NLQEG]$","", rownames(lib)) # Only because gsubed in source
      if(any(!gen%in%rownames(lib))){
        tmp = lib
        if(!gen[1]%in%rownames(tmp)){
          tmp= rbind(tmp, rep(0, ncol(tmp)))
          rownames(tmp)[nrow(tmp)] = gen[1]
          print(paste(x, "Missing", gen[1]))
        }
        if(!gen[2]%in%rownames(tmp)){
          tmp= rbind(tmp, rep(0, ncol(tmp)))
          rownames(tmp)[nrow(tmp)] = gen[2]
           print(paste(x, "Missing", gen[2]))
        }
       # print(paste(x, "Missing", gen))

        out = cbind(out, as.matrix(tmp[gen,]))
      }else{
        out = cbind(out, as.matrix(lib[gen,]))
      }  }
    
    out = apply(out,2,function(x){x=toupper(x);
    if(type=="nuc"){
      x=sort(x);
    }
    x=paste(x,collapse="")
    x[grep("0",x)]="00"
    
    return(x)})
    return(out)}, mc.cores=detectCores(), mc.preschedule=F, mc.silent=F)
  
  out = do.call(rbind,out)

  if(type=="nuc"){
    out = out[,order(as.numeric(gsub(".*_|I","",colnames(out))))]
    out=  out[,as.numeric(gsub(".*_|I","",colnames(out)))>=25*10^6]  
  }
  
  out  = t(out)
  colnames(out) = id
  
  out = out[,match(fam$IID, colnames(out))]  
  out = cbind(rownames(out),rownames(out),out)
  out[,1] = paste0("imputed_",type,"_", out[,1])

  return(out)
 
}

assign.SNP.PROT=function(tmp){   
  # MAKE TO DOSAGE
  if(any(-grep("[I]",out[,1]))){
    tmp = out[-grep("[I]",out[,1]),] 
  }
  #
  row.names=tmp[,1]
  
  tmp = tmp[,-(1:2)]
#  print(head(tmp)[,1:10])
  ###################################
  # Get the major and minor allele
  ##################################
  maj=apply(tmp,1,function(x){x=(unlist(strsplit(x,""))); 
  
  x[x==0]=NA; 
  x=(sort(table(x), decreasing=T))
  x=names(x)
  if(length(x)==0){ #If monomorphic set 0
            return(c("0","0"))
  }
  if(length(x)==1){ #If monomorphic set 0
    return(c(x,"0"))
  }
  if(length(x)==2){ #If biallelic return both alleles sorted
    return(x)
  }
  if(length(x)>2){ #Set InDels to missing
    return(c(NA,NA))
  }})                                    
  
  
  maj=t(maj)
  
  rownames(maj) =row.names
  
  ##################################
  # Calculate dosage on minor allele
  #################################
  
  data = apply(cbind(row.names,tmp),1,function(x){
    al = maj[x[1], ]
    res = c()
    if(!is.na(al[1]) & al[1]!="0"){
      res = c(6,x[1],gsub(".*_", "",x[1]),al,to.dosage(x[-1],al[1]))
      return(res)
    }else{
      return(NULL)
    }
  }    
  )
  data = do.call(rbind,data)
  
  
  data = cbind(data[,1:5],NA,0,"HIBAG_SNPs/AA", data[,-(1:5)])
  colnames(data)[1:8] = c("chr","id","pos","REF","ALT","INFO","TYPE","source")
  
  return(data)
} 

translate= function(ALT, REF, dos){
  dos[dos==2] = paste0(REF, REF)
  dos[dos==1] = paste0(ALT, REF)
  dos[dos==0] = paste0(ALT, ALT)
  dos[is.na(dos)]=paste0("00")
  return(dos)
}

format_to_plink = function(fam, data){
  tr = apply(data, 1, function(x){
    translate(x[4], x[5], x[-(1:8)])
  })
  ped = cbind(fam, tr[match(rownames(tr),fam$IID),])
  map = cbind(data[,c(1,2)], 0, data[,3])
  map[,2] = gsub("\\*|:","_", map[,2])
  return(list(ped=ped, map=map))
}
save(args, file="args.RData")
################################
# MAIN
################################

### A: SINGLE ALLELE LOCI

## Concatanate Gene Loci 
load(args[1])
data = lapply(pred, function(x){x$value})
data = cbind(rep(names(pred), times = unlist(lapply(data, nrow))), do.call(rbind, data))
data = data.frame(data)
colnames(data)=c("locus","id","A","B","postprob")

## Match *RData to fam file
fam= read.table(paste0(args[2],".fam"), h=F, col.names = c("FID","IID","","","","PHENO"))

data$A.name= paste(data$locus,data$A,sep="*")
data$B.name= paste(data$locus,data$B,sep="*")



position=read.table(args[3],h=T)


head(data)
out = assign.HLA(data, fam, position)
HLA = rbind(out, convert(out))


save = out
write.table(HLA, paste("imputation_",args[2],".data", sep=""), quote=F, sep="\t", row.names=F)

data$A=unlist(sapply(strsplit(data$A, ":"), function(x){x=paste(x[1], x[2],sep=":")})) 
data$B=unlist(sapply(strsplit(data$B, ":"), function(x){x=paste(x[1], x[2],sep=":")}))
data$A.name = unlist(sapply(strsplit(data$A.name, ":"), function(x){x=paste(x[1], x[2],sep=":")})) #remove 3(4) -field for following analysis
data$B.name = unlist(sapply(strsplit(data$B.name, ":"), function(x){x=paste(x[1], x[2],sep=":")}))

print(head(data))

load(args[4])
if(any(-grep("E", data$A.name))){
  data = data[-grep("E", data$A.name),]
}




out = assign.SNP.PROT_prepare(data, nuc2dig, fam, "nuc")

out = assign.SNP.PROT(out)
NUC=out
write.table(out, paste("imputation_",args[2],".nuc.data", sep=""), quote=F, sep="\t", row.names=F)
  
out = assign.SNP.PROT_prepare(data, prot, fam, "prot")
out = assign.SNP.PROT(out)
PROT=out
  
write.table(out, paste("imputation_",args[2],".prot.data", sep=""), quote=F, sep="\t", row.names=F)
  
plink =  rbind(HLA, NUC, PROT)

plink = format_to_plink(fam, plink)

write.table(plink$map, paste("imputation_",args[2],".map", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
write.table(plink$ped, paste("imputation_",args[2],".ped", sep=""), quote=F, sep="\t", row.names=F, col.names=F)

cmd = paste("plink --file", paste0("imputation_",args[2]),
            "--make-bed --out", paste0("imputation_",args[2]))
system(cmd)
                                    
