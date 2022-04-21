###################################################
# MAKE IMGT REFERENCE
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
###################################################

###################################################
# SETTINGS
####################################################
args = commandArgs(T)
util.dir=args[1]
alignments.dir=args[2]
script.dir=args[3]

options(stringsAsFactors=F)
library(readxl)

####################################################
# FUNCTIONS
####################################################


prepare = function(file, reference, locus){
# file: alignment file
# reference: reference string
# locus: HLA locus 

## READ IN ALLELE STATUS (uncomment next line if you only want to include
## confirmed alleles

 conf = read.table(file.path(dir, "Allele_status.txt"), sep=",", h=T)
# conf = conf[conf$Confirmed=="Confirmed",]

  cmd=paste("paste <(awk '{print $1}'", file,  
  " ) <(awk '{for (i=2; i<NF; i++) printf $i; print $NF}'", file, 
  ") >", paste0(file,".prep"))
  write.table(cmd, "bash.sh", col.names=F,row.names=F, quote=F)
  system("bash bash.sh")
  

# READ IN ALIGNMENT FILE
  tmp = read.table(paste0(file, ".prep"), h=F, sep="\t")
  tmp = tmp[tmp[,1]%in%conf[,1],]
  ind = grep("[0-9]", tmp[,2])
  if(any(ind)){
    tmp = tmp[-ind,]
  }
  # GET ORDER OF SEQUENCES
  order = tmp[grep(gsub("\\*","\\\\*", locus), tmp[,1]),1]

  # MAKE MATRIX OF NUC/PROT PER POSITION
  tmp = tapply(tmp[,2], tmp[,1], function(x){do.call(c,strsplit(paste(x, collapse=""),""))})
  names = names(tmp)
  l = unlist(lapply(tmp, length))
  n = max(l)

  # FILL UP PROT AFTER X
  tmp = lapply(tmp, function(x){if("X"%in%x){x = c(x, rep("*",n-length(x)))};return(x)}) # 

  # FILL UP IF ANY LENGTH != maximal length
  l = unlist(lapply(tmp, length))
  if(length(table(l))!=1){
    tmp = lapply(tmp, function(x){if(length(x)!=n){x = c(x, rep(".",n-length(x)))};return(x)}) 
  }
 # DELETE ROWS THAT DO NOT CONTAIN ALLELES 
  names = names[grep(gsub("\\*","\\\\*", locus), names)]
  tmp = tmp[names]
  tmp = do.call(rbind, tmp)
  print(order[!order%in%names])
  # ORDER SEQUENCES
  tmp = tmp[order,]
  
  # FILL UP MATRIX WITH THE FIRST ENTRY
  tmp = apply(tmp,2, function(x){x = gsub("-",x[1],x); x = gsub("\\*","0",x); return(x)})
  names(tmp) = names
  
  print(paste("Read file", file,"...."))
  print(paste("Reference is", reference))
  
  # REDUCE MATRIX TO POSITIONS OF REFERENCE STRING
  ref_string = tmp[reference, ]
  ind = ref_string=="." |ref_string=="|"
  tmp = tmp[,!ind]

  return(tmp)
}

flip=function(x){
  x=unlist(sapply(x,function(x){x=switch(x,A="T",T="A",C="G",G="C","0"="0","."="."); return(as.character(x))}))
  return(x)
}
## ASSIGN POSITIONS TO NUCLEOTIDES
nuc = function(data,locus, strand="+"){
  print(locus)
  start_end=read.table(file.path(util.dir, "PGF_start_end.txt"), h=T)
  if(locus=="DRB"){locus="DRB1*"}
  start_end = start_end[start_end$allele==locus,]
  bp= c()
  for(i in 1:nrow(start_end)){
    bp = c(bp, start_end[i,"start"]:start_end[i,"end"])
  }
  colnames(data)=bp
  # FLIP ALLELES THAT ARE ON THE - STRAND
  if(strand=="-"){
    data = apply(data,2,flip)
  }
  data[data=="."]="I/D"
  print(dim(data))
  return(data)
}

## ASSIGN POSITIONS TO AMINO ACIDS
prot= function(data,  reference, string){
  string = unlist(strsplit(string,""))
  print(string)
  i=1
  pos=1
  ref_string = data[reference, ]
  while(sum(ref_string[i:(i+length(string)-1)]==string)!=length(string)){
    i=i+1;
  }
  colnames(data)=c(-(i-1):-1, 1:(length(ref_string)-i+1))
  print(dim(data))
return(data)
}
   

####################################################
# MAIN
####################################################

# PGF REFERENCE ASSIGNMENTS, start (i.e. where is start of sequence);STRINGS and strands 
PGF=read_excel(file.path(util.dir, "TE_PGF.xls"), 1)[1:8,]
strings =c("GSHSMRYFFT","GSHSMRYFYT","CSHSMRYFDT","IKADHVSTYA",
           "RATPENYLFQ","EDIVADHVAS","RDSPEDFVFQ","GDTRPRFLWQ")

strand = c("+","-","-","-","+","+","-","-")

list_prot = list()
list_nuc = list()

for(i in (1:nrow(PGF))){
  locus = as.character(PGF[i,1])
  if(i==8){locus="DRB"}
  reference = paste0(PGF[i,1], PGF[i,2])
  for(suffix in c("nuc","prot")){
    # READ IN ALIGNMENT
    file = file.path(alignments.dir, "alignments", paste(gsub("\\*", "", locus),"_",suffix,".txt",sep=""))
    print(file)
    out = prepare(file, reference, locus)
    names = rownames(out)
    
    if(suffix=="nuc"){
      b = sapply(strsplit(names,":"), 
                 function(x){
                   if(length(x)==2){x=paste0(x[1],":",x[2])}; 
                   if(length(x)==3){x=paste0(x[1],":",x[2],":",x[3])}; 
                   if(length(x)==4){x=paste0(x[1],":",x[2],":",x[3])}; 
                   return(x)})
      ind = !duplicated(b)
      out = out[ind,]
      rownames(out) = b[ind]
      data = nuc(out,locus, strand = strand[i])
      data = data[,as.numeric(colnames(data))>29*10^6 & as.numeric(colnames(data))<35*10^6]
      print(dim(data))
      list_nuc[[paste(sub("\\*", "", locus),"_",suffix,sep="")]]= data
    } 
    if(suffix=="prot"){
      b = sapply(strsplit(names,":"), 
                 function(x){
                   if(length(x)==2){x=paste0(x[1],":",x[2])}; 
                   if(length(x)==3){x=paste0(x[1],":",x[2])}; 
                   if(length(x)==4){x=paste0(x[1],":",x[2])}; 
                   return(x)})
      ind = !duplicated(b)
      out = out[ind,]
      rownames(out) = b[ind]
      reference = sapply(strsplit(reference,":"), function(x){paste0(x[1],":",x[2])})
      data = prot(out, reference, strings[i])
      
      list_prot[[paste(sub("\\*", "", locus),"_",suffix,sep="")]]= data
    }
  }
}


prot = list_prot
nuc3dig = list_nuc

pdf("analysis_missingness_prot_full.pdf") 
for(i in names(prot)){
  ind = apply(prot[[i]],2,function(x){sum=length(x); x=length(which(x=="0"))/sum;} )
  plot(names(ind), ind*100, xlab="position", ylab="% missing", main = paste0(i) )
}
dev.off()

pdf("analysis_missingness_nuc3dig_full.pdf") 
for(i in names(nuc3dig)){
ind = apply(nuc3dig[[i]],2,function(x){sum=length(x); x=length(which(x=="0"))/sum;} )
plot(as.numeric(names(ind)), ind*100, xlab="position", ylab="% missing", main = paste0(i) )
}
dev.off()


# !!NOTE THAT ONLY PARTIAL SEQUENCES EXIST FOR SOME LOCI
source(file.path(script.dir, "make_groups.R"))
# FOR FULL: Always choose the most complete sequence (usually first sequence)
translate_to_lower = function(data,group,from,to){
  group = group[match(rownames(data), paste0(group$locus,"*",group[,from])),c("locus",from,to)]
  group[,from]=paste0(group[,"locus"],"*", group[,from])
  group[,to]=paste0(group[,"locus"],"*",group[,to])
  pos = colnames(data)
  data = split(data.frame(data), group[,to])
 a = lapply(data, function(x){
   if(nrow(x)>1){
     ind = apply(x,2,function(x){
       x=x[x!="0"]
       return(length(table(x))>1)}); x[,ind]="0"; x= x[1,]}
   return(x)
   })
 data = do.call(rbind, a)
 colnames(data) = pos
 
 ind = apply(data,2,function(x){sum=length(x); x=length(which(x=="0"))/sum;} )
 #PLOT

 plot(names(ind), ind*100, xlab="position", ylab="% missing" )

 return(data)
}

# FOR G: Always choose the most complete sequence (usually first sequence)
translate_to_g = function(data,group,from,to){
  group = group[match(rownames(data), paste0(group$locus,"*",group[,from])),c("locus",from,to)]
  group[,from]=paste0(group[,"locus"],"*", group[,from])
  group[,to]=paste0(group[,"locus"],"*",group[,to])
  pos = colnames(data)
  data = split(data.frame(data), group[,to])
  a = lapply(data, function(x){
    if(nrow(x)>1){
      ind = apply(x,2,function(x){
        return(length(table(x))>1)}); x[,ind]="0"; x= unique(x)}
    return(x)
  })
  
  data = do.call(rbind, a)
  colnames(data) = pos
  #PLOT
  ind = apply(data,2,function(x){sum=length(x); x=length(which(x=="0"))/sum;} )
  plot(names(ind), ind*100, xlab="position", ylab="% missing" )
  data = data[,ind< 0.05]
  return(data)
}

load("g_groups.RData")

pdf("analysis_missingness_nuc_full_2dig.pdf") 
nuc2dig = lapply(list_nuc, function(data){print(head(data[,1:10])); translate_to_lower(data,out,"X3field","X2field")})
dev.off()

save(nuc3dig,nuc2dig,prot,file="impute_SNPs_AA_full.RData")

####### TRANSLATE TO G GROUPS
pdf("analysis_missingness_nuc_3digG.pdf") 
nuc3dig =  lapply(list_nuc, function(data){translate_to_g(data,out,"X3field","X3Gfield")})
dev.off()
pdf("analysis_missingness_nuc_2digG.pdf") 
nuc2dig = lapply(list_nuc, function(data){translate_to_g(data,out,"X3field","X2Gfield")})
dev.off()
pdf("analysis_missingness_prot_protG.pdf") 
prot = lapply(prot, function(data){translate_to_g(data,out,"X2field","X2Gfield")})
dev.off()

save(nuc3dig,nuc2dig,prot,file="impute_SNPs_AA_G.RData")

