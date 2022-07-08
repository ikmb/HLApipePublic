###################################################
# MAKE IMGT REFERENCE
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
###################################################

###################################################
# SETTINGS
####################################################
options(stringsAsFactors=F)

library(readxl)

args = command.args(T)

aa_snp.reference.dir = args[1]

alignment.dir = args[2] ## Download latest alignment form the IDP-IMGT/HLA database

PGF = read_excel(file.path(aa_snp.reference.dir, "PGF.xls"), 1)[1:8,]
start_end = read.table(file.path(aa_snp.reference.dir,"PGF_start_end.txt"), h=T)

source(file.path(aa_snp.reference.dir,"make_groups.R"))
groups=out
####################################################
# FUNCTIONS
####################################################

prepare = function(file, reference, locus, status.dir=NULL){
  ## Optional comment in: Take only alleles with confirmed HLA allele status
  if(!is.null(status.dir)){
    conf = read.table(file.path(status.dir,"Allele_status.txt", sep=",", h=T))
    conf = conf[conf$Confirmed=="Confirmed",]
  }

  ## Prepare for alignment
  cmd=paste("paste <(awk '{print $1}'", file,  
  " ) <(awk '{for (i=2; i<NF; i++) printf $i; print $NF}'", file, 
  ") >", paste0(file,".prep"))
  write.table(cmd, "bash.sh", col.names=F,row.names=F, quote=F)
  system("bash bash.sh")

  ## Read in allignment file
  tmp = read.table(paste0(file, ".prep"), h=F, sep="\t")
  if(!is.null(status.dir)){
  tmp = tmp[tmp[,1]%in%conf[,1],]}
    
  ## Delete empty rows
  ind = grep("[0-9]", tmp[,2])
  if(any(ind)){
    tmp = tmp[-ind,]
  }
  
  ## Get the order of sequences
  order = tmp[grep(gsub("\\*","\\\\*", locus), tmp[,1]),1]

  ##  MAKE MATRIX OF NUC/PROT PER POSITION
  tmp = tapply(tmp[,2], tmp[,1], function(x){do.call(c,strsplit(paste(x, collapse=""),""))})
  names = names(tmp)
  l = unlist(lapply(tmp, length))
  n = max(l)
  ## FILL UP PROT AFTER X
  tmp = lapply(tmp, function(x){if("X"%in%x){x = c(x, rep("*",n-length(x)))};return(x)}) # 

  # FILL UP IF ANY LENGTH != maximal length
  l = unlist(lapply(tmp, length))
  if(length(table(l))!=1){
    tmp = lapply(tmp, function(x){if(length(x)!=n){x = c(x, rep(".",n-length(x)))};return(x)}) 
  }
  ## DELETE ROWS THAT DO NOT CONTAIN ALLELES 
  names = names[grep(gsub("\\*","\\\\*", locus), names)]
  tmp = tmp[names]
  tmp = do.call(rbind, tmp)
  print(order[!order%in%names])
  ## ORDER SEQUENCES
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

nuc = function(data,locus, start_end,strand="+"){
  print(locus)

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
  data[data=="."]="I/D"
return(data)
}
   
## Translate to different resolution or group
# !!NOTE THAT ONLY PARTIAL SEQUENCES EXIST FOR SOME LOCI

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
  ind = apply(data,2,function(x){sum=length(x); x=length(which(x=="0"))/sum;} )
  data = data[,ind< 0.05]
  return(data)
}


plotfunction = function(data, title){
  pdf(paste0(title, ".pdf"))
  for(i in names(data)){
    ind = apply(data[[i]],2,function(x){sum=length(x); x=length(which(x=="0"))/sum;} )
    plot(names(ind), ind*100, xlab="position", ylab="% missing", main = paste0(i) )
  }
  dev.off()
}

####################################################
# MAIN: Create libraries
####################################################


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
    
    # read in alignment
    file = file.path(alignment.dir, paste(gsub("\\*", "", locus),"_",suffix,".txt",sep=""))
    out = prepare(file, reference, locus)
    names = rownames(out)
    
    ## prepare nucleotide library
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
      data = nuc(out,locus, start_end, strand = strand[i])
      data = data[,as.numeric(colnames(data))>29*10^6 & as.numeric(colnames(data))<35*10^6]
      print(dim(data))
      list_nuc[[paste(sub("\\*", "", locus),"_",suffix,sep="")]]= data
    } 
    
    ## prepare protein library
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

###############################################################################
# MAIN: Modify libraries
################################################################################

## full context
prot = list_prot
nuc3dig = list_nuc
nuc2dig = lapply(list_nuc, function(data){print(head(data[,1:10])); translate_to_lower(data,groups,"X3field","X2field")})
save(nuc3dig,nuc2dig,prot,file="impute_SNPs_AA_full.RData")

plotfunction(prot, "prot")
plotfunction(nuc3dig, "nuc3dig")
plotfunction(nuc2dig, "nuc2dig")

## g groups
load("g_groups.RData")
nuc3dig =  lapply(list_nuc, function(data){translate_to_g(data,groups,"X3field","X3Gfield")})
nuc2dig = lapply(list_nuc, function(data){translate_to_g(data,groups,"X3field","X2Gfield")})
prot = lapply(prot, function(data){translate_to_g(data,groups,"X2field","X2Gfield")})


plotfunction(prot, "protG")
plotfunction(nuc3dig, "nuc3digG")
plotfunction(nuc2dig, "nuc2digG")

save(nuc3dig,nuc2dig,prot,file="impute_SNPs_AA_G.RData")

# 
# ### Sanity check
# info = nuc2dig
# lengths = c()
# for(i in 1:8){
#   tmp = info[[i]]
#   tmp = tmp[1:ceiling(length(tmp)/2),]
#   check = apply(tmp,2,function(x){x[x=="0"]=NA; return(table(as.character(x)))})
#   lengths = rbind(lengths,
#                   cbind(names(info)[i],unlist(lapply(check,length)),
#                         unlist(lapply(check, function(x){paste(names(x),collapse=",")}))))
# }
# 
# lengths=data.frame(lengths)
# indel = lengths[grep("I", rownames(lengths)),]
# lengths = lengths[as.numeric(lengths[,2])>2,]
# colnames(lengths)=c("loc", "no.alleles","alleles")
# lengths = data.frame(cbind(pos=rownames(lengths), lengths))
# lengths$alleles=unlist(lapply(strsplit(as.character(lengths$alleles),","), function(x){x=gsub("I/D|D/I","-", x); x=paste(unique(sort(x)), collapse=",")}))
# library(data.table)
# mult=fread("hgTables_snpdb150.txt", h=T,sep="\t")
# mult[mult$strand=="-","observed"]=sapply(mult[mult$strand=="-","observed"], function(x){chartr("ATCG/","TAGC/", x)})
# mult= tapply(mult$observed,mult$chromEnd, function(x){paste(x,collapse="/")})
# 
# mult = mult[match(lengths$pos,paste0(names(mult)))]
# mult = unlist(lapply(strsplit(mult,"/"), function(x){sort(x); paste(unique(sort(x)), collapse=",")}))
# 
# 
# lengths$known_in_UCSC = mult
# lengths = lengths[!is.na(lengths$known_in_UCSC) &!lengths$known_in_UCSC=="",]
# table(lengths$alleles==lengths$known_in_UCSC)
# write.table(lengths, "multallelic_variants_HLA_reference.csv", row.names=T, col.names=T,quote=T,sep="\t")




