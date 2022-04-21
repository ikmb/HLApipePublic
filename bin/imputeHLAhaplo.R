#!/usr/bin/env Rscript
############################################
# HAPLOTYPE OVERLAP HLA DATA
# AUTHOR: Frauke Degenhardt
# Mai 2020
# email: f.degenhardt@ikmb.uni-kiel.de
############################################
############################################
# SETTINGS
############################################
args = commandArgs(T)
gene =args[1]
model = args[2]
name=args[3]
print(args)
############################################
# FUNCTIONS
############################################
loadObject = function(x){
  model = local({
    load(x)
    stopifnot(length(ls())==1)
    environment()[[ls()]]
  } 
  )
  return(model)}


############################################
#  MAIN
############################################

HLAModelList = loadObject(model)

tmp = HLAModelList[[gene]]$classifiers

# GET ALL POSSIBLE HLA ALLELES
hla = c()
for(i in 1:100){
  out=tmp[[i]]$haplo
  hla = c(hla, out$hla)
}
hla = unique(hla)
library(parallel)



# GET POSSIBLE COMBINATIONS OF ALLELES
res = list()
counter=0
used = c()
ch = choose(length(hla), 2)
print(paste("There are",ch, "possible combinations for", gene ))
pos.model = HLAModelList[[gene]]$snp.position
if(length(args)==4){
bim = read.table(paste0(args[4],".bim"),h=F, col.names=c("chr","","","pos","A","B"))
pos = bim$pos
print(table(unique(pos.model)%in%pos))
write.table(pos.model[pos.model%in%pos],paste(gene,"_", name,"_pos_overlap.txt", sep=""),sep="\t", quote=F, row.names=F, col.names=F)
}
# GET OVERLAP BETWEEN HAPLOTYPES FOR THESE ALLELES OVER 100 CLASSIFIERS
# IF POS is defined, exclude POS not present in POS
for(one in hla){
  used = c(one, used)
  for(two in hla){
    if(!two%in%used){
      counter = counter+1
      if(counter%%100==0){
        print(counter)
      }
      list=mclapply(as.list(1:100), function(x){
        out=tmp[[x]]$haplos
        
        # GET OVERLAPPING POSITIONS IF POS IS DEFINED
        if(length(args)==4){
          haplo = do.call(rbind, strsplit(out$haplo,""))
          colnames(haplo) = pos.model[tmp[[x]]$snpidx]
          haplo = apply(haplo[,colnames(haplo)%in%pos],1,function(x){paste(x, collapse="")})
         out$haplo=haplo 
        }
        
        out = out[out$hla%in%c(one,two),]
         
        # COMBINE ALLELS WITH SAME HAPLOTYPE
        out= tapply(out$hla,out$haplo,function(x){paste(x,collapse=",")})
        return(out)}, 
                    mc.cores = detectCores(), mc.preschedule=F)
      res[[paste(one,two)]] =  lapply(list, function(x){length(grep(",",x))/length(x)})  
    }
  }
}

# INTERMEDIATE SAVE
save(res, file= paste(gene,"_",name,".RData",sep=""))
 

# CALCULATE STATISTICS OF OVERLAP
out = lapply(res, function(x){
     x=unlist(x)
     x= c(quantile(x, na.rm=T)*100, length(x[x>0]))
     x=round(x,3)
     if(x[6]>0){ return(x)}
   
    })
    
out=do.call(rbind,out)
names = colnames(out)
out = data.frame(gene=gene,rownames(out), out)
colnames(out) = c("gene", "alleles", c("0%","25%","50%","75%","100%", "NClassifiers"))
# OUTPUT

write.table(out, paste(gene,"_", name,".txt", sep=""),sep="\t",row.names=F, col.names=T,quote=F)

