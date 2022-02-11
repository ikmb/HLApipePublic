#!/usr/bin/env Rscript

##############################################
# Compare study to reference
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################
################################
# Settings
################################


args = commandArgs(T)
model.name = args[1]
params.rootname = args[2]
print(args)
print(length(args))
library(ggplot2)
library(grid)
library(gridExtra)

options(stringsAsFactors=F )
loadObject = function(x){
  model = local({
    load(x)
    stopifnot(length(ls())==1)
    environment()[[ls()]]
    
  } 
  )
  return(model)}

##############################
# GET FREQUENCIES IN MODEL
##############################

model = loadObject(model.name) # Load Model
reference = do.call(rbind,
             lapply(model, 
                    function(x){return(cbind(x$snp.position, do.call(rbind, strsplit(x$snp.allele,"/")), 
                                                                      x$snp.allele.freq))}))
# FOR MULTIPLE POS (USED MORE T 1 ACROSS ALL LOCI) CALCULATE MEAN OF MAF
mean_MAF = tapply(reference[,4], paste(reference[,1], reference[,2], reference[,3],sep="_"), function(x){ mean(as.numeric(as.matrix(x)))})
reference = cbind(do.call(rbind, strsplit(names(mean_MAF), "_")), unlist(mean_MAF))

#####################################################
# IF SPECIFIC POPULATION: GET FREQUENCY IN POPULATION
#####################################################

if(length(args)>2){
  tmp = reference

  reference.file = args[3]
  reference = read.table(reference.file,h=T)
  reference = cbind(reference[,c(3,5,6)], reference[,args[4]])
  reference = reference[reference[,1]%in%tmp[,1],]
  
}

colnames(reference) = c("pos","A1","A2","MAF")

reference = data.frame(reference)
reference$MAF = as.numeric(as.matrix(reference$MAF))
ind = reference$MAF >0.5
reference[ind, c("A1","A2")] = reference[ind, c("A2","A1")]
reference[ind, c("MAF")] = 1-   reference[ind, c("MAF")] 
print(table(ind))
colnames(reference) = c("pos","A1","A2","MAF")


print(length(reference$MAF))
#####################################################
# GET FREQUENCY IN STUDY
#####################################################

cmd = paste("plink --bfile", params.rootname, "--freq --out tmp")
system(cmd)

tmp = cbind(read.table(paste(params.rootname,".bim", sep=""),                     
                        col.names=c("chr","id","","pos","A1","A2"), colClasses="character"),
             read.table("tmp.frq", h=T, colClasses="character"))


tmp = tmp[tmp[,1]==6,]
save(tmp, reference, file="before.RData")
############################
# MATCH
###########################
print(head(reference))

reference = reference[reference$pos%in%tmp$pos,]
print(dim(reference))
tmp       = tmp[match(reference$pos, tmp$pos),]
table(!is.na(tmp$pos))

save(tmp, reference, file="test.RData")
##############################
# FUNCTIONS
##############################


flip = function(x){
  x=as.character(x)
  x= sapply(x, function(x){
    switch(x,
           A="T",
           T="A",
           C="G",
           G="C")})
  return(x)
} # flip a SNP

is.ATCG = function(x,y){
  return((x=="A" & y=="T") |(x=="T" & y =="A") | (x=="C" & y=="G") | (x=="G" & y=="C"))
} # c

nonATCG = function(x,y){
  same = (x$A1==y$A1 & x$A2==y$A2) |  (x$A2==y$A1 & x$A1==y$A2)
  flip = (x$A1==flip(y$A1) & x$A2==flip(y$A2)) | (x$A2==flip(y$A1) & x$A1==flip(y$A2))
  monomorph = (x$A1 == 0 & x$A2==y$A2) |   (y$A1 == 0 & y$A2==x$A2)
  flip_monomorph = (x$A1 == 0 & x$A2==flip(y$A2)) |   (y$A1 == 0 & y$A2==flip(x$A2) )
  flip = flip  | flip_monomorph
  exclude = !(any(same | flip | monomorph | flip_monomorph))
  x$action[same]="no.action"
  x$action[flip]="flipped"
  x$action[exclude]="excluded"
  x$type="nonATCG"
  x$action = factor(x$action, levels=c("no.action", "flipped","excluded"))
  return(cbind(x,y))
} # check if a non ATCG has to be flipped



# TODO
ATCG = function(x,y){
  # x = file,  y =reference
  diff = abs(as.numeric(y$MAF) - as.numeric(x$MAF))
  same = (x$A1==y$A1 & x$A2==y$A2)
  flip = (x$A1==y$A2 & x$A2==y$A1)   
  exclude =  (((x$MAF > 0.4)  |(y$MAF > 0.4 )))
  # SNP is an ATCG SNP with MAF > 0.4 cannot be assigned
  x$action[same]="no.action"
  x$action[flip]="flipped"
  x$action[exclude]="excluded"
  x$type="ATCG"
  x$action = factor(x$action, levels=c("no.action", "flipped","excluded"))
  return(cbind(x,y))
  return(cbind(x,y))
} #

############################
# DO ANALYSIS
############################


print("Searching for nonATCG (alleles different) to exclude and ATCG to exclude (alleles different but AF similar)")
print(head(reference))
ind = is.ATCG(reference[,2], reference[,3])
print(head(nonATCG(tmp[!ind,],reference[!ind,])))

if (any(ind)){
  print(head(ATCG( tmp[ind,], reference[ind,])))
  out = rbind(nonATCG(tmp[!ind,],reference[!ind,]),
              ATCG( tmp[ind,], reference[ind,]))
} else {
  out = nonATCG(tmp[!ind,],reference[!ind,])
}



write.table(out$id[grep("excluded", out$action)], "exclude.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(out$id[grep("no.action|flipped", out$action)], "extract.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(out$id[grep("flipped", out$action)],    "flip.txt", sep="\t", row.names=F, col.names=F, quote=F)


#######################################
# Exclude & FLIP
#######################################
cmd = paste("plink --bfile",
            params.rootname, 
            "--extract extract.txt",
            "--make-bed --out",
            paste0(params.rootname,".refchecked"))
system(cmd)
cmd = paste("plink --bfile",
            paste0(params.rootname,".refchecked"), 
            "--flip flip.txt",
            "--make-bed --out",
            paste0(params.rootname,".refchecked"))
system(cmd)

############################
# GET FREQUENCY IN NEW DATA
###########################
cmd = paste("plink --bfile", paste0(params.rootname,".refchecked"),"--freq --out tmp")
system(cmd)

tmp = cbind(read.table(paste0(params.rootname,".refchecked.bim"),                     
                       col.names=c("chr","id","","pos","A1","A2"), colClasses="character"),
            read.table("tmp.frq", h=T, colClasses="character"))

############################
# MATCH
###########################
reference = reference[reference$pos%in%tmp$pos,]
tmp       = tmp[match(reference$pos, tmp$pos),]

ind = is.ATCG(reference[,2], reference[,3])
if (any(ind)){
  out.plot = rbind(nonATCG(tmp[!ind,],reference[!ind,]),
            ATCG( tmp[ind,], reference[ind,]))
} else {
  out.plot = rbind(nonATCG(tmp[!ind,],reference[!ind,]))
}
save(out.plot,reference, tmp, file="tmp.RData")

#############################
# Plot data
#############################


out = data.frame(out)
out.plot = data.frame(out.plot)

p = ggplot(out,aes(x=as.numeric(MAF), y=as.numeric(MAF.1), col=action,pch=type)) +
  geom_point(size=2) + xlab("AF reference") + ylab("AF input data") + theme_bw() + 
  scale_color_manual(values = c("darkgreen", "blue", "red"),name="action", drop=F) +
  scale_shape_manual(values =c(16,17), name="type", drop=F) + ggtitle("BEFORE") + theme(legend.position="none")

q = ggplot(out.plot,aes(x=as.numeric(MAF), y=as.numeric(MAF.1), col=action,pch=type)) +
  geom_point(size=2) + xlab("AF reference") + ylab("AF input data") + theme_bw() + 
  scale_color_manual(values = c("darkgreen", "blue", "red"),name="action", drop=F) +
  scale_shape_manual(values =c(16,17), name="type", drop=F) + ggtitle("AFTER")

get_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
leg = get_legend(q)

q = q +  theme(legend.position="none")

################################################
# CALCULATE OVERLAP BETWEEN STUDY AND REFERENCE
################################################


bim = read.table(paste0(params.rootname,".refchecked.bim"), h=F, 
                 col.names=c("chr","id","","pos","A1","A2"), colClasses="character")
pos = lapply(model, function(x){pos = x$snp.position;
                                       unlist(lapply(x$classifiers,function(x){x=pos[x$snpidx];x=sum(x%in%bim$pos)/length(x) }))})
loc=unique(unlist(lapply(model, function(x){pos = x$snp.position})))
pos = cbind(rep(names(pos), times=unlist(lapply(pos, length))),do.call(c, pos))


pos = data.frame(pos)

save(pos, loc, bim, model, file="tmp.RData")
r= ggplot(aes(x=X1, y=as.numeric(as.matrix(X2))), data= pos) + geom_boxplot() +
  ggtitle(paste("Study contains", sum(bim$pos%in%loc),"positions of", length(loc), "present in the used model.")) + labs(x="Locus",y="%-available SNPS used in panel") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0,1))
          
png("refchecked_A.png")
print(grid.arrange(p,q, leg,ncol=3))
dev.off()

png("refchecked_B.png")
print(r)
dev.off()



################################################
# WRITE OUTPUT FOR REPORT
################################################

write.table(out, 
            paste0(params.rootname, ".summary.refchecked.txt",sep=""), 
            col.names=T, row.names=F, quote=F, sep="\t")


summary = table(out$action, out$type)

to.output = function(x){
  a = colnames(x)
  b = rownames(x)
  x = apply(x, 2, as.numeric)
  x = cbind(b,x)
  colnames(x) = c("",a)
  return(x)
}
summary = to.output(summary)


write.table(summary,file="report.refchecked.txt",quote=F,sep=",", col.names=T, row.names=T)
