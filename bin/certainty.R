#!/usr/bin/env Rscript

#############################################
# Calculate phasing uncertainty
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################

################################
# Settings
################################
library(parallel)


args = commandArgs(T)

name = args[1]
args = args[-1]

################################
#  Main
################################

tmp = mclapply(as.list(args),function(x){
  print(paste("Reading", x))
  return(read.table(x, h=F))},mc.cores=8)

samples = paste(name,".S1.sample", sep="")
samples = read.table(samples, h=T)[-1,]
id = paste(samples$ID_1, samples$ID_2,sep="_")

tmp = do.call(rbind,tmp)
tmp = tmp[tmp[,6]%in%c(0,1),]

certainty= apply(tmp[,-(1:5)],2,function(x){return(tapply(as.numeric(as.matrix(x)),tmp$V2, mean))})

certainty= apply(certainty, 2, function(x){x[x<0.5]= 1- x[x< 0.5]; return(x)})
pos = unlist(strsplit(rownames(certainty), ":", ))[c(F,T)]

colnames(certainty)=rep(id, each=2)



print(paste("Plotting certainty"))
pdf(paste(name,".certainty.pdf",sep=""))
hist(apply(certainty, 2, mean), breaks="FD", ylab="frequency", xlab="mean certainty sample",  main="")
hist(apply(certainty, 2, median), breaks="FD", ylab="frequency", xlab="median certainty sample",  main="")
plot(as.numeric(pos)/10^6, apply(certainty, 1, mean),ylab="mean certainty",xlab="position chr6 in Mb")
plot(as.numeric(pos)/10^6, apply(certainty, 1, median),ylab="median certainty",xlab="position chr6 in Mb")
dev.off()


med_sample = apply(certainty, 2, median)
med_snp = apply(certainty, 1, median)

if(nrow(samples)>10){
  out = med_snp<0.8
}else{
  out = med_snp<0
}


print(paste("Writing files..."))
write.table(names(med_snp)[out], paste(name,".certainty.exclude",sep=""),sep="\t", row.names=F, col.names=F)

med_sample = apply(certainty[!out, ], 1,mean)
med_snp = med_snp[!out]

write.table(certainty[!out,], paste(name,".certainty.all",sep=""),sep="\t", row.names=T, col.names=T)
write.table(med_snp, paste(name,".certainty.snp",sep=""),sep="\t", row.names=T, col.names=T)
write.table(med_sample, paste(name,".certainty.sample",sep=""),sep="\t", row.names=T, col.names=T)
