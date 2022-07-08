##############################################
# PLOT SNPS/HLA alleles
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# May 2022
##############################################
################################
# SETTINGS
################################
options(stringsAsFactors=F)

args = commandArgs(T)

stats = args[1] # single association file
loc = args[2] # location file (e.g. location_HLA_genes_hg19.txt)


if(length(args)>2){
  meta = args[3] # file with meta-analysis
  info = read.table(meta, h=T, sep="\t", fill =T)
  colnames(info)[colnames(info)%in%c("BP","POS", "pos")]="POS"
  colnames(info)[colnames(info)%in%c("id","rsid", "ID","RSID","SNP")]="SNPID"
  colnames(info)[colnames(info)%in%c("PVALUE_RE2","marker.p")]="P"
  colnames(info)[colnames(info)%in%c("I2")]="I_SQUARE"
}

## FORMAT LOCATION FILES
loc = read.table(loc,h=T,sep="\t",fill=T)
loc = loc[loc$locus%in%c("A","B","C","DPA1","DPB1","DRB1","DQA1","DQB1","DRB3","DRB4","DRB5"),]
loc  = data.frame(locus=names(tapply(loc$start, gsub("[AB]1","", loc$locus),invisible)),
                  start=tapply(loc$start, gsub("[AB]1","", loc$locus), min),
                  end=tapply(loc$end, gsub("[AB]1","", loc$locus), max))
loc$start=as.numeric(loc$start)
loc$end=as.numeric(loc$end)
loc$mid=(loc$start+ (loc$end-loc$start)/2)/10^6

loc_new= data.frame(rbind(loc,
                          c(locus="DR/DQ", start=min(loc[loc$locus%in%c("DR", "DQ","DRB3","DRB4","DRB5"),2]),
                            end = max(loc[loc$locus%in%c("DR", "DQ","DRB3","DRB4","DRB5"),3]),
                            mid = mean(loc[loc$locus%in%c("DR", "DQ","DRB3","DRB4","DRB5"),4])),
                          c(locus="C/B", start=min(loc[loc$locus%in%c("C","B"),2]),
                            end = max(loc[loc$locus%in%c("C","B"),3]),
                            mid = mean(loc[loc$locus%in%c("C","B"),4]))))
loc_new = loc_new[!loc_new$locus%in%c("DR","DQ","C","B","DRB3","DRB4","DRB5"),]
loc = loc_new
loc[,2:4] = apply(loc[,2:4],2,as.numeric)

########################################
# PLOT
########################################
stats = read.table(stats, h=T)
colnames(stats)[colnames(stats)%in%c("BP","POS", "pos")]="POS"
colnames(stats)[colnames(stats)%in%c("id","ID", "rsid","RSID","SNP")]="SNPID"
colnames(stats)[colnames(stats)%in%c("PVALUE_RE2","marker.p")]="P"
colnames(stats)[colnames(stats)%in%c("AF","FRQ","MAF")]="AF_ALL"
stats=stats[!is.na(stats$P),]
print(head(stats))
print(quantile(-log10(stats$P)))


if("AF_ALL"%in%colnames(stats)){
  ## Get only more frequent data
  stats = stats[stats$AF_ALL > 0.01 & stats$AF_ALL < 0.99,]
}


if(length(args)>2){
  info = info[match(stats$SNPID, info$SNPID),]
  info$POS=stats$POS
}

cex=2
########################################
# PLOT
########################################

png(paste0(args[1],".png"),width=1180,height=420)
par(lwd=3)
par(fig= c(0,1,0.85,1), mar=c(0,5.1,4.1,5.5) )
plot(1, xlab="", xlim=c(29,34), col="transparent", xaxt="n", ylab="", yaxt="n")
loc$col="black"

# Make box with allele location
for(i in 1:nrow(loc)){
  abline(v=loc$start[i]/10^6, las =loc$class[i], col=loc$col[i])
  abline(v=loc$end[i]/10^6, las =loc$class[i], col = loc$col[i])
}  

axis(3,labels=loc_new$locus, at=loc_new$mid,las=1, cex.axis=cex, lwd=3)

## Plot SNPS in grey
par(fig=c(0,1,0,0.85),mar=c(5.1,5.1,0,5.5), new=T)
type=rep("snp", nrow(stats))
print(head(stats))
type[grep("\\*|HLA",stats$SNPID)]="hla"
type[grep("prot",stats$SNPID)]="prot"
print(table(type))
type=factor(type,c("hla","snp","prot"))
col = c("gold","grey","red")[factor(type)]
plot(stats$POS/10^6, -log10(stats$P), col=col,  ylab="-log10(P)", xlab="chromosome 6 [Mb]", xlim=c(29,34), 
     cex.axis=cex,cex.lab=cex, pch=16, cex=1.5, ylim=c(0, max(-log10(stats$P))+10))

## Plot HLA alleles in golden
if(any(type=="hla")){
  points(cbind(stats$POS/10^6, -log10(stats$P))[type=="hla",],
         col = "gold", pch=16, cex=1.5)
  axis(side=1, cex=cex, cex.axis=cex, lwd=3)
  axis(side=2, cex=cex, cex.axis=cex, lwd=3)
  
  legend("topright",c("SNPs","HLA alleles"),col=c("grey","gold"),pch=16,cex=cex)
}

# ## Plot amino acids in red
# if(any(type=="prot")){
# points(cbind(stats$POS/10^6, -log10(stats$P))[type=="prot",],
#        col = "red", pch=16, cex=1.5)
# axis(side=1, cex=cex, cex.axis=cex, lwd=3)
# axis(side=2, cex=cex, cex.axis=cex, lwd=3)
# 
# legend("topright",c("SNPs","HLA alleles","prot"),col=c("grey","gold","red"),pch=16,cex=cex)
# }

## Plot lines of significance
abline(h=-log10(5*10^-8), lty=2, col="black")
abline(h=-log10(10^-5), lty=2, col="darkgrey")

if(length(args)>2){

  par(fig= c(0,1,0,0.85),mar=c(5.1,5.1,0,5.5),new=T)
  
  # Plot smooth line accross positions 29-34Mb with smooth line
  mean = tapply(info$I_SQUARE, cut(info$POS/10^6, seq(29,34,0.005)), function(x){mean(x,na.rm=T)}) 
  names(mean) = seq(29,33.995,0.005)
  mean = mean[!is.na(mean)]
  lo = smooth.spline(mean~names(mean))
  plot("", ylim=c(0, 100),
         type="l", xlab="", ylab="", xlim=c(29,34), xaxt="n", yaxt="n",cex.axis=cex, cex=cex)
  lines(lo, lwd=3)
  axis(side=4, ylim=c(0,100),cex=cex, cex.axis=cex, lwd=3, at = pretty(range(mean)))
  mtext("I_SQUARE", side=4, line=4, cex=cex)

}

dev.off()


