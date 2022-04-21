options(stringsAsFactors=F)
tmp = read.table(file.path(dir, "hla_nom_g.txt"), h=F)

# FORMAT LOCUS ALLELES GROUP
tmp = do.call(rbind,lapply(strsplit(tmp$V1,";"), 
                           function(x){if(length(x)==2){
                             x = c(x[1],x[2], paste0(x[2],"G")); 
}
                             return(x)}))
# SPLIT ALLELES OF A GROUP
out = apply(tmp,1, function(x){
y = unlist(strsplit(x[2], "/"));
x = cbind(x[1], y, x[3]); return(x)})

out = do.call(rbind, out)
out[,2] = sapply(out[,2], function(x){x= unlist(strsplit(x,":")); 
if(length(x)>2){x= paste(x[1], x[2], x[3], sep=":") 
}else{x= paste(x[1],x[2], sep=":")}
return(x)})
out = unique(out)

twofield = sapply(out[,2], function(x){x= unlist(strsplit(x,":")); 
x= paste(x[1], x[2], sep=":") 
return(x)})

twofieldG = sapply(out[,3], function(x){x= unlist(strsplit(x,":")); 
x= paste(x[1], x[2], sep=":") 
return(x)})
twofieldG =gsub("GG", "G",paste0(twofieldG, "G"))


out = cbind(out,twofield, twofieldG)
colnames(out) = c("locus","3field","3Gfield", "2field", "2Gfield") 
out[,2] = gsub("[QLN]", "", out[,2])
out = data.frame(out)
out[,1] = gsub("\\*", "", out[,1])
save(out, file="g_groups.RData")
