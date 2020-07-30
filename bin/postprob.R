#!/usr/bin/env Rscript

##############################################
# Plot posterior probability of imputation
# author: Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################


################################
# Settings
################################

args = commandArgs(T)


library(ggplot2)
library(reshape)


png("postprob.png")

################################
#Posterior probability
################################

load(args[1])

post = sapply(pred,function(x) x$value[,4])
post = melt(post)

colnames(post) <- c("sample","locus","postprob")

ggplot(aes(x=postprob,color=locus), data=post) + 
  stat_ecdf(geom="step")+
  xlab("post prob") + ylab("proportion of posterior propability <= x") +
  scale_y_continuous(breaks=seq(0,1,0.2),expand = c(0,0)) + scale_x_continuous(breaks=seq(0,1,0.2),expand = c(0,0)) +
  ylim(0,1)+ theme_bw()

dev.off()

