#### lm of per gene nucleotide diversity and outliers
library(ggplot2)
library(magrittr)
library(dplyr)
library(gplots)
library(ggthemes)
library(reshape2)
library(tidyverse)
library(ggsci)
library(plyr)

#### Distribution of nucleotide diversity in sliding windows
genes.inter.intra.df <- read.table("190214_interIntra_perGenePi.txt", 
                                   header = T, 
                                   sep = "\t",
                                   na.strings = NA,
                                   stringsAsFactors = F)

genes.inter.intra.df$Patient <- factor(genes.inter.intra.df$Patient, 
                                       levels=c("Patient 2","Patient 3","Patient 4","Patient 5",
                                                "Patient 7","Patient 8","Patient 9","Patient 10",
                                                "Patient 11","Patient 14","Patient 21","Patient 22",
                                                "Patient 23","Patient 35-Sample 1","Patient 35-Sample 2","Patient 35-Sample 3"))

genes.inter.intra.df$Culture <- as.numeric(genes.inter.intra.df$Culture)
genes.inter.intra.df$Sputum <- as.numeric(genes.inter.intra.df$Sputum)

patients <- c(unique(as.character(genes.inter.intra.df$Patient)))

log.df <- mutate(genes.inter.intra.df, log.pi.s = log10(Sputum), log.pi.c = log10(Culture))
log.df <- log.df[apply(log.df[,c("log.pi.c","log.pi.s")], 1, function(x) all(is.finite(x))),]

lm.df <- NULL

pvals <- c()

patients <- c(unique(as.character(genes.inter.intra.df$Patient)))

for (x in patients) {
  log.df.patient <- log.df[which(as.character(log.df$Patient) == x),]
  fit <- lm(log.pi.c ~ log.pi.s, data=log.df.patient)
  pval <- anova(fit)$'Pr(>F)'[1]
  cooksd <- cooks.distance(fit)
  cooks.df <- as.data.frame(cooksd)
  cutoff <- 4*(sum(cooks.df$cooksd)/length(cooks.df$cooksd))
  cooks.outliers <- which(cooks.df > cutoff)
  log.df.patient$Cooks <- cooks.df$cooksd
  log.df.patient$outlier <- FALSE
  log.df.patient$outlier[cooks.outliers] <- TRUE
  pvals <- c(pvals, pval)
  lm.df <- rbind(lm.df, log.df.patient)
}

pvals

n.occur.lm <- data.frame(table(lm.df$Gene))
n.occur.lm <- n.occur.lm[n.occur.lm$Freq > 1,]
n.occur.lm <- n.occur.lm[order(n.occur.lm$Var1, decreasing = T),]