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

### read in pi per gene
genes.inter.intra.df <- read.table("interIntra_perGenePi.txt", 
                                   header = T, 
                                   sep = "\t",
                                   na.strings = NA,
                                   stringsAsFactors = F)

### patient as factor
genes.inter.intra.df$Patient <- factor(genes.inter.intra.df$Patient, 
                                       levels=c("Patient 2","Patient 3","Patient 4","Patient 5",
                                                "Patient 7","Patient 8","Patient 9","Patient 10",
                                                "Patient 11","Patient 14","Patient 21","Patient 22",
                                                "Patient 23","Patient 35-Sample 1","Patient 35-Sample 2","Patient 35-Sample 3"))
### pi as numeric
genes.inter.intra.df$Culture <- as.numeric(genes.inter.intra.df$Culture)
genes.inter.intra.df$Sputum <- as.numeric(genes.inter.intra.df$Sputum)

### list of patients
patients <- c(unique(as.character(genes.inter.intra.df$Patient)))


### log transform pi values for sputum and culture, remove Inf and NA
log.df <- mutate(genes.inter.intra.df, log.pi.s = log10(Sputum), log.pi.c = log10(Culture))
log.df <- log.df[apply(log.df[,c("log.pi.c","log.pi.s")], 1, function(x) all(is.finite(x))),]

lm.df <- NULL

pvals <- c()

### fit to lm, calculate F-test p-value for each model and Cook's distance per gene; bind to data frame
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

### p-values for F-test of overall significance 
pvals

### data frame of outliers only
outlier.genes <- d[which(d$outlier==TRUE),]

### table of outlier genes and their frequency, ordered from high to low frequency
n.occur.lm <- data.frame(table(lm.df$Gene))
n.occur.lm <- n.occur.lm[order(n.occur.lm$Var1, decreasing = T),]