### z score for fold change in nucleotide diversity
library(ggplot2)
library(magrittr)
library(dplyr)
library(gplots)
library(viridis)
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

### vector of patients
patients <- c(unique(as.character(genes.inter.intra.df$Patient)))

genes.inter.intra.df <- as.tibble(genes.inter.intra.df)

### calculate fold change and remove 0s,Inf and NAs
genes.inter.intra.df$fold.change <- genes.inter.intra.df$Sputum/genes.inter.intra.df$Culture
genes.inter.intra.df.fin <- genes.inter.intra.df[apply(genes.inter.intra.df[,c("fold.change")], 1, function(x) all(is.finite(x))),]
genes.inter.intra.df.fin <- genes.inter.intra.df.fin[which(genes.inter.intra.df.fin$fold.change > 0),]

sig.fold <- NULL

### calculate z score and fdr adjusted p-value per gene in each patient; bind to data frame 
for (x in patients) {
  df.patient <- genes.inter.intra.df.fin[which(as.character(genes.inter.intra.df.fin$Patient) == x),]
  z <- scale(df.patient$fold.change, center = TRUE, scale = TRUE)
  pval = pnorm(-abs(z))
  pval = as.vector(2 * pval)
  pval.adj <- p.adjust(pval, method = "fdr", n = length(pval))
  pval.adj <- as.vector(pval.adj)
  if (any(pval.adj < 0.05, na.rm = T)) {
    sig.val <- data.frame(gene = df.patient$Gene[which(pval.adj < 0.05)],
                          fold.change = df.patient$fold.change[which(pval.adj < 0.05)],
                          pval.adj = pval.adj[which(pval.adj < 0.05)])
    sig.val$patient <- x
    sig.fold <- rbind(sig.fold, sig.val)
  }
}

### table of outlier genes and their frequency, ordered from high to low frequency
sig.fold$gene <- as.character(sig.fold$gene)
n.occur.z <- data.frame(table(sig.fold$gene))
n.occur.z <- n.occur.z[order(n.occur.z$Var1, decreasing = TRUE),]
