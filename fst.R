library(ggplot2)
library(magrittr)
library(dplyr)
library(gplots)
library(ggthemes)
library(reshape2)
library(tidyverse)
library(ggsci)
library(plyr)

### read in p-values from Fisher's exact test per gene
fet.df <- read.table("pvalsFstFishersExact.txt", 
                                   header = T, 
                                   sep = "\t",
                                   na.strings = NA,
                                   stringsAsFactors = F)


fet.df <- as.tibble(fet.df)

### set col names
colnames(fet.df) <- c("Patient","Gene","neglog10pval")

sig.fet.df <- NULL

### vector of patients
patients <- c(unique(as.character(fet.df$Patient)))

### calculate p-value from -log10() transformed value and fdr adjusted p-value per gene in each patient
for (x in patients) {
  p.fet.df <- fet.df[which(as.character(fet.df$Patient) == x),]
  correction.l <- length(p.fet.df$Gene)
  p.fet.df <- p.fet.df[-which(is.na(p.fet.df$neglog10pval)),]
  p.fet.df <- p.fet.df[-which(p.fet.df$neglog10pval == 0),]
  p.fet.df$pval <- as.numeric(10**(-1*p.fet.df$neglog10pval))
  pval.adj <- p.adjust(p.fet.df$pval, method = "fdr", n = correction.l)
  if (any(pval.adj < 0.05, na.rm = T)) {
    sig.p.fet.df <- p.fet.df[which(pval.adj < 0.05),]
    sig.fet.df <- rbind(sig.fet.df, sig.p.fet.df)
  }
}

### table of outlier genes and their frequency, ordered from high to low frequency
sig.fet.df$Gene <- as.character(sig.fet.df$Gene) 
n.occur.fst <- data.frame(table(sig.fet.df$Gene))
n.occur.fst <- n.occur.fst[order(n.occur.fst$Freq, decreasing = T),]
