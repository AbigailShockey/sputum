######### Overlap in windows of extreme nucleotide diversity
library(ggplot2)
library(magrittr)
library(dplyr)
library(gplots)
library(ggthemes)
library(reshape2)
library(tidyverse)
library(ggsci)
library(plyr)

### Read in slinding-window tables (long format)
inter.intra.df.l <- read.table("/Users/abbas/Documents/analysis_05.07.18/190212_interIntra_windowsLong.txt", 
                               header = T, 
                               sep = "\t",
                               na.strings = NA,
                               stringsAsFactors = F)

inter.intra.df.l <- as.tibble(inter.intra.df.l)

### stat as character
inter.intra.df.l$Stat <- as.character(inter.intra.df.l$Stat)

### patient as factor
inter.intra.df.l$Patient <- factor(inter.intra.df.l$Patient, 
                                   levels=c("Patient 2","Patient 3","Patient 4","Patient 5",
                                            "Patient 7","Patient 8","Patient 9","Patient 10",
                                            "Patient 11","Patient 14","Patient 21","Patient 22",
                                            "Patient 23","Patient 35-Sample 1","Patient 35-Sample 2","Patient 35-Sample 3"))

### list of patients
patients <- c(unique(as.character(inter.intra.df.l$Patient)))

### inter-host samples
mtb.patients <- patients[1:13]

### intra-host samples
mbovis.patients <- patients[14:16]


### finite values only
inter.intra.df.l <- inter.intra.df.l[apply(inter.intra.df.l[,c("sput.val", "cult.val")], 1, function(x) all(is.finite(x))),]

### inter-host analysis, sputum
sig.win.sput <- NULL

### calculate z-score for nucleotide diversity in each window per patient
for (x in mtb.patients) {
  df.patient <- inter.intra.df.l[which(as.character(inter.intra.df.l$Patient) == x & inter.intra.df.l$Stat == "pi"),]
  z <- scale(df.patient$sput.val, center = TRUE, scale = TRUE)
  pval = pnorm(-abs(z))
  pval.adj <- p.adjust(pval, method = "fdr", n = length(pval))
  pval.adj <- as.vector(pval.adj)
  if (any(pval.adj < 0.05, na.rm = T)) {
    sig.val <- data.frame(window = df.patient$window[which(pval.adj < 0.05)], sput.val = df.patient$sput.val[which(pval.adj < 0.05)],pval.adj = pval.adj[which(pval.adj < 0.05)])
    sig.val$patient <- x
    sig.win.sput <- rbind(sig.win.sput, sig.val)
  }
}

n.occur.s <- data.frame(table(as.numeric(sig.win.sput$window)))
n.occur.sig.s <- n.occur.s[n.occur.s$Freq > 1,]
n.occur.sig.s <- n.occur.sig.s[order(n.occur.sig.s$Var1),]

### inter-host analysis, culture
sig.win.cult <- NULL

for (x in mtb.patients) {
  df.patient <- inter.intra.df.l[which(as.character(inter.intra.df.l$Patient) == x & inter.intra.df.l$Stat == "pi"),]
  z <- scale(df.patient$cult.val, center = TRUE, scale = TRUE)
  pval = pnorm(-abs(z))
  pval.adj <- p.adjust(pval, method = "fdr", n = length(pval))
  pval.adj <- as.vector(pval.adj)
  if (any(pval.adj < 0.05, na.rm = T)) {
    sig.val <- data.frame(window = df.patient$window[which(pval.adj < 0.05)],
                          cult.val = df.patient$cult.val[which(pval.adj < 0.05)],
                          pval.adj = pval.adj[which(pval.adj < 0.05)])
    sig.val$patient <- x
    sig.win.cult <- rbind(sig.win.cult, sig.val)
  }
}

n.occur.c <- data.frame(table(as.numeric(sig.win.cult$window)))
n.occur.sig.c <- n.occur.c[n.occur.c$Freq > 1,]
n.occur.sig.c <- n.occur.sig.c[order(n.occur.sig.c$Var1),]
length(n.occur.c$Var1[n.occur.c$Freq > 1])

### intra-host analysis, sputum
n.occur.s <- NULL
sig.win.sput <- NULL

for (x in mbovis.patients) {
  df.patient <- inter.intra.df.l[which(as.character(inter.intra.df.l$Patient) == x & inter.intra.df.l$Stat == "pi"),]
  z <- scale(df.patient$sput.val, center = TRUE, scale = TRUE)
  pval = pnorm(-abs(z))
  pval.adj <- p.adjust(pval, method = "fdr", n = length(pval))
  pval.adj <- as.vector(pval.adj)
  if (any(pval.adj < 0.05, na.rm = T)) {
    sig.val <- data.frame(window = df.patient$window[which(pval.adj < 0.05)], sput.val = df.patient$sput.val[which(pval.adj < 0.05)],pval.adj = pval.adj[which(pval.adj < 0.05)])
    sig.val$patient <- x
    sig.win.sput <- rbind(sig.win.sput, sig.val)
  }
}

n.occur.s <- data.frame(table(as.character(sig.win.sput$window)))
n.occur.sig.s <- n.occur.s[n.occur.s$Freq > 1,]
n.occur.sig.s <- n.occur.sig[order(n.occur.sig.s$Var1),]

### intra-host analysis, culture
n.occur.c <- NULL
sig.win.cult <- NULL

for (x in mbovis.patients) {
  df.patient <- log.df[which(as.character(log.df$Patient) == x & log.df$Stat == "pi"),]
  z <- scale(df.patient$cult.val, center = TRUE, scale = TRUE)
  pval = pnorm(-abs(z))
  pval.adj <- p.adjust(pval, method = "fdr", n = length(pval))
  pval.adj <- as.vector(pval.adj)
  if (any(pval.adj < 0.05, na.rm = T)) {
    sig.val <- data.frame(window = df.patient$window[which(pval.adj < 0.05)], cult.val = df.patient$cult.val[which(pval.adj < 0.05)],pval.adj = pval.adj[which(pval.adj < 0.05)])
    sig.val$patient <- x
    sig.win.cult <- rbind(sig.win.cult, sig.val)
  }
}

n.occur.c <- data.frame(table(as.character(sig.win.cult$window)))
n.occur.sig.c <- n.occur.c[n.occur.c$Freq > 1,]
n.occur.sig.c <- n.occur.sig.c[order(n.occur.sig.c$Var1),]
