##### Pairwise differences in genome-wide means

library(ggplot2)
library(magrittr)
library(dplyr)
library(gplots)
library(ggthemes)
library(reshape2)
library(tidyverse)
library(ggsci)
library(plyr)

inter.intra.df.l <- read.table("/Users/abbas/Documents/analysis_05.07.18/190212_interIntra_windowsLong.txt", 
                               header = T, 
                               sep = "\t",
                               na.strings = NA,
                               stringsAsFactors = F)

inter.intra.df.l$Patient <- factor(inter.intra.df.l$Patient, 
                                   levels=c("Patient 2","Patient 3","Patient 4","Patient 5",
                                            "Patient 7","Patient 8","Patient 9","Patient 10",
                                            "Patient 11","Patient 14","Patient 21","Patient 22",
                                            "Patient 23","Patient 35-Sample 1","Patient 35-Sample 2","Patient 35-Sample 3"))

inter.intra.df.l <- as.tibble(inter.intra.df.l)

group.means.s <- ddply(inter.intra.df.l[which(inter.intra.df.l$Stat == "pi"),],c("Patient"),summarise,mean=mean(sput.val, na.rm = T))
group.means.c <- ddply(inter.intra.df.l[which(inter.intra.df.l$Stat == "pi"),],c("Patient"),summarise,mean=mean(cult.val, na.rm = T))
group.means <- data.frame(Patient = group.means.s$Patient, Sputum = group.means.s$mean, Culture = group.means.c$mean)

mtb.s <- as.data.frame(t(combn(group.means.s$mean[1:13], 2)))
mtb.c <- as.data.frame(t(combn(group.means.c$mean[1:13], 2)))
mtb.s.c <- rbind(mtb.s,mtb.c)
mtb.s.c$diff <- abs(mtb.s.c$V1 - mtb.s.c$V2)
mtb.s.c$Sample <- "Inter-patient"

mb.s <- as.data.frame(t(combn(group.means.s$mean[14:16], 2)))
mb.c <- as.data.frame(t(combn(group.means.c$mean[14:16], 2)))
mb.s.c <- rbind(mb.s,mb.c)
mb.s.c$diff <- abs(mb.s.c$V1 - mb.s.c$V2)
mb.s.c$Sample <- "Intra-patient"

mb.mtb.diffs <- rbind(mtb.s.c,mb.s.c)

