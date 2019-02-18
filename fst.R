fet.df <- read.table("pvalsFstFishersExact.txt", 
                                   header = T, 
                                   sep = "\t",
                                   na.strings = NA,
                                   stringsAsFactors = F)


fet.df <- as.tibble(fet.df)
colnames(fet.df) <- c("Patient","Gene","neglog10pval")

sig.fet.df <- NULL

patients <- c(unique(as.character(fet.df$Patient)))

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

sig.fet.df$Gene <- as.character(sig.fet.df$Gene) 
n.occur.fst <- data.frame(table(sig.fet.df$Gene))

n.occur.fst$Freq <- as.numeric(n.occur.fst$Freq) 
n.occur.fst <- n.occur.fst[n.occur.fst$Freq > 1,]
n.occur.fst <- n.occur.fst[order(n.occur.fst$Freq, decreasing = T),]
