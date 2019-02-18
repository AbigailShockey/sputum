####Fst analsyes for brown data

file_list <- list.files(path="C:/Users/abbas/Documents/analysis_05.07.18/genesFst/",pattern = ".fet")

fet.df <- NULL

for (file in file_list){
  # if the merged dataset does exist, append to it
  if (exists("fet.df")){
    temp_dataset <-read.table(paste("C:/Users/abbas/Documents/analysis_05.07.18/genesFst/",file,sep = ""), header=F, sep="\t")
    temp_dataset$Patient <- file
    temp_dataset$Patient <- gsub("_genes.fet","",temp_dataset$Patient)
    temp_dataset$Patient <- gsub("patient","",temp_dataset$Patient)
    temp_dataset$V6 <- gsub("1:2=","",temp_dataset$V6)
    temp_dataset$V6 <- gsub("na",NA,temp_dataset$V6)
    temp_dataset$V6 <- as.numeric(temp_dataset$V6)
    temp_dataset$Patient <- as.character(temp_dataset$Patient)
    fet.df <-rbind(fet.df, temp_dataset)
    rm(temp_dataset)
  }
}

fet.df <- as.tibble(fet.df)

sig.fet.df <- NULL

patients <- c(unique(as.character(fet.df$Patient)))

for (x in patients) {
  p.fet.df <- fet.df[which(as.character(fet.df$Patient) == x),]
  correction.l <- length(p.fet.df$V1)
  p.fet.df <- p.fet.df[-which(is.na(p.fet.df$V6)),]
  p.fet.df <- p.fet.df[-which(p.fet.df$V6 == 0),]
  p.fet.df$V7 <- as.numeric(10**(-1*p.fet.df$V6))
  pval.adj <- p.adjust(p.fet.df$V7, method = "fdr", n = correction.l)
  if (any(pval.adj < 0.05, na.rm = T)) {
    sig.p.fet.df <- p.fet.df[which(pval.adj < 0.05),]
    sig.fet.df <- rbind(sig.fet.df, sig.p.fet.df)
  }
}

sig.fet.df$V1 <- as.character(sig.fet.df$V1) 
n.occur.fst <- data.frame(table(sig.fet.df$V1))

n.occur.fst$Freq <- as.numeric(n.occur.fst$Freq) 
n.occur.fst <- n.occur.fst[n.occur.fst$Freq > 1,]
n.occur.fst <- n.occur.fst[order(n.occur.fst$Freq, decreasing = T),]
