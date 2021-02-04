# data frame with all the characteristics
clinical_feat <- c("n patients", "p miRNA", "p mRNA", "Censoring rate", "Survival - 3 years", 
                   "C-index - EN", "C-index - RF", "Library size - miRNA (thousands of reads)", 
                   "Library size - mRNA (thousands of reads)")
data_ch <- data.frame(matrix(ncol = length(clinical_feat), nrow = length(cancers_all)))
row.names(data_ch) <- cancers_all
colnames(data_ch) <- clinical_feat

med_seq_depth <- data.frame(matrix(ncol = 2, nrow = length(cancers_all)))
row.names(med_seq_depth) <- cancers_all
colnames(med_seq_depth) <- c("miRNA", "mRNA")

lib_size_mRNA = lib_size_miRNA <- list()

load(file = "data_fit/pred_ch_cancers.RData")

for(cancer in cancers_all){

  print(paste("*** start learning for", cancer, "***"))
  
  # load the data ---
  source(file = "load_data/load_data_final.R")
  
  # librayr size (number of reads for each patients)
  lib_size_miRNA[[cancer]] <- apply(count, 1, sum)
  lib_size_mRNA[[cancer]] <- apply(count_mRNA, 1, sum)

  # compute the sequencing depth
  med_seq_depth[cancer, "miRNA"] <- median(lib_size_miRNA[[cancer]])
  med_seq_depth[cancer, "mRNA"] <- median(lib_size_mRNA[[cancer]])
  
  data_ch[cancer, "Library size - miRNA (thousands of reads)"] <- signif(median(lib_size_miRNA[[cancer]]) / 10^6, 1)
  data_ch[cancer, "Library size - mRNA (thousands of reads)"] <- signif(median(lib_size_mRNA[[cancer]]) / 10^6, 1)

  # number of patients
  data_ch[cancer, "n patients"] <- nrow(clin)

  # censoring rate
  data_ch[cancer, "Censoring rate"] <- signif(sum(clin$status == 0)/nrow(clin), digits = 2)

  # survival at 3 years
  km_fit <- survfit(y_cox~1)
  km_fit_summary <- summary(km_fit, times = 3)
  data_ch[cancer, "Survival - 3 years"] <- signif(km_fit_summary$surv[1],2)

  # C-index
  data_ch[cancer, "C-index - EN"] <- signif(median(C_df_cancers_EN[, cancer], na.rm = T), digits = 2)
  data_ch[cancer, "C-index - RF"] <- signif(median(C_df_cancers_RF[, cancer], na.rm = T), digits = 2)

  # p miRNA / mRNA
  data_ch[cancer, "p miRNA"] <- ncol(count)
  data_ch[cancer, "p mRNA"] <- ncol(count_mRNA)

}

data_ch
data_ch[cancers_final,]

write.xlsx(data_ch[cancers_final,], file = "tables/data_ch.xlsx", rowNames = T)
print("Table saved in: tables/data_ch.xlsx")

n_patients <- data_ch[, "n patients"]
names(n_patients) <- row.names(data_ch)
save(n_patients, med_seq_depth, lib_size_mRNA, lib_size_miRNA, 
     file = "data_fit/seq_depth.RData")

