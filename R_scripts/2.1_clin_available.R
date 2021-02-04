# Clinical data available -------------------------------------------------

# data frame with all the clinical features
clin_available_miRNA <- data.frame(matrix(ncol = 6, nrow = length(cancers_final)))
colnames(clin_available_miRNA) <- c("age", "gender", "grade", "T_stage", "N_stage", "M_stage")
row.names(clin_available_miRNA) <- cancers_final

# vector for the plot
clin_text_plot <- rep(NA, length(cancers_final))
names(clin_text_plot) <- cancers_final

for(cancer in cancers_final){
  
  print(paste("*** start learning for", cancer, "***"))
  
  # load the data ---
  source(file = "load_data/load_data_final.R")
  
  clin_available_miRNA[cancer,] <- as.numeric(c("age", "gender", "grade", 
                                                "T_stage", "N_stage", "M_stage") %in%
                                                colnames(clin_data_cox))
  
  text_tmp <- ""
  if("grade" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "G")
  }
  if("T_stage" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "T")
  }
  if("N_stage" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "N")
  }
  if("M_stage" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "M")
  }
  
  clin_text_plot[cancer] <- text_tmp
  
}

clin_available_miRNA

write.xlsx(clin_available_miRNA, file = "tables/clin_available_miRNA.xlsx",
           row.names = T)
print("Table saved in: tables/clin_available_miRNA.xlsx")

save(clin_text_plot, file = "data_fit/clin_text_plot.RData")
