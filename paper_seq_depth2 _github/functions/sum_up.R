

# p-values : differences between without and with degradation -------------

p_val_deg_func <- function(df_final, metric, deg_type){
  
  # data frame
  if(deg_type == "seq"){
    prop_df <- df_final[, , "0.8"]
  }else if(deg_type == "patients"){
    prop_df <- df_final[, "1", ]
  }
  
  # compute the one-sided p-value
  if(metric == "IBS"){
    alternative <- "less"
  }else{
    alternative <- "greater"
  }
  
  # compute the p-value
  p_val_prop_test <- rep(NA, ncol(prop_df) - 1)
  for(i in 1:(ncol(prop_df)-1)){
    
    if(sum(is.na(prop_df[,i])) < 45){
      p_val_prop_test[i] <- wilcox.test(prop_df[,ncol(prop_df)], prop_df[,i], 
                                        alternative = alternative)$p.value
    }else{
      p_val_prop_test[i] <- NA
    }
    
  }
  
  return(p_val_prop_test)
}

# miRNA - sequencing depth -----------------------------------------------------

sum_up_miRNA_func_seq <- function(method, metric, signif_level){
  
  # table with the sum up
  table_sum_up <- data.frame(matrix(ncol = length(cancers_final), nrow = 2))
  colnames(table_sum_up) <- cancers_final
  row.names(table_sum_up) <- c("d_factor", "med_lib_size")
  
  for(cancer in cancers_final){
    
    print(paste0("Start for cancer:", cancer))
    
    # compute the p-values
    if(method == "EN"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_EN.RData"))
    }else if(method == "RF"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_RF.RData"))
    }
    
    if(metric == "C"){
      p_val <- p_val_deg_func(C_df_final, "C", "seq")
    }else if(metric == "IBS"){
      p_val <- p_val_deg_func(IBS_df_final, "IBS", "seq")
    }else if(metric == "both"){
      p_val_C <- p_val_deg_func(C_df_final, "C", "seq")
      p_val_IBS <- p_val_deg_func(IBS_df_final, "IBS", "seq")
    }
    
    # degradation factor from which the metric is degradated
    if(metric == "C" | metric == "IBS"){
      id <- p_val < signif_level  
    }else if(metric == "both"){
      id_C <- p_val_C < signif_level
      id_IBS <- p_val_IBS < signif_level
      id <- id_C | id_IBS
    }
    
    d_factor <- as.numeric(dimnames(C_df_final)[[2]][max(which(id)) + 1]) 
    
    # # if it is impossible to degrade by a factor 10, load the data for a degradation of 5 and 2
    # if(!is.na(d_factor) & d_factor == 1){
    #   
    #   if(method == "EN"){
    #     load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_5_2_EN.RData"))
    #   }else if(method == "RF"){
    #     load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_5_2_RF.RData"))
    #   }
    #   
    #   
    #   if(metric == "C"){
    #     ref_C <- C_df_final[, "1", "0.8"]
    #     p_val_C <- apply(C_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_C, x, alternative = "greater")$p.value)
    #     id <- p_val_C < signif_level
    #   }else if(metric == "IBS"){
    #     ref_IBS <- IBS_df_final[, "1", "0.8"]
    #     p_val_IBS <- apply(IBS_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_IBS, x, alternative = "less")$p.value)
    #     id <- p_val_IBS < signif_level
    #   }else if(metric == "both"){
    #     ref_C <- C_df_final[, "1", "0.8"]
    #     p_val_C <- apply(C_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_C, x, alternative = "greater")$p.value)
    #     id_C <- p_val_C < signif_level
    #     
    #     ref_IBS <- IBS_df_final[, "1", "0.8"]
    #     p_val_IBS <- apply(IBS_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_IBS, x, alternative = "less")$p.value)
    #     id_IBS <- p_val_IBS < signif_level
    #     
    #     id_C <- p_val_C < signif_level
    #     id_IBS <- p_val_IBS < signif_level
    #     id <- id_C | id_IBS
    #   }
    #   
    #   
    #   if(length(which(id)) == 0){
    #     d_factor <- as.numeric(colnames(C_df_final_5_2)[1]) 
    #   }else if(max(which(id)) == 1){
    #     d_factor <- as.numeric(colnames(C_df_final_5_2)[2]) 
    #   }else if(max(which(id)) == 2){
    #     d_factor <- 1
    #   }
    # }
    
    table_sum_up["d_factor", cancer] <- 1/d_factor
    
    # median library size ---
    table_sum_up["med_lib_size", cancer] <- 
      round(median(med_seq_depth[cancer, "miRNA"] * d_factor))
  }
  
  table_sum_up_round <- table_sum_up
  table_sum_up_round["med_lib_size", ] <- signif(ceiling(as.numeric(table_sum_up_round["med_lib_size", ])/1000), 1)
  table_sum_up_round
  
  return(table_sum_up_round)
}

# miRNA - % of patients -----------------------------------------------------

sum_up_miRNA_func_patients <- function(method, metric, signif_level){
  
  # table with the sum up
  table_sum_up <- data.frame(matrix(ncol = length(cancers_final), nrow = 2))
  colnames(table_sum_up) <- cancers_final
  row.names(table_sum_up) <- c("perc", "n_patients")
  
  for(cancer in cancers_final){
    
    print(paste0("Start for cancer:", cancer))
    
    # compute the p-values
    if(method == "EN"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_EN.RData"))
    }else if(method == "RF"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_RF.RData"))
    }
    
    if(metric == "C"){
      p_val <- p_val_deg_func(C_df_final, "C", "patients")
    }else if(metric == "IBS"){
      p_val <- p_val_deg_func(IBS_df_final, "IBS", "patients")
    }else if(metric == "both"){
      p_val_C <- p_val_deg_func(C_df_final, "C", "patients")
      p_val_IBS <- p_val_deg_func(IBS_df_final, "IBS", "patients")
    }
    
    # degradation factor from which the metric is degradated
    if(metric == "C" | metric == "IBS"){
      id <- p_val < signif_level  
    }else if(metric == "both"){
      id_C <- p_val_C < signif_level
      id_IBS <- p_val_IBS < signif_level
      id <- id_C | id_IBS
    }
    
    if(sum(id, na.rm = T) == 0){
      perc <- NA
    }else{
      perc <- as.numeric(dimnames(C_df_final)[[3]][max(which(id)) + 1]) 
    }
    
    table_sum_up["perc", cancer] <- perc
    
    # number of patients ---
    table_sum_up["n_patients", cancer] <- round(perc * n_patients[cancer])
  }
  
  return(table_sum_up)
}

# mRNA - sequencing depth -------------------------------------------------------------------

sum_up_mRNA_func_seq <- function(method, metric, signif_level){
  
  # table with the sum up
  table_sum_up <- data.frame(matrix(ncol = length(cancers_final), nrow = 2))
  colnames(table_sum_up) <- cancers_final
  row.names(table_sum_up) <- c("d_factor", "med_lib_size")
  
  for(cancer in cancers_final){
    
    print(paste0("Start for cancer:", cancer))
    
    # compute the p-values
    if(method == "EN"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_EN.RData"))
    }else if(method == "RF"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_RF.RData"))
    }
    
    if(metric == "C"){
      p_val <- p_val_deg_func(C_df_final, "C", "seq")
    }else if(metric == "IBS"){
      p_val <- p_val_deg_func(IBS_df_final, "IBS", "seq")
    }else if(metric == "both"){
      p_val_C <- p_val_deg_func(C_df_final, "C", "seq")
      p_val_IBS <- p_val_deg_func(IBS_df_final, "IBS", "seq")
    }
    
    # degradation factor from which the metric is degradated
    if(metric == "C" | metric == "IBS"){
      id <- p_val < signif_level  
    }else if(metric == "both"){
      id_C <- p_val_C < signif_level
      id_IBS <- p_val_IBS < signif_level
      id <- id_C | id_IBS
    }
    
    if(sum(id) == 0){
      d_factor <- NA
    }else{
      d_factor <- as.numeric(dimnames(C_df_final_mRNA)[[2]][max(which(id)) + 1])  
    }
    
    # # if it is impossible to degrade by a factor 10, load the data for a degradation of 5 and 2
    # if(!is.na(d_factor) & d_factor == 1){
    #   
    #   if(method == "EN"){
    #     load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_5_2_EN.RData"))
    #   }else if(method == "RF"){
    #     load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_5_2_RF.RData"))
    #   }
    #   
    #   if(metric == "C"){
    #     ref_C <- C_df_final_mRNA[, "1", "0.8"]
    #     p_val_C <- apply(C_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_C, x, alternative = "greater")$p.value)
    #     id <- p_val_C < signif_level
    #   }else if(metric == "IBS"){
    #     ref_IBS <- IBS_df_final_mRNA[, "1", "0.8"]
    #     p_val_IBS <- apply(IBS_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_IBS, x, alternative = "less")$p.value)
    #     id <- p_val_IBS < signif_level
    #   }else if(metric == "both"){
    #     ref_C <- C_df_final[, "1", "0.8"]
    #     p_val_C <- apply(C_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_C, x, alternative = "greater")$p.value)
    #     id_C <- p_val_C < signif_level
    #     
    #     ref_IBS <- IBS_df_final[, "1", "0.8"]
    #     p_val_IBS <- apply(IBS_df_final_5_2, 2, function(x)
    #       wilcox.test(ref_IBS, x, alternative = "less")$p.value)
    #     id_IBS <- p_val_IBS < signif_level
    #     
    #     id_C <- p_val_C < signif_level
    #     id_IBS <- p_val_IBS < signif_level
    #     id <- id_C | id_IBS
    #   }
    #   
    #   
    #   if(length(which(id)) == 0){
    #     d_factor <- as.numeric(colnames(C_df_final_5_2)[1]) 
    #   }else if(max(which(id)) == 1){
    #     d_factor <- as.numeric(colnames(C_df_final_5_2)[2]) 
    #   }else if(max(which(id)) == 2){
    #     d_factor <- 1
    #   }
    # }
    
    table_sum_up["d_factor", cancer] <- 1/d_factor
    
    # median library size ---
    table_sum_up["med_lib_size", cancer] <- 
      round(median(med_seq_depth[cancer, "mRNA"] * d_factor))
  }
  
  table_sum_up_round <- table_sum_up
  table_sum_up_round["med_lib_size", ] <- signif(ceiling(as.numeric(table_sum_up_round["med_lib_size", ])/1000), 1)
  table_sum_up_round
  
  return(table_sum_up_round)
}


# mRNA - % of patients -----------------------------------------------------

sum_up_mRNA_func_patients <- function(method, metric, signif_level){
  
  # table with the sum up
  table_sum_up <- data.frame(matrix(ncol = length(cancers_final), nrow = 2))
  colnames(table_sum_up) <- cancers_final
  row.names(table_sum_up) <- c("perc", "n_patients")
  
  for(cancer in cancers_final){
    
    print(paste0("Start for cancer:", cancer))
    
    # compute the p-values
    if(method == "EN"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_EN.RData"))
    }else if(method == "RF"){
      load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_RF.RData"))
    }
    
    if(metric == "C"){
      p_val <- p_val_deg_func(C_df_final_mRNA, "C", "patients")
    }else if(metric == "IBS"){
      p_val <- p_val_deg_func(IBS_df_final_mRNA, "IBS", "patients")
    }else if(metric == "both"){
      p_val_C <- p_val_deg_func(C_df_final_mRNA, "C", "patients")
      p_val_IBS <- p_val_deg_func(IBS_df_final_mRNA, "IBS", "patients")
    }
    
    # degradation factor from which the metric is degradated
    if(metric == "C" | metric == "IBS"){
      id <- p_val < signif_level  
    }else if(metric == "both"){
      id_C <- p_val_C < signif_level
      id_IBS <- p_val_IBS < signif_level
      id <- id_C | id_IBS
    }
    
    
    if(sum(id, na.rm = T) == 0){
      perc <- NA
    }else{
      perc <- as.numeric(dimnames(C_df_final_mRNA)[[3]][max(which(id)) + 1]) 
    }
    
    table_sum_up["perc", cancer] <- perc
    
    # number of patients ---
    table_sum_up["n_patients", cancer] <- round(perc * n_patients[cancer])
  }
  
  return(table_sum_up)
}
