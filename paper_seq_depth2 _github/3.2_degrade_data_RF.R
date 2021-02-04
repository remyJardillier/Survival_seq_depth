
method <- "RF"

# learn models - RF --------------------------------------------------------


if(learn_new_models){
  
  # learn models for random forest ("RF").
  
  for(cancer in cancers_final[2:3]){
    
    print(paste0("*** Start learning for cancer: ", cancer, " ***"))
    
    # load the data
    source(file = "load_data/load_data_final.R")
    
    C_df_final = IBS_df_final = n_genes_df_final = n_patients_df_final <- 
      array(dim = c(n_rep * K_folds,
                    length(proportions), 
                    length(perc_train)),
            dimnames = list(1:(n_rep * K_folds), 
                            proportions,
                            perc_train))
    
    C_df_final_mRNA = IBS_df_final_mRNA = n_genes_df_final_mRNA = n_patients_df_final_mRNA <- 
      array(dim = c(n_rep * K_folds,
                    length(proportions), 
                    length(perc_train)),
            dimnames = list(1:(n_rep * K_folds), 
                            proportions,
                            perc_train))
    
    ind <- 1
    
    for(i in 1:n_rep){
      
      print(paste0("*** Start learning for repetition number: ", i, " (", cancer, ") ***"))
      
      flds <- createFolds(1:nrow(count), k = K_folds, list = TRUE, returnTrain = FALSE)
      
      for(j in 1:length(proportions)){
        
        print(paste0("Start learning for proportions: ", proportions[j], " (", cancer, ") ***"))
        
        for(p in 1:length(perc_train)){
          
          if(perc_train[p] == 0.8 | proportions[j] == 1){
            
            print(paste0("Start learning for percentage of patients: ", 
                         perc_train[p], " (", cancer, ") ***"))
            
            # subsample the patients - miRNA
            count_tmp <- t(generateSubsampledMatrix(counts = t(count),
                                                    proportion = proportions[j], 
                                                    seed = rnorm(1)))
            # remove genes with NA values - miRNA
            id_NA_gene <- apply(count_tmp, 2, function(x) sum(which(is.na(x))))
            id_NA_gene <- id_NA_gene > 0
            count_tmp <- count_tmp[, !id_NA_gene]
            
            # subsample the patients - mRNA
            count_tmp_mRNA <- t(generateSubsampledMatrix(counts = t(count_mRNA), 
                                                         proportion = proportions[j], 
                                                         seed = rnorm(1)))
            # remove genes with NA values - mRNA
            id_NA_gene <- apply(count_tmp_mRNA, 2, function(x) sum(which(is.na(x))))
            id_NA_gene <- id_NA_gene > 0
            sum(id_NA_gene)
            count_tmp_mRNA <- count_tmp_mRNA[, !id_NA_gene]
            
            for(k in 1:K_folds){
              
              print(paste0("Start learning for fold: ", k, " ***"))
              
              # build training and testing dataset
              id_test <- flds[[k]]
              
              # id of the training dataset
              id_train_all <- 1:nrow(count_tmp)
              id_train_all <- id_train_all[-id_test]
              id_train <- sample(id_train_all, 
                                 size = min(as.integer(perc_train[p]*nrow(count_tmp)),
                                            length(id_train_all)))
              
              # testing and training dataset
              count_train <- count_tmp[id_train,]
              count_test <- count_tmp[id_test,]
              clin_train <- clin[id_train,]
              clin_test <- clin[id_test,]
              dim(count_train)
              dim(count_test)
              
              # remove patients with a library size of 0
              id_train_0 <- apply(count_train, 1, function(x) sum(x) == 0)
              count_train <- count_train[!id_train_0, ]
              clin_train <- clin_train[!id_train_0, ]
              id_test_0 <- apply(count_test, 1, function(x) sum(x) == 0)
              count_test <- count_test[!id_test_0, ]
              clin_test <- clin_test[!id_test_0, ]
              
              # filtering ang logCPM in the training data
              logCPM_list <- log.cpm.cv(count_train, count_test)
              
              logCPM_train <- logCPM_list$logCPM_train
              sd_train <- apply(logCPM_train, 2, sd)
              logCPM_train <- logCPM_train[, sd_train != 0]
              
              logCPM_test <- logCPM_list$logCPM_test
              logCPM_test <- logCPM_test[, colnames(logCPM_train)]
              
              n_patients_df_final[ind + k - 1, j, p] <- nrow(logCPM_train)
              n_genes_df_final[ind + k - 1, j, p] <- ncol(logCPM_train)
              
              
              dim(clin_train)
              dim(logCPM_train)
              dim(clin_test)
              dim(logCPM_test)
              
              # learn a model 
              res <- try(fit <- learn_models(logCPM_train, clin_train, method = "RF"))
              
              # compute the p-value of the univariate Cox
              if(!inherits(res, "try-error")){
                
                # compute the C-index and the IBS 
                C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_train, 
                                        logCPM_test, method = "RF")
                
                C_df_final[ind + k - 1, j, p] <- C_IBS_tmp[1]
                IBS_df_final[ind + k - 1, j, p] <- C_IBS_tmp[2]
              }
              
              # mRNA
              print("Learning for mRNA")
              
              # testing and training dataset
              count_train <- count_tmp_mRNA[id_train,]
              count_test <- count_tmp_mRNA[id_test,]
              clin_train <- clin[id_train,]
              clin_test <- clin[id_test,]
              
              dim(count_train)
              dim(count_test)
              
              # remove patients with a library size of 0
              id_train_0 <- apply(count_train, 1, function(x) sum(x) == 0)
              count_train <- count_train[!id_train_0, ]
              clin_train <- clin_train[!id_train_0, ]
              id_test_0 <- apply(count_test, 1, function(x) sum(x) == 0)
              count_test <- count_test[!id_test_0, ]
              clin_test <- clin_test[!id_test_0, ]
              
              y_cox_train <- Surv(clin_train$time, clin_train$status)
              y_cox_test <- Surv(clin_test$time, clin_test$status)
              
              # filtering ang logCPM in the training data
              logCPM_list <- log.cpm.cv(count_train, count_test)
              
              logCPM_train <- logCPM_list$logCPM_train
              sd_train <- apply(logCPM_train, 2, sd)
              logCPM_train <- logCPM_train[, sd_train != 0]
              
              logCPM_test <- logCPM_list$logCPM_test
              logCPM_test <- logCPM_test[, colnames(logCPM_train)]
              
              n_patients_df_final_mRNA[ind + k - 1, j, p] <- nrow(logCPM_train)
              n_genes_df_final_mRNA[ind + k - 1, j, p] <- ncol(logCPM_train)
              
              # standardization
              logCPM_train_std <- std_train(logCPM_train)
              logCPM_test_std <- std_test(logCPM_train, logCPM_test)
              
              # 2500 genes with the lowest p-values
              p_val <- p_val_univCox_func(logCPM_train_std, y_cox_train)
              id_low_p_val_gene <- which.minn(p_val, 2500)
              
              logCPM_train <- logCPM_train[, id_low_p_val_gene]
              logCPM_test <- logCPM_test[, id_low_p_val_gene]
              
              # learn a model 
              res <- try(fit <- learn_models(logCPM_train, clin_train, method = "RF"))
              
              # compute the p-value of the univariate Cox
              if(!inherits(res, "try-error")){
                
                # compute the C-index and the IBS - choose if the data have to
                # be used standardized or not
                C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_train, 
                                        logCPM_test, method = "RF")
                
                C_df_final_mRNA[ind + k - 1, j, p] <- C_IBS_tmp[1]
                IBS_df_final_mRNA[ind + k - 1, j, p] <- C_IBS_tmp[2]
                
              }
              
            }
          }
        }
      }
      
      ind <- ind + K_folds
    }
    
    C_df_final
    IBS_df_final
    n_genes_df_final
    
    save(C_df_final, IBS_df_final, n_genes_df_final, n_patients_df_final, 
         file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_RF.RData"))
    print(paste0("Data saved in: data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_RF.RData"))
    
    save(C_df_final_mRNA, IBS_df_final_mRNA, n_genes_df_final_mRNA, n_patients_df_final_mRNA, 
         file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_RF.RData"))
    print(paste0("Data saved in: data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_RF.RData"))
  }
}


# plot for one cancer -----------------------------------------------------

cancer <- "ACC"

source(file = "load_data/load_data_final.R")

load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_RF.RData"))
load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_RF.RData"))

main_fig_C <- my_plot(C_df_final, C_df_final_mRNA, n_genes_df_final, count, count_mRNA, y_pval_mRNA = 0.95, 
                      y_pval_miRNA = 0.92, y_n_pat = 0.46, main = "C-index", legend_tick = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), 
                      ylim = c(0.45, 0.95), na.rm = T, title = cancer)
print(main_fig_C)

main_fig_IBS <- my_plot(IBS_df_final, IBS_df_final_mRNA, n_genes_df_final, count, count_mRNA, y_pval_mRNA = 0.35, 
                        y_pval_miRNA = 0.32, y_n_pat = 0, main = "IBS", legend_tick = seq(0, 0.4, by = 0.025),
                        ylim = c(0, 0.35), na.rm = T, title = cancer)
print(main_fig_IBS)
