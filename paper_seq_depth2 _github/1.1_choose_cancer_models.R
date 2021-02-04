

# learn the models and compute prediction measures ------------------------

if(learn_new_models){
  
  C_df_cancers_EN = IBS_df_cancers_EN = 
    C_df_cancers_RF = IBS_df_cancers_RF <- 
    data.frame(matrix(ncol = length(cancers_all), nrow = n_rep * K_folds))
  colnames(C_df_cancers_EN) = colnames(IBS_df_cancers_EN) =
    colnames(C_df_cancers_RF) = colnames(IBS_df_cancers_RF) <- 
    cancers_all
  
  
  for(cancer in cancers_all){
    
    print(paste0("*** Start learning for cancer : ", cancer, " ***"))
    
    # load the data ---
    source(file = "load_data/load_data_final.R")
    
    ind <- 1
    
    for(i in 1:n_rep){
      
      print(paste0("Start learning for repetition number: ", i, " ***"))
      
      flds <- createFolds(1:nrow(count), k = K_folds, list = TRUE, returnTrain = FALSE)
      
      for(k in 1:K_folds){
        
        print(paste0("Start learning for fold: ", k, " ***"))
        
        # build training and testing dataset
        id_test <- flds[[k]]
        
        # testing and training dataset
        count_train <- count[-id_test,]
        count_test <- count[id_test,]
        
        clin_train <- clin[-id_test,]
        clin_test <- clin[id_test,]
        
        dim(count_train)
        dim(count_test)
        
        # filtering ang logCPM in the training data
        logCPM_list <- log.cpm.cv(count_train, count_test)
        
        logCPM_train <- logCPM_list$logCPM_train
        sd_train <- apply(logCPM_train, 2, sd)
        logCPM_train <- logCPM_train[, sd_train != 0]
        
        logCPM_test <- logCPM_list$logCPM_test
        logCPM_test <- logCPM_test[, colnames(logCPM_train)]
        
        # standardization
        logCPM_train_std <- std_train(logCPM_train)
        logCPM_test_std <- std_test(logCPM_train, logCPM_test)
        
        # prediction measures - EN ---
        res <- try(fit <- learn_models(logCPM_train_std, clin_train, "EN"))
        
        # compute the prediction accuracy
        if(!inherits(res, "try-error")){
          
          # compute the C-index and the IBS - choose if the data have to
          # be used standardized or not
          C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_train_std, 
                                  logCPM_test_std, "EN")
          
          C_df_cancers_EN[ind, cancer] <- C_IBS_tmp["C"]
          IBS_df_cancers_EN[ind, cancer] <- C_IBS_tmp["IBS"]
          
        }
        
        # prediction measures - RF ---
        res <- try(fit <- learn_models(logCPM_train, clin_train, "RF"))
        
        # compute the prediction accuracy
        if(!inherits(res, "try-error")){
          
          # compute the C-index and the IBS - choose if the data have to
          # be used standardized or not
          C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_train, 
                                  logCPM_test, "RF")
          
          C_df_cancers_RF[ind, cancer] <- C_IBS_tmp["C"]
          IBS_df_cancers_RF[ind, cancer] <- C_IBS_tmp["IBS"]
          
        }
        
        ind <- ind + 1
      }
    }
  }
  
  save(C_df_cancers_EN, IBS_df_cancers_EN, 
       C_df_cancers_RF, IBS_df_cancers_RF,
       file = "data_fit/pred_ch_cancers.RData")
  print("Data saved in: data_fit/pred_ch_cancers.RData")
}



# Choose cancer -----------------------------------------------------------

# Compute p-values of wilcoxon tests for both Cox and random forest (RF)

# EN ---
load(file = "data_fit/pred_ch_cancers.RData")

# compute the p-values
p_val_EN <- rep(NA, ncol(C_df_cancers_EN))
names(p_val_EN) <- colnames(C_df_cancers_EN)
for(i in 1:ncol(C_df_cancers_EN)){
  p_val_EN[i] <- wilcox.test(C_df_cancers_EN[,i], mu = 0.6, 
                             alternative = "greater")$p.value
}
p_val_EN
p_val_EN_BH <- p.adjust(p_val_EN, method = "BH")
sum(p_val_EN_BH <= signif_level)
length(p_val_EN_BH)
p_val_EN_BH[p_val_EN_BH <= signif_level]


# RF ---
load(file = "data_fit/pred_ch_cancers.RData")

# compute the p-values
p_val_RF <- rep(NA, ncol(C_df_cancers_RF))
names(p_val_RF) <- colnames(C_df_cancers_RF)
for(i in 1:ncol(C_df_cancers_RF)){
  p_val_RF[i] <- wilcox.test(C_df_cancers_RF[,i], mu = 0.6, 
                             alternative = "greater")$p.value
}
p_val_RF
p_val_RF_BH <- p.adjust(p_val_RF, method = "BH")
sum(p_val_RF_BH <= signif_level)
length(p_val_RF_BH)
p_val_RF_BH[p_val_RF_BH <= signif_level]

# sort by decreasing order of the median C-index - EN
med_C_EN <- apply(C_df_cancers_EN, 2, function(x) median(x, na.rm = T))
med_C_RF <- apply(C_df_cancers_RF, 2, function(x) median(x, na.rm = T))

cancers_all <- names(sort(med_C_EN, decreasing = T))
p_val_EN_BH <- p_val_EN_BH[cancers_all]
p_val_RF_BH <- p_val_RF_BH[cancers_all]

# cancer choosed 
cancers_final <- cancers_all[p_val_EN_BH <= signif_level | p_val_RF_BH <= signif_level]
print(cancers_final)
length(cancers_final)

