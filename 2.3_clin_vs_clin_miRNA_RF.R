

# learn models ------------------------------------------------------------

if(learn_new_models){
  
  method <- "RF"
  
  for(cancer in cancers_final){
    
    print(paste0("Start learning for: ", cancer))
    
    # load the data ---
    source(file = "load_data/load_data_final.R")
    
    n_patients <- nrow(clin_data_cox)
    n_patients
    
    C_df_final = IBS_df_final <- 
      data.frame(matrix(ncol = 3, nrow = K_folds * n_rep))
    
    colnames(C_df_final) = colnames(IBS_df_final) <- 
      c("clin", "miRNA", "both")
    
    i=1
    k=1
    
    ind <- 1
    
    for(i in 1:n_rep){
      
      print(paste0("*** Start learning for repetition number: ", i, " ***"))
      
      flds <- createFolds(1:nrow(count), k = K_folds, list = TRUE, returnTrain = FALSE)
      
      for(k in 1:K_folds){
        
        print(paste0("Start learning for fold: ", k, " ***"))
        
        # build training and testing dataset
        id_test <- flds[[k]]
        
        # id of the training dataset
        id_train <- 1:nrow(count)
        id_train <- id_train[-id_test]
        
        # testing and training dataset
        count_train <- count[id_train,]
        count_test <- count[id_test,]
        
        clin_train <- clin[id_train,]
        clin_test <- clin[id_test,]
        
        clin_cox_train <- clin_data_cox[id_train,]
        clin_cox_test <- clin_data_cox[id_test,]
        
        dim(count_train)
        dim(count_test)
        
        # filtering ang logCPM in the training data
        logCPM_list <- log.cpm.cv(count_train, count_test)
        
        logCPM_train <- logCPM_list$logCPM_train
        sd_train <- apply(logCPM_train, 2, sd)
        logCPM_train <- logCPM_train[, sd_train != 0]
        
        logCPM_test <- logCPM_list$logCPM_test
        logCPM_test <- logCPM_test[, colnames(logCPM_train)]
        
        
        # miRNA only ---
        # learn a model 
        res_miRNA <- try(fit <- learn_models(logCPM_train, clin_train, "RF"))
      
        
        # compute the p-value of the univariate Cox
        if(!inherits(res_miRNA, "try-error")){
          
          data_ranger_train <- cbind(time = clin_train$time, status = clin_train$status, logCPM_train)
          data_ranger_test <- cbind(time = clin_test$time, status = clin_test$status, logCPM_test)
          
          colnames(data_ranger_train) <- str_replace_all(colnames(data_ranger_train), "-", "_")
          colnames(data_ranger_train) <- str_replace_all(colnames(data_ranger_train), "\\|", "_")
          colnames(data_ranger_train) <- str_replace_all(colnames(data_ranger_train), "\\?", "a")
          
          colnames(data_ranger_test) <- str_replace_all(colnames(data_ranger_test), "-", "_")
          colnames(data_ranger_test) <- str_replace_all(colnames(data_ranger_test), "\\|", "_")
          colnames(data_ranger_test) <- str_replace_all(colnames(data_ranger_test), "\\?", "a")
          
          # risck score (mean of the cumulative hazard function) for the C-index ---
          PI_test_miRNA <- predict(fit, data = data_ranger_test)
          PI_test_miRNA <- rowMeans(PI_test_miRNA$chf)
          
          PI_train_miRNA <- predict(fit, data = data_ranger_train)
          PI_train_miRNA <- rowMeans(PI_train_miRNA$chf)
          
          # compute the C-index and the IBS 
          C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_train, 
                                  logCPM_test, "RF")
          
          
          C_df_final[ind + k - 1, "miRNA"] <- C_IBS_tmp[1]
          IBS_df_final[ind + k - 1, "miRNA"] <- C_IBS_tmp[2]
          
        }
        
        # clin only ---
        # learn a model 
        res <- try(fit <- learn_models(clin_cox_train, clin_train, "RF"))
        
        # compute the p-value of the univariate Cox
        if(!inherits(res, "try-error")){
          
          # compute the C-index and the IBS 
          C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, clin_cox_train, 
                                  clin_cox_test, "RF")
          
          C_df_final[ind + k - 1, "clin"] <- C_IBS_tmp[1]
          IBS_df_final[ind + k - 1, "clin"] <- C_IBS_tmp[2]
        }
        
        # clin + miRNA ---
        # learn a model 
        if(!inherits(res_miRNA, "try-error")){
          
          miRNA_clin_train <- cbind(PI = PI_train_miRNA, clin_cox_train)
          miRNA_clin_test <- cbind(PI = PI_test_miRNA, clin_cox_test)
          
          res <- try(fit <- learn_models(miRNA_clin_train, clin_train, "RF"))
          
          # compute the p-value of the univariate Cox
          if(!inherits(res, "try-error")){
            
            # compute the C-index and the IBS 
            C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, miRNA_clin_train, 
                                    miRNA_clin_test, "RF")
            
            C_df_final[ind + k - 1, "both"] <- C_IBS_tmp[1]
            IBS_df_final[ind + k - 1, "both"] <- C_IBS_tmp[2]
          }
        }
        
        
      }
      
      ind <- ind + K_folds
    }
    
    boxplot(C_df_final, main = cancer)
    boxplot(IBS_df_final, main = cancer)
    
    save(C_df_final, IBS_df_final, n_patients,
         file = paste0("data_fit/", cancer, "/pred_clin_miRNA_RF.RData"))
    print(paste0("Data saved in: data_fit/", cancer, "/pred_clin_miRNA_RF.RData"))
  }
  
}



# p-values - Mix > clinical ? ---------------------------------------------

p_val_vect_C = p_val_vect_IBS <- rep(NA, length(cancers_final))
names(p_val_vect_C) = names(p_val_vect_IBS) <- cancers_final

for(cancer in cancers_final){
  
  load(file = paste0("data_fit/", cancer, "/pred_clin_miRNA_RF.RData"))
  
  p_val_vect_C[cancer] <- wilcox.test(C_df_final$both, C_df_final$clin, 
                                    alternative = "greater", paired = T)$p.value
  p_val_vect_IBS[cancer] <- wilcox.test(IBS_df_final$both, IBS_df_final$clin, 
                                      alternative = "less", paired = T)$p.value
}

# C-index
p_val_vect_C_BH <- p.adjust(p_val_vect_C, method = "BH")
stars_C <- sapply(p_val_vect_C_BH, my_stars.pval)

# IBS
p_val_vect_IBS_BH <- p.adjust(p_val_vect_IBS, method = "BH")
stars_IBS <- sapply(p_val_vect_IBS_BH, my_stars.pval)


# boxplot - C-index -----------------------------------------------------------------

# build the data frame
C_list = IBS_list <- list()

for(cancer in cancers_final){
  
  load(file = paste0("data_fit/", cancer, "/pred_clin_miRNA_RF.RData"))
  
  C_list[[cancer]] <- stack(C_df_final)
  IBS_list[[cancer]] <- stack(IBS_df_final)
}

C_final_ggplot <- melt(C_list)
IBS_final_ggplot <- melt(IBS_list)

# median C-index (order)
load(file = "data_fit/pred_ch_cancers.RData")
C_df_cancers <- C_df_cancers_EN[, cancers_final]
med_C <- apply(C_df_cancers, 2, function(x) median(x, na.rm = T))

# clinical data available for each cancer
load(file = "data_fit/clin_text_plot.RData")

col_cancer <- rep("black", length(cancers_final))
col_cancer[cancers_final %in% c("CESC", "PRAD", "UCEC")] <- "darkgray"

# ggplot - C-index
plot_clin_miRNA_C <- ggplot() + 
  geom_boxplot(data = C_final_ggplot, position = position_dodge(0.75), 
               aes(x = factor(L1, levels = names(sort(med_C, decreasing = T))), 
                   y = value, fill = ind)) + 
  theme_Publication_legend_bottom() + xlab("Cancer") + ylab("C-index")+
  theme(legend.title=element_blank(), legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = col_cancer)) + ylim(NA, 1) +
  scale_fill_manual(values=c("brown3","#386cb0", "orchid"),
                      labels = c("Clinical", "miRNA - RF", "Both")) +
  geom_text(aes(x = 1:length(cancers_final), y = rep(1, length(cancers_final)), 
                label = stars_C), col = "darkorchid", size = 5)  +
  geom_text(aes(x = 1:length(cancers_final), y = rep(0.15, length(cancers_final)), 
                label = clin_text_plot), col = "brown3", size = 4) +
  geom_bracket(xmin = length(cancers_final)-0.25, xmax = length(cancers_final)+0.25,
               y.position = 0.92, label = "", col = "orchid", size = 1)
  
print(plot_clin_miRNA_C)

ggsave(ggarrange(plot_clin_miRNA_C, labels = "A"), 
       filename = "pdf/clin_miARN_C_RF.pdf")
print("Figure saved in: pdf/clin_miARN_C_RF.pdf")

# ggplot - IBS
plot_clin_miRNA_IBS <- ggplot() + 
  geom_boxplot(data = IBS_final_ggplot, position = position_dodge(0.75), 
               aes(x = factor(L1, levels = names(sort(med_C, decreasing = T))), 
                   y = value, fill = ind)) + 
  theme_Publication() + xlab("Cancer") + ylab("IBS")+
  theme(legend.title=element_blank(), legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = col_cancer)) + ylim(NA, 0.4) +
  scale_fill_manual(values=c("brown3","#386cb0", "orchid"),
                    labels = c("Clinical", "miRNA - EN", "Mixed")) +
  geom_text(aes(x = 1:length(cancers_final), y = rep(0.4, length(cancers_final)), 
                label = stars_IBS), col = "darkorchid", size = 5)  +
  geom_text(aes(x = 1:length(cancers_final), y = rep(-0.05, length(cancers_final)), 
                label = clin_text_plot), col = "brown3", size = 4) +
  geom_bracket(xmin = length(cancers_final)-0.25, xmax = length(cancers_final)+0.25,
               y.position = 0.35, label = "", col = "orchid", size = 1)

print(plot_clin_miRNA_IBS)

ggsave(ggarrange(plot_clin_miRNA_IBS, labels = "B"), 
       filename = "pdf/clin_miARN_IBS_RF.pdf")
print("Figure saved in: pdf/clin_miARN_IBS_RF.pdf")