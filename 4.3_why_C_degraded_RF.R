


# test hypothesis --------------------------------------------------------

if(learn_new_models){
  
  method <- "RF"
  
  C_list_testH = IBS_list_testH = n_genes_list_testH <- list()
  
  for(cancer in cancers_final){
    
    print(paste0("*** Start learning for : ", cancer))
    
    source(file = "load_data/load_data_final.R")
    
    C_df_testH = IBS_df_testH =  n_genes_df_testH <- 
      data.frame(matrix(ncol = 3, nrow = n_rep * K_folds))
    colnames(C_df_testH) = colnames(IBS_df_testH) = colnames(n_genes_df_testH) <- 
      c("sub", "MEG", "original") # MEG = most expressed miRNA
    
    ind <- 1
    
    for(i in 1:n_rep){
      
      print(paste0("*** Start learning for repetition number: ", i, " ***"))
      
      flds <- createFolds(1:nrow(count), k = K_folds, list = TRUE, returnTrain = FALSE)
      
      # subsample the patients - miRNA
      count_tmp <- t(generateSubsampledMatrix(counts = t(count),
                                              proportion = 10^(-4), 
                                              seed = rnorm(1)))
      # remove genes with NA values - miRNA
      id_NA_gene <- apply(count_tmp, 2, function(x) sum(which(is.na(x))))
      id_NA_gene <- id_NA_gene > 0
      count_tmp <- count_tmp[, !id_NA_gene]
      
      for(k in 1:K_folds){
        
        print(paste0("Start learning for fold: ", k, " ***"))
        
        # build training and testing dataset
        id_test <- flds[[k]]
        
        # testing and training dataset
        count_sub_train <- count_tmp[-id_test,]
        count_sub_test <- count_tmp[id_test,]
        count_train <- count[-id_test,]
        count_test <- count[id_test,]
        clin_train <- clin[-id_test,]
        clin_test <- clin[id_test,]
        dim(count_train)
        dim(count_test)
        dim(count_sub_train)
        dim(count_sub_test)
        
        # filtering ang logCPM in the training data - subsample
        logCPM_sub_list <- log.cpm.cv(count_sub_train, count_sub_test)
        
        logCPM_sub_train <- logCPM_sub_list$logCPM_train
        sd_train <- apply(logCPM_sub_train, 2, sd)
        logCPM_sub_train <- logCPM_sub_train[, sd_train != 0]
        
        logCPM_sub_test <- logCPM_sub_list$logCPM_test
        logCPM_sub_test <- logCPM_sub_test[, colnames(logCPM_sub_train)]
        
        dim(logCPM_sub_test)
        dim(logCPM_sub_train)
        
        # filtering ang logCPM in the training data
        logCPM_list <- log.cpm.cv(count_train, count_test)
        
        logCPM_train <- logCPM_list$logCPM_train
        sd_train <- apply(logCPM_train, 2, sd)
        logCPM_train <- logCPM_train[, sd_train != 0]
        
        logCPM_test <- logCPM_list$logCPM_test
        logCPM_test <- logCPM_test[, colnames(logCPM_train)]
        
        dim(logCPM_test)
        dim(logCPM_train)
        
        # filtering ang logCPM in the training data - most expressed genes (MEG)
        logCPM_MEG_list <- log.cpm.cv(count_train[, colnames(logCPM_sub_train)], 
                                      count_test[, colnames(logCPM_sub_train)])
        
        logCPM_MEG_train <- logCPM_MEG_list$logCPM_train
        sd_train <- apply(logCPM_MEG_train, 2, sd)
        logCPM_MEG_train <- logCPM_MEG_train[, sd_train != 0]
        
        logCPM_MEG_test <- logCPM_MEG_list$logCPM_test
        logCPM_MEG_test <- logCPM_MEG_test[, colnames(logCPM_MEG_train)]
        
        dim(logCPM_MEG_test)
        dim(logCPM_MEG_train)
        
        # standardization
        logCPM_train_std <- std_train(logCPM_train)
        logCPM_test_std <- std_test(logCPM_train, logCPM_test)
        
        logCPM_sub_train_std <- std_train(logCPM_sub_train)
        logCPM_sub_test_std <- std_test(logCPM_sub_train, logCPM_sub_test)
        
        logCPM_MEG_train_std <- std_train(logCPM_MEG_train)
        logCPM_MEG_test_std <- std_test(logCPM_MEG_train, logCPM_MEG_test)
        
        dim(logCPM_train_std)
        dim(logCPM_sub_train_std)
        dim(logCPM_MEG_train_std)
        
        n_genes_df_testH[ind + k - 1, ] <- c(ncol(logCPM_sub_train_std),
                                             ncol(logCPM_MEG_train_std),
                                             ncol(logCPM_train_std))
        
        # compute the predictive measure - sub ---
        res <- try(fit <- learn_models(logCPM_sub_train, clin_train, "RF"))
        
        
        if(!inherits(res, "try-error")){
          
          # compute the C-index and the IBS - choose if the data have to
          # be used standardized or not
          C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_sub_train, 
                                    logCPM_sub_test, "RF")
          
          C_df_testH[ind + k - 1, "sub"] <- C_IBS_tmp[1]
          IBS_df_testH[ind + k - 1, "sub"] <- C_IBS_tmp[2]
          
        }
        
        # compute the predictive measure - MEG ---
        res <- try(fit <- learn_models(logCPM_MEG_train, clin_train, "RF"))
        
        if(!inherits(res, "try-error")){
          
          # compute the C-index and the IBS - choose if the data have to
          # be used standardized or not
          C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_MEG_train, 
                                    logCPM_MEG_test, "RF")
          
          C_df_testH[ind + k - 1, "MEG"] <- C_IBS_tmp[1]
          IBS_df_testH[ind + k - 1, "MEG"] <- C_IBS_tmp[2]
          
        }
        
        # compute the predictive measure - original ---
        res <- try(fit <- learn_models(logCPM_train, clin_train, "RF"))
        
        if(!inherits(res, "try-error")){
          
          # compute the C-index and the IBS - choose if the data have to
          # be used standardized or not
          C_IBS_tmp <- C_IBS_func(fit, clin_train, clin_test, logCPM_train, 
                                    logCPM_test, "RF")
          
          C_df_testH[ind + k - 1, "original"] <- C_IBS_tmp[1]
          IBS_df_testH[ind + k - 1, "original"] <- C_IBS_tmp[2]
        }
        
      }
      
      ind <- ind + K_folds
      
    }
    
    C_list_testH[[cancer]] <- melt(C_df_testH)
    IBS_list_testH[[cancer]] <- melt(IBS_df_testH)
    n_genes_list_testH[[cancer]] <- melt(n_genes_df_testH)
    # boxplot(C_df_testH)
  }
  
  
  save(C_list_testH, IBS_list_testH, n_genes_list_testH,
       file = "data_fit/testH_deg_RF.RData")
  print("Data saved in: data_fit/testH_deg_RF.RData")
}



# Test if the differences are significative -------------------------------

load(file = "data_fit/testH_deg_RF.RData")

test_H1_C = test_H2_C = test_H1_IBS = test_H2_IBS <- rep(NA, length(cancers_final))
names(test_H1_C) = names(test_H2_C) = names(test_H1_IBS) = names(test_H2_IBS) <- cancers_final

for(cancer in cancers_final){
  
  # C-index
  MEG <- C_list_testH[[cancer]]$value[C_list_testH[[cancer]]$variable == "MEG"]
  sub <- C_list_testH[[cancer]]$value[C_list_testH[[cancer]]$variable == "sub"] 
  original <- C_list_testH[[cancer]]$value[C_list_testH[[cancer]]$variable == "original"] 
  
  test_H1_C[cancer] <- wilcox.test(MEG, sub, alternative = "greater", paired = T)$p.value
  test_H2_C[cancer] <- wilcox.test(original, MEG,alternative = "greater", paired = T)$p.value
  
  # IBS
  MEG <- IBS_list_testH[[cancer]]$value[IBS_list_testH[[cancer]]$variable == "MEG"]
  sub <- IBS_list_testH[[cancer]]$value[IBS_list_testH[[cancer]]$variable == "sub"] 
  original <- IBS_list_testH[[cancer]]$value[IBS_list_testH[[cancer]]$variable == "original"] 
  
  test_H1_IBS[cancer] <- wilcox.test(MEG, sub, alternative = "less", paired = T)$p.value
  test_H2_IBS[cancer] <- wilcox.test(original, MEG,alternative = "less", paired = T)$p.value
}

# Benjamini-Hochberg correction
test_H1_C_BH <- p.adjust(test_H1_C, method = "BH")
test_H2_C_BH <- p.adjust(test_H2_C, method = "BH")

test_H1_IBS_BH <- p.adjust(test_H1_IBS, method = "BH")
test_H2_IBS_BH <- p.adjust(test_H2_IBS, method = "BH")

# stars for ggplot
stars_H1_C <- sapply(test_H1_C_BH, my_stars.pval)
stars_H2_C <- sapply(test_H2_C_BH, my_stars.pval)

stars_H1_IBS <- sapply(test_H1_IBS_BH, my_stars.pval)
stars_H2_IBS <- sapply(test_H2_IBS_BH, my_stars.pval)


# plot --------------------------------------------------------------------

n_med_genes <- rep(NA, length(cancers_final))
names(n_med_genes) <- cancers_final

n_genes_sub_df <- data.frame(matrix(ncol = length(cancers_final), nrow = 50))
colnames(n_genes_sub_df) <- cancers_final

n_genes_all_df <- data.frame(matrix(ncol = length(cancers_final), nrow = 50))
colnames(n_genes_all_df) <- cancers_final

# median number of genes detected
for(cancer in cancers_final){
  n_genes_df_tmp <- n_genes_list_testH[[cancer]]
  n_genes_all_df[, cancer] <- n_genes_df_tmp$value[n_genes_df_tmp$variable == "original"]
  
  n_genes_df_tmp <- n_genes_df_tmp$value[n_genes_df_tmp$variable == "sub"]
  n_genes_sub_df[, cancer] <- n_genes_df_tmp
  
  n_med_genes[cancer] <- round(median(n_genes_df_tmp, na.rm = T))
}

mean(as.matrix(n_genes_sub_df))
mean(as.matrix(n_genes_all_df))

# median C-index (order)
load(file = "data_fit/pred_ch_cancers.RData")
C_df_cancers <- C_df_cancers_EN[, cancers_final]
med_C <- apply(C_df_cancers, 2, function(x) median(x, na.rm = T))

# C-index
C_list_testH_ggplot <- melt(C_list_testH)

why_C_deg_plot <-  ggplot() +
  geom_boxplot(data = C_list_testH_ggplot, aes(x = factor(L1, levels = names(sort(med_C, decreasing = T))), 
                                               y = value, fill = variable),
               position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("C-index") + ylim(0.3, 1) +
  theme_Publication() + theme(axis.text.x = element_text(size = 12)) +
  geom_text(aes(x = 1:length(C_list_testH), y = rep(1, length(C_list_testH)), 
                label = stars_H1_C), col = "#7fc97f", size = 4) +
  geom_text(aes(x = 1:length(C_list_testH), y = rep(0.97, length(C_list_testH)), 
                label = stars_H2_C), col = "#386cb0", size = 5) +
  geom_text(aes(x = 1:length(C_list_testH), y = rep(0.3, length(C_list_testH)), 
                label = n_med_genes), size = 4) +
  scale_fill_manual(values=c("#7fc97f","#fdb462","#386cb0")) +
  geom_bracket(xmin = 11.05, xmax = 11.25, y.position = 0.92, label = "", col = "#386cb0", size = 1) +
  geom_bracket(xmin = 10.75, xmax = 10.95, y.position = 0.92, label = "", col = "#7fc97f", size = 1)
print(why_C_deg_plot)

ggsave(ggarrange(why_C_deg_plot, labels = "A"), 
       filename = "pdf/why_C_deg_plot_miRNA_RF.pdf")
print("Figure saved in: pdf/why_C_deg_plot_miRNA_RF.pdf"

# IBS
IBS_list_testH_ggplot <- melt(IBS_list_testH)

why_IBS_deg_plot <- ggplot() +
  geom_boxplot(data = IBS_list_testH_ggplot, aes(x = factor(L1, levels = names(sort(med_C, decreasing = T))), 
                                                 y = value, fill = variable),
               position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("IBS") + ylim(0,0.35) + 
  theme_Publication() + theme(axis.text.x = element_text(size = 13)) +
  geom_text(aes(x = 1:length(IBS_list_testH), y = rep(0.35, length(IBS_list_testH)), 
                label = stars_H1_IBS), col = "#7fc97f", size = 4) +
  geom_text(aes(x = 1:length(C_list_testH), y = rep(0.32, length(C_list_testH)), 
                label = stars_H2_IBS), col = "#386cb0", size = 5) +
  geom_text(aes(x = 1:length(C_list_testH), y = rep(0, length(C_list_testH)), 
                label = n_med_genes), size = 4) +
  scale_fill_manual(values=c("#7fc97f","#fdb462","#386cb0")) +
  geom_bracket(xmin = 11.05, xmax = 11.25, y.position = 0.3, label = "", col = "#386cb0", size = 1) +
  geom_bracket(xmin = 10.75, xmax = 10.95, y.position = 0.3, label = "", col = "#7fc97f", size = 1)

print(why_IBS_deg_plot)

ggsave(ggarrange(why_IBS_deg_plot, labels = "B"), 
       filename = "pdf/why_IBS_deg_plot_miRNA_RF.pdf")
print("Figure saved in: pdf/why_IBS_deg_plot_miRNA_RF.pdf")
