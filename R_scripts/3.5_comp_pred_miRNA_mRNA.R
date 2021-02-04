
# data - miRNA ------------------------------------------------------------

# load data from "3.1_degrade_data_EN.R" and "3.2_degrade_data_RF.R"
C_df_cancers_EN_miRNA = IBS_df_cancers_EN_miRNA = 
  C_df_cancers_RF_miRNA = IBS_df_cancers_RF_miRNA <- 
  data.frame(matrix(ncol = length(cancers_final), nrow = n_rep * K_folds))
colnames(C_df_cancers_EN_miRNA) = colnames(IBS_df_cancers_EN_miRNA) =
  colnames(C_df_cancers_RF_miRNA) = colnames(IBS_df_cancers_RF_miRNA) <- 
  cancers_final

for(cancer in cancers_final){
  
  # Cox
  load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_EN.RData"))
  C_df_cancers_EN_miRNA[, cancer] <- C_df_final[, length(proportions), length(perc_train)]
  IBS_df_cancers_EN_miRNA[, cancer] <- IBS_df_final[, length(proportions), length(perc_train)]
  
  # random survival forest
  load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_miRNA_RF.RData"))
  C_df_cancers_RF_miRNA[, cancer] <- C_df_final[, length(proportions), length(perc_train)]
  IBS_df_cancers_RF_miRNA[, cancer] <- IBS_df_final[, length(proportions), length(perc_train)]
}

# data - mRNA -------------------------------------------------------------

# load data from "3.1_degrade_data_EN.R" and "3.2_degrade_data_RF.R"
C_df_cancers_EN_mRNA = IBS_df_cancers_EN_mRNA = 
  C_df_cancers_RF_mRNA = IBS_df_cancers_RF_mRNA <- 
  data.frame(matrix(ncol = length(cancers_final), nrow = n_rep * K_folds))
colnames(C_df_cancers_EN_mRNA) = colnames(IBS_df_cancers_EN_mRNA) =
  colnames(C_df_cancers_RF_mRNA) = colnames(IBS_df_cancers_RF_mRNA) <- 
  cancers_final

for(cancer in cancers_final){
  
  # Cox
  load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_EN.RData"))
  C_df_cancers_EN_mRNA[, cancer] <- C_df_final_mRNA[, length(proportions), length(perc_train)]
  IBS_df_cancers_EN_mRNA[, cancer] <- IBS_df_final_mRNA[, length(proportions), length(perc_train)]
  
  # random survival forest
  load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_RF.RData"))
  C_df_cancers_RF_mRNA[, cancer] <- C_df_final_mRNA[, length(proportions), length(perc_train)]
  IBS_df_cancers_RF_mRNA[, cancer] <- IBS_df_final_mRNA[, length(proportions), length(perc_train)]
}


# plot - Cox --------------------------------------------------------------

# p-values ---
p_val_miRNA_mRNA_C = p_val_miRNA_mRNA_IBS <- rep(NA, length(cancers_final))
col_p_val_C = col_p_val_IBS <- rep(NA, length(cancers_final))
names(p_val_miRNA_mRNA_C) = names(p_val_miRNA_mRNA_IBS) = 
  names(col_p_val_C) = names(col_p_val_IBS) <- cancers_final

med_C_miRNA <- apply(C_df_cancers_EN_miRNA, 2, function(x) median(x, na.rm = T))
med_C_mRNA <- apply(C_df_cancers_EN_mRNA, 2, function(x) median(x, na.rm = T))

med_IBS_miRNA <- apply(IBS_df_cancers_EN_miRNA, 2, function(x) median(x, na.rm = T))
med_IBS_mRNA <- apply(IBS_df_cancers_EN_mRNA, 2, function(x) median(x, na.rm = T))

for(cancer in cancers_final){
  
  # C-index
  p_val_miRNA_mRNA_C[cancer] <- wilcox.test(C_df_cancers_EN_miRNA[, cancer], 
                                     C_df_cancers_EN_mRNA[, cancer], paired = T)$p.value
  if(med_C_miRNA[cancer] > med_C_mRNA[cancer]){
    col_p_val_C[cancer] <-  "royalblue3"
  }else{
    col_p_val_C[cancer] <-  "lightblue"
  }
  
  # IBS
  p_val_miRNA_mRNA_IBS[cancer] <- wilcox.test(IBS_df_cancers_EN_miRNA[, cancer], 
                                         IBS_df_cancers_EN_mRNA[, cancer], paired = T)$p.value
  if(med_IBS_miRNA[cancer] < med_IBS_mRNA[cancer]){
    col_p_val_IBS[cancer] <-  "royalblue3"
  }else{
    col_p_val_IBS[cancer] <-  "lightblue"
  }  
}

p_val_miRNA_mRNA_C_BH <- p.adjust(p_val_miRNA_mRNA_C, method = "BH")
stars_p_val_miRNA_mRNA_C_BH <- sapply(p_val_miRNA_mRNA_C_BH, my_stars.pval) 

p_val_miRNA_mRNA_IBS_BH <- p.adjust(p_val_miRNA_mRNA_IBS, method = "BH")
stars_p_val_miRNA_mRNA_IBS_BH <- sapply(p_val_miRNA_mRNA_IBS_BH, my_stars.pval)

# data ---
C_df_ggplot <- rbind(cbind(melt(C_df_cancers_EN_miRNA), data_type = "miRNA") ,
                     cbind(melt(C_df_cancers_EN_mRNA), data_type = "mRNA") )

IBS_df_ggplot <- rbind(cbind(melt(IBS_df_cancers_EN_miRNA), data_type = "miRNA") ,
                     cbind(melt(IBS_df_cancers_EN_mRNA), data_type = "mRNA") )

# plot ---
C_all_cancers_plot <- ggplot() + 
  geom_boxplot(data = C_df_ggplot, aes(x = factor(variable, levels = cancers_final), 
                                       y = value, fill = data_type), position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("C-index") + ylim(NA, 1) +
  theme_Publication_legend_right() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("#386cb0", "lightblue")) +
  geom_text(aes(x = 1:length(cancers_final), y = rep(1, length(cancers_final)), 
                label = stars_p_val_miRNA_mRNA_C_BH), col = col_p_val_C, size = 4)
print(C_all_cancers_plot)

ggsave(C_all_cancers_plot, filename = "pdf/C_miRNA_mRNA_EN.pdf")
print("Figure saved in 'pdf/C_miRNA_mRNA_EN.pdf'")


IBS_all_cancers_plot <- ggplot() + 
  geom_boxplot(data = IBS_df_ggplot, aes(x = factor(variable, levels = cancers_final), 
                                       y = value, fill = data_type), position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("IBS") + ylim(0, 0.5) +
  theme_Publication_legend_right() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("#386cb0", "lightblue")) +
  geom_text(aes(x = 1:length(cancers_final), y = rep(0.5, length(cancers_final)), 
                label = stars_p_val_miRNA_mRNA_IBS_BH), col = col_p_val_IBS, size = 4)
print(IBS_all_cancers_plot)

ggsave(IBS_all_cancers_plot, filename = "pdf/IBS_miRNA_mRNA_EN.pdf")
print("Figure saved in 'pdf/IBS_miRNA_mRNA_EN.pdf'")


# plot - RF --------------------------------------------------------------

# p-values ---
p_val_miRNA_mRNA_C = p_val_miRNA_mRNA_IBS <- rep(NA, length(cancers_final))
col_p_val_C = col_p_val_IBS <- rep(NA, length(cancers_final))
names(p_val_miRNA_mRNA_C) = names(p_val_miRNA_mRNA_IBS) = 
  names(col_p_val_C) = names(col_p_val_IBS) <- cancers_final

med_C_miRNA <- apply(C_df_cancers_RF_miRNA, 2, function(x) median(x, na.rm = T))
med_C_mRNA <- apply(C_df_cancers_RF_mRNA, 2, function(x) median(x, na.rm = T))
median(med_C_mRNA - med_C_miRNA)

med_IBS_miRNA <- apply(IBS_df_cancers_RF_miRNA, 2, function(x) median(x, na.rm = T))
med_IBS_mRNA <- apply(IBS_df_cancers_RF_mRNA, 2, function(x) median(x, na.rm = T))

for(cancer in cancers_final){
  
  # C-index
  p_val_miRNA_mRNA_C[cancer] <- wilcox.test(C_df_cancers_RF_miRNA[, cancer], 
                                            C_df_cancers_RF_mRNA[, cancer], paired = T)$p.value
  if(med_C_miRNA[cancer] > med_C_mRNA[cancer]){
    col_p_val_C[cancer] <-  "orange2"
  }else{
    col_p_val_C[cancer] <-  "yellow2"
  }
  
  # IBS
  p_val_miRNA_mRNA_IBS[cancer] <- wilcox.test(IBS_df_cancers_RF_miRNA[, cancer], 
                                              IBS_df_cancers_RF_mRNA[, cancer], paired = T)$p.value
  if(med_IBS_miRNA[cancer] < med_IBS_mRNA[cancer]){
    col_p_val_IBS[cancer] <-  "orange2"
  }else{
    col_p_val_IBS[cancer] <-  "yellow2"
  }  
}

p_val_miRNA_mRNA_C_BH <- p.adjust(p_val_miRNA_mRNA_C, method = "BH")
stars_p_val_miRNA_mRNA_C_BH <- sapply(p_val_miRNA_mRNA_C_BH, my_stars.pval) 

p_val_miRNA_mRNA_IBS_BH <- p.adjust(p_val_miRNA_mRNA_IBS, method = "BH")
stars_p_val_miRNA_mRNA_IBS_BH <- sapply(p_val_miRNA_mRNA_IBS_BH, my_stars.pval)

# data ---
C_df_ggplot <- rbind(cbind(melt(C_df_cancers_RF_miRNA), data_type = "miRNA") ,
                     cbind(melt(C_df_cancers_RF_mRNA), data_type = "mRNA") )

IBS_df_ggplot <- rbind(cbind(melt(IBS_df_cancers_RF_miRNA), data_type = "miRNA") ,
                       cbind(melt(IBS_df_cancers_RF_mRNA), data_type = "mRNA") )

# plot ---
C_all_cancers_plot <- ggplot() + 
  geom_boxplot(data = C_df_ggplot, aes(x = factor(variable, levels = cancers_final), 
                                       y = value, fill = data_type), position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("C-index") + ylim(NA, 1) +
  theme_Publication_legend_right() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("orange", "yellow")) +
  geom_text(aes(x = 1:length(cancers_final), y = rep(1, length(cancers_final)), 
                label = stars_p_val_miRNA_mRNA_C_BH), col = col_p_val_C, size = 6)
print(C_all_cancers_plot)

ggsave(C_all_cancers_plot, filename = "pdf/C_miRNA_mRNA_RF.pdf")
print("Figure saved in 'pdf/C_miRNA_mRNA_RF.pdf'")

IBS_all_cancers_plot <- ggplot() + 
  geom_boxplot(data = IBS_df_ggplot, aes(x = factor(variable, levels = cancers_final), 
                                         y = value, fill = data_type), position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("IBS") + ylim(0, 0.35) +
  theme_Publication_legend_right() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("orange", "yellow")) +
  geom_text(aes(x = 1:length(cancers_final), y = rep(0.35, length(cancers_final)), 
                label = stars_p_val_miRNA_mRNA_IBS_BH), col = col_p_val_IBS, size = 6)
print(IBS_all_cancers_plot)

ggsave(IBS_all_cancers_plot, filename = "pdf/IBS_miRNA_mRNA_RF.pdf")
print("Figure saved in 'pdf/IBS_miRNA_mRNA_RF.pdf'")
