# p-values - compare distributions EN - RF --------------------------------

load(file = "data_fit/pred_ch_cancers.RData")

# C-index
p_val_EN_RF <- rep(NA, length(cancers_all))
names(p_val_EN_RF) <- names(sort(med_C_EN, decreasing = T))

col_p_val <- rep(NA, length(cancers_all))
names(col_p_val) <- names(sort(med_C_EN, decreasing = T))

# IBS
p_val_EN_RF_IBS <- rep(NA, length(cancers_all))
names(p_val_EN_RF_IBS) <- names(sort(med_C_EN, decreasing = T))

col_p_val_IBS <- rep(NA, length(cancers_all))
names(col_p_val_IBS) <- names(sort(med_C_EN, decreasing = T))

med_IBS_EN <- apply(IBS_df_cancers_EN, 2, function(x) median(x, na.rm = T))
med_IBS_RF <- apply(IBS_df_cancers_RF, 2, function(x) median(x, na.rm = T))

for(cancer in cancers_all){
  
  # C-index
  p_val_EN_RF[cancer] <- wilcox.test(C_df_cancers_EN[, cancer], 
                                     C_df_cancers_RF[, cancer], paired = T)$p.value
  if(med_C_EN[cancer] > med_C_RF[cancer]){
    col_p_val[cancer] <-  "royalblue3"
  }else{
    col_p_val[cancer] <-  "orange2"
  }
  
  # IBS
  p_val_EN_RF_IBS[cancer] <- wilcox.test(IBS_df_cancers_EN[, cancer], 
                                         IBS_df_cancers_RF[, cancer], paired = T)$p.value
  if(med_IBS_EN[cancer] < med_IBS_RF[cancer]){
    col_p_val_IBS[cancer] <-  "royalblue3"
  }else{
    col_p_val_IBS[cancer] <-  "orange3"
  }
}

p_val_EN_RF_BH <- p.adjust(p_val_EN_RF, method = "BH")
stars_p_val_EN_RF_BH <- sapply(p_val_EN_RF_BH, my_stars.pval) 

p_val_EN_RF_IBS_BH <- p.adjust(p_val_EN_RF_IBS, method = "BH")
stars_p_val_EN_RF_IBS_BH <- sapply(p_val_EN_RF_IBS_BH, my_stars.pval) 

# plot --------------------------------------------------------------------

load(file = "data_fit/pred_ch_cancers.RData")

# C-index --- 
C_df_ggplot <- rbind(cbind(melt(C_df_cancers_RF), method = "RF") ,
                     cbind(melt(C_df_cancers_EN), method = "EN") )

col_cancer_name <- rep("black", length(cancers_all))
names(col_cancer_name) <- names(sort(med_C_EN, decreasing = T))
col_cancer_name[cancers_final] <- "red"

C_all_cancers_plot <- ggplot() + 
  geom_boxplot(data = C_df_ggplot, aes(x = factor(variable, levels = names(sort(med_C_EN, decreasing = T))), 
                                       y = value, fill = method), position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("C-index") + ylim(NA, 1) +
  theme_Publication_legend_right() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("#386cb0", "orange")) +
  geom_abline(slope = 0, intercept = 0.6, col = "red", lwd = 1, lty = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = col_cancer_name, size = 12)) +
  geom_text(aes(x = 1:length(cancers_all), y = rep(1, length(cancers_all)), 
                label = stars_p_val_EN_RF_BH), col = col_p_val, size = 4)
print(C_all_cancers_plot)

ggsave(C_all_cancers_plot, filename = "pdf/C_all_cancers_miRNA.pdf")
print("Figure saved in: pdf/C_all_cancers_miRNA.pdf")

# IBS ---
IBS_df_ggplot <- rbind(cbind(melt(IBS_df_cancers_RF), method = "RF") ,
                       cbind(melt(IBS_df_cancers_EN), method = "EN") )


IBS_all_cancers_plot <- ggplot() + 
  geom_boxplot(data = IBS_df_ggplot, aes(x = factor(variable, levels = names(sort(med_C_EN, decreasing = T))), 
                                         y = value, fill = method), position = position_dodge(0.75)) + 
  xlab("Cancer") + ylab("IBS") + 
  theme_Publication_legend_right() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("#386cb0", "orange")) + ylim(NA, 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = col_cancer_name, size = 12)) +
  geom_text(aes(x = 1:length(cancers_all), y = rep(0.7, length(cancers_all)), 
                label = stars_p_val_EN_RF_IBS_BH), col = col_p_val_IBS, size = 4)
print(IBS_all_cancers_plot)

ggsave(IBS_all_cancers_plot, filename = "pdf/IBS_all_cancers_miRNA.pdf")
print("Figure saved in: pdf/IBS_all_cancers_miRNA.pdf")
