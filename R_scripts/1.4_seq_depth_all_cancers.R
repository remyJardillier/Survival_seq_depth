

# median seq depth betwenn mRNA and miRNA ---------------------------------

load(file = "data_fit/seq_depth.RData")

col_cancers <- rep(NA, length(cancers_all))
names(col_cancers) <- cancers_all
col_cancers[cancers_final] <- "red"

med_SD_plot <- ggplot(med_seq_depth, aes(x=miRNA, y=mRNA, col = col_cancers)) +
  geom_point() + 
  geom_text_repel(label=rownames(med_seq_depth)) +
  theme_Publication()
print(med_SD_plot)

ggsave(med_SD_plot, filename = "pdf/med_seq_depth_miRNA_mRNA.pdf")
print("Figure saved in: pdf/med_seq_depth_miRNA_mRNA.pdf")

# boxplot -----------------------------------------------------------------

# load seq depth
load(file = "data_fit/seq_depth.RData")

col_cancers <- rep("black", length(cancers_all))
names(col_cancers) <- cancers_all
col_cancers[cancers_final] <- "red"

# ggplot miRNA
lib_size_miRNA_df <- melt(lib_size_miRNA)
SD_boxplot_miRNA <- ggplot(data = lib_size_miRNA_df, aes(x = factor(L1, levels = cancers_all), 
                                     y = value)) +
  geom_boxplot(fill = "lightblue") + xlab("Cancer") + ylab("Library size (miRNA)") +
  theme_Publication() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = col_cancers, size = 12)) +
  ylim(0, 2*10^7)  #+ ggtitle("miRNA") 

print(SD_boxplot_miRNA)

ggsave(SD_boxplot_miRNA, file = "pdf/SD_boxplot_miRNA.pdf")
print("Figure saved in: pdf/SD_boxplot_miRNA.pdf")

# ggplot mRNA
lib_size_mRNA_df <- melt(lib_size_mRNA)
SD_boxplot_mRNA <- ggplot(data = lib_size_mRNA_df, aes(x = factor(L1, levels = names(sort(med_C_EN, decreasing = T))), 
                                     y = value)) +
  geom_boxplot(fill = "lightblue") + xlab("Cancer") + ylab("Library size (mRNA)") +
  theme_Publication() + theme(axis.text.x = element_text(size = 12)) +
  ylim(0, 1*10^8) 

print(SD_boxplot_mRNA)

ggsave(SD_boxplot_mRNA, file = "pdf/SD_boxplot_mRNA.pdf")
print("Figure saved in: pdf/SD_boxplot_mRNA.pdf")

