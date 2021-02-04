# load the data
source(file = "load_data/load_data_final.R")

# number of genes detected ------------------------------------------------

n_genes_detect <- rep(NA, length(proportions))
names(n_genes_detect) <- proportions

for(i in 1:length(proportions)){
  
  count_sub <- t(generateSubsampledMatrix(counts = t(count), 
                                          proportion = proportions[i], 
                                          seed = rnorm(1)))
  
  # remove genes with NA values
  id_NA_gene <- apply(count_sub, 2, function(x) sum(which(is.na(x))))
  id_NA_gene <- id_NA_gene > 0
  sum(id_NA_gene)
  count_sub <- count_sub[, !id_NA_gene]
  
  # filtering and log-CPM normalization
  logCPM_sub <- log.cpm(count_sub)
  n_genes_detect[i] <- ncol(logCPM_sub)
}

n_genes_detect_df <- data.frame(id = 1:length(proportions), n_genes = n_genes_detect)

n_genes_detect_ggplot <- ggplot(data = n_genes_detect_df, aes(x = id, y = n_genes)) + 
  geom_point(col = "red", shape = 4, size = 1, stroke = 2) + geom_line(col = "red", size = 1) + 
  xlab("Fold reduction") + ylab("N genes detected") +
  theme_Publication() + 
  scale_x_discrete(limits = as.character(1 / proportions))



# expression level of the genes detected ----------------------------------


# degrade the data and identify the genes detected
count_sub <- t(generateSubsampledMatrix(counts = t(count), 
                                        proportion = 10^(-4), 
                                        seed = rnorm(1)))

# remove genes with NA values
id_NA_gene <- apply(count_sub, 2, function(x) sum(which(is.na(x))))
id_NA_gene <- id_NA_gene > 0
sum(id_NA_gene)
count_sub <- count_sub[, !id_NA_gene]

# filtering and log-CPM normalization
logCPM_sub <- log.cpm(count_sub)
genes_detect_sub <- colnames(logCPM_sub)

# logCPM values of all the genes and only the genes detected
logCPM_all <- log.cpm(count)
logCPM_sub <- logCPM_all[, intersect(genes_detect_sub, colnames(logCPM_all))]

logCPM_all_vect <- as.vector(logCPM_all)
names(logCPM_all_vect) <- rep("No subsampling", length(logCPM_all_vect))

logCPM_sub_vect <- as.vector(logCPM_sub)
names(logCPM_sub_vect) <- rep(paste(expression(delta), "=10000"), length(logCPM_sub_vect))

# ggplot
genes_detect_df <- rbind(stack(logCPM_all_vect), stack(logCPM_sub_vect))
head(genes_detect_df)
colnames(genes_detect_df)[2] <- "p"
 
hist_expr_level_ggplot <- ggplot(data = genes_detect_df, aes(x = values, 
                                                             col = p, 
                                                             fill = p)) + 
  geom_density(alpha = 0.6) + 
  theme_Publication_legend_bottom() + 
  xlab("Level of expression (log2-CPM)") + ylab("Density") + 
  scale_color_manual(values = c("#386cb0", "#fdb462"))+ 
  scale_fill_manual(values = c("#386cb0", "#fdb462")) +
  theme(legend.title = element_blank())


# arrange the two plots
n_genes_detect_final_plot <- ggarrange(n_genes_detect_ggplot, hist_expr_level_ggplot,
          ncol = 2, labels = c("A", "B"))
n_genes_detect_final_plot <- annotate_figure(n_genes_detect_final_plot, top = cancer)
print(n_genes_detect_final_plot)

ggsave(n_genes_detect_final_plot, filename = "pdf/n_genes_and_distrib.pdf")
print("Figure saved in: pdf/n_genes_and_distrib.pdf")


# Same for mRNA -----------------------------------------------------------

# number of genes detected ---

n_genes_detect <- rep(NA, length(proportions))
names(n_genes_detect) <- proportions

for(i in 1:length(proportions)){
  
  count_sub <- t(generateSubsampledMatrix(counts = t(count_mRNA), 
                                          proportion = proportions[i], 
                                          seed = rnorm(1)))
  # remove genes with NA values
  id_NA_gene <- apply(count_sub, 2, function(x) sum(which(is.na(x))))
  id_NA_gene <- id_NA_gene > 0
  sum(id_NA_gene)
  count_sub <- count_sub[, !id_NA_gene]
  
  # filtering and log-CPM normalization
  logCPM_sub <- log.cpm(count_sub)
  n_genes_detect[i] <- ncol(logCPM_sub)
}

n_genes_detect_df <- data.frame(id = 1:length(proportions), n_genes = n_genes_detect)

n_genes_detect_ggplot <- ggplot(data = n_genes_detect_df, aes(x = id, y = n_genes)) + 
  geom_point(col = "red", shape = 4, size = 1, stroke = 2) + geom_line(col = "red", size = 1) + 
  xlab("Fold reduction") + ylab("N genes detected") +
  theme_Publication() + 
  scale_x_discrete(limits = as.character(1 / proportions))


# expression level of the genes detected ---

# degrade the data and identify the genes detected
count_sub <- t(generateSubsampledMatrix(counts = t(count_mRNA), 
                                        proportion = 10^(-5), 
                                        seed = rnorm(1)))

# remove genes with NA values
id_NA_gene <- apply(count_sub, 2, function(x) sum(which(is.na(x))))
id_NA_gene <- id_NA_gene > 0
sum(id_NA_gene)
count_sub <- count_sub[, !id_NA_gene]

# filtering and log-CPM normalization
logCPM_sub <- log.cpm(count_sub)
genes_detect_sub <- colnames(logCPM_sub)

# logCPM values of all the genes and only the genes detected
logCPM_all <- log.cpm(count_mRNA)
logCPM_sub <- logCPM_all[, intersect(genes_detect_sub, colnames(logCPM_all))]

logCPM_all_vect <- as.vector(logCPM_all)
names(logCPM_all_vect) <- rep("No subsampling", length(logCPM_all_vect))

logCPM_sub_vect <- as.vector(logCPM_sub)
names(logCPM_sub_vect) <- rep(paste(expression(delta), "=10000"), length(logCPM_sub_vect))

# ggplot
genes_detect_df <- rbind(stack(logCPM_all_vect), stack(logCPM_sub_vect))
head(genes_detect_df)
colnames(genes_detect_df)[2] <- "p"

hist_expr_level_ggplot <- ggplot(data = genes_detect_df, aes(x = values, 
                                                             col = p, 
                                                             fill = p)) + 
  geom_density(alpha = 0.6) + 
  theme_Publication_legend_bottom() + 
  xlab("Niveau d'expression (log-CPM)") + ylab("Density") + 
  scale_color_manual(values = c("#386cb0","#fdb462"))+ 
  scale_fill_manual(values = c("#386cb0", "#fdb462")) +
  theme(legend.title = element_blank())

# arrange the two plots
n_genes_detect_final_plot <- ggarrange(n_genes_detect_ggplot, hist_expr_level_ggplot,
                                       ncol = 2, labels = c("A", "B"))
n_genes_detect_final_plot <- annotate_figure(n_genes_detect_final_plot, top = cancer)
print(n_genes_detect_final_plot)

ggsave(n_genes_detect_final_plot, filename = "pdf/n_genes_and_distrib_mRNA.pdf")
print("Figure saved in: pdf/n_genes_and_distrib_mRNA.pdf")