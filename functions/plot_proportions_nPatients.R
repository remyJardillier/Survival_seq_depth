pal <- wes_palette("Zissou1", 100, type = "continuous")

my_stars.pval <- function(p_val){
  
  if(p_val <= 0.001){ 
    return("***")
  }else if(p_val <= 0.01 & p_val > 0.001){
    #return("\U1F7B6\U1F7B6")
    return("**")
  }else if(p_val <= 0.05 & p_val > 0.01){
    return("*")
  }else if(p_val <= 0.1 & p_val > 0.05){
    return("+")
  }else{
    return("n.s.")
  }
}

# function to plot the main figure (Fig. 2)

my_plot <- function(df_final, df_final_mRNA, n_genes_df_final, count, count_mRNA, title = NULL,
                    y_pval_mRNA, y_pval_miRNA, y_n_pat, main = "", ylim = c(0.5, 0.8), legend_tick, na.rm = F){
  
  # grid ---
  if(na.rm){
    df_grid <- melt(apply(df_final, c(2, 3), function(x) median(x, na.rm = na.rm)))
    n_genes_grid <- melt(apply(n_genes_df_final, c(2, 3), function(x) median(x, na.rm = na.rm)))
  }else{
    df_grid <- melt(apply(df_final, c(2, 3), median))
    n_genes_grid <- melt(apply(n_genes_df_final, c(2, 3), median))
  }
  df_grid$Var2 <- as.character(as.numeric(df_grid$Var2) * 100)
  n_genes_grid$value <- round(n_genes_grid$value)
  n_genes_grid$Var2 <- as.character(as.numeric(n_genes_grid$Var2) * 100)
  
  df_grid[, 2] <- as.character(df_grid[, 2])
  df_grid[,1] <- factor(df_grid[,1], levels = dimnames(df_final)[[2]])
  
  n_genes_grid[, 2] <- as.character(n_genes_grid[, 2])
  n_genes_grid[,1] <- factor(n_genes_grid[,1], levels = dimnames(df_final)[[2]])
  
  grid_ggplot <- ggplot() + 
    geom_tile(data = df_grid, aes(x=Var1, y=Var2, fill=value))  + 
    xlab("Fold reduction") + 
    ylab("Perc. of patients (%)") + theme_Publication_legend_bottom() +
    scale_fill_gradientn(colours = pal, name = main,
                         breaks = legend_tick) +
    #geom_text(data = n_genes_grid, aes(x=Var1, y=Var2, label=value), color = "white", size = 4)  +
    geom_rect(aes(ymin = 0.5, ymax = length(unique(df_grid$Var2)) + 0.5, 
                  xmin = length(unique(df_grid$Var1)) - 0.5, 
                  xmax = length(unique(df_grid$Var1)) + 0.5), fill = NA, 
              colour = "lightgray", size = 1.5)+
    geom_rect(aes(xmin = 0.5, xmax = length(unique(df_grid$Var1)) + 0.5, 
                  ymin = length(unique(df_grid$Var2)) - 0.5, 
                  ymax = length(unique(df_grid$Var2)) + 0.5), fill = NA, 
              colour = "darkslategray", size = 1.5) +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 14)) +
    scale_x_discrete(labels= 1 / proportions)
  # print(grid_ggplot)
  
  # C-index as a function of the proportions ---
  # boxplot for miRNA
  prop_df <- df_final[, , "0.8"]
  prop_df_ggplot <- melt(prop_df)
  
  # median for mRNA
  prop_df_mRNA <- df_final_mRNA[, , "0.8"]
  prop_mRNA_med <- apply(prop_df_mRNA, 2, function(x) median(x, na.rm = T))
  prop_mRNA_med <- stack(prop_mRNA_med)
  
  # compute the one-sided p-value
  if(main == "IBS"){
    alternative <- "less"
  }else{
    alternative <- "greater"
  }
  
  p_val_prop_test <- rep(NA, ncol(prop_df) - 1)
  for(i in 1:(ncol(prop_df)-1)){
    p_val_tmp <- wilcox.test(prop_df[,ncol(prop_df)], prop_df[,i], 
                        alternative = alternative)$p.value
    
    p_val_prop_test[i] <- my_stars.pval(p_val_tmp)
  }
  
  p_val_prop_test_mRNA <- rep(NA, ncol(prop_df_mRNA) - 1)
  for(i in 1:(ncol(prop_df_mRNA)-1)){
    p_val_tmp <- wilcox.test(prop_df_mRNA[,ncol(prop_df_mRNA)], prop_df_mRNA[,i], 
                        alternative = alternative)$p.value
    
    p_val_prop_test_mRNA[i] <- my_stars.pval(p_val_tmp)
  }
  
  prop_ggplot <- ggplot() +
    geom_boxplot(data = prop_df_ggplot, aes(x = factor(Var2), y = value),
                 fill = "darkslategray", position = position_dodge(0.75)) +
    geom_point(data = prop_mRNA_med, aes(x = factor(ind), y = values, 
                                         col = "red"), shape = 4, size = 1, stroke = 2) +
    geom_line(data = prop_mRNA_med, aes(x = factor(ind), y = values, group = 1), 
                                        col = "red", size = 1) +
    geom_text(aes(x = 1:length(p_val_prop_test_mRNA), y = rep(y_pval_mRNA, length(p_val_prop_test_mRNA)), 
                  label = p_val_prop_test_mRNA), col = "red", size = 5) +
    geom_text(aes(x = 1:length(p_val_prop_test), y = rep(y_pval_miRNA, length(p_val_prop_test)), 
                  label = p_val_prop_test), size = 5) +
    theme_Publication() + xlab("Fold reduction") + ylab(main) + ylim(ylim) +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 15))  +
    scale_x_discrete(labels= 1 / proportions)
  
  # print(prop_ggplot)
  
  # C-index as a function of the percentage of patients ---
  # boxplot for miRNA
  pat_df <- df_final[, "1", ]
  colnames(pat_df) <- as.character(as.numeric(colnames(pat_df)) * 100)
  pat_df_ggplot <- melt(pat_df)
  
  # median for mRNA
  pat_df_mRNA <- df_final_mRNA[, "1", ]
  colnames(pat_df_mRNA) <- as.character(as.numeric(colnames(pat_df_mRNA)) * 100)
  pat_mRNA_med <- apply(pat_df_mRNA, 2, function(x) median(x, na.rm = T))
  pat_mRNA_med <- stack(pat_mRNA_med)
  
  # compute the one-sided p-value
  p_val_pat_test <- rep(NA, ncol(pat_df) - 1)
  for(i in 1:(ncol(pat_df)-1)){
    p_val_tmp <- wilcox.test(pat_df[,ncol(pat_df)], pat_df[,i], 
                             alternative = alternative)$p.value
    
    p_val_pat_test[i] <- my_stars.pval(p_val_tmp)
  }
  
  p_val_pat_test_mRNA <- rep(NA, ncol(pat_df_mRNA) - 1)
  for(i in 1:(ncol(pat_df_mRNA)-1)){
    p_val_tmp <- wilcox.test(pat_df_mRNA[,ncol(pat_df_mRNA)], pat_df_mRNA[,i], 
                             alternative = alternative)$p.value
    
    p_val_pat_test_mRNA[i] <- my_stars.pval(p_val_tmp)
  }
  
  # number of patients
  n_pat_vect <- rep(NA, 3)
  ind <- 1
  for(i in c(0.1, 0.4, 0.8)){
    n_pat_vect[ind] <- round(i*nrow(count))
    ind <- ind + 1
  }
  
  pat_ggplot <- ggplot() +
    geom_boxplot(data = pat_df_ggplot, aes(x = factor(Var2), y = value),
                 fill = "lightgray", position = position_dodge(0.75)) +
    geom_point(data = pat_mRNA_med, aes(x = factor(ind), y = values, 
                                         col = "red"), shape = 4, size = 1, stroke = 2) +
    geom_line(data = pat_mRNA_med, aes(x = factor(ind), y = values, group = 1), 
              col = "red", size = 1) +
    geom_text(aes(x = 1:length(p_val_pat_test_mRNA), y = rep(y_pval_mRNA, length(p_val_pat_test_mRNA)), 
                  label = p_val_pat_test_mRNA), col = "red", size = 5) +
    geom_text(aes(x = 1:length(p_val_pat_test), y = rep(y_pval_miRNA, length(p_val_pat_test)), 
                  label = p_val_pat_test), size = 5) +
    geom_text(aes(x = c(1, 4, 8), y = rep(y_n_pat, length(n_pat_vect)), 
                  label = n_pat_vect), size = 5) +
    theme_Publication() + xlab("Perc. of patients (%)") + ylab(main) + ylim(ylim) +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 15)) 
  
  # print(pat_ggplot)
  
  figure <- ggarrange(grid_ggplot,                                                 
            ggarrange(prop_ggplot, pat_ggplot, nrow = 2, labels = c("B", "C")), 
            ncol = 2, 
            labels = "A"
  ) 
  if(!is.null(title)){
    figure <- annotate_figure(figure,
                              top = text_grob(title, face = "bold"))
  }
  
  return(figure)
}




# boxplot - proportion ----------------------------------------------------

# function to plot only Figure 2B
prop_boxplot <- function(df_final, df_final_mRNA, n_genes_df_final, count, count_mRNA, title = NULL,
                    y_pval_mRNA, y_pval_miRNA, y_n_pat, main = "", ylim = c(0.5, 0.8), legend_tick, na.rm = F){
  
  # C-index as a function of the proportions ---
  # boxplot for miRNA
  prop_df <- df_final[, , "0.8"]
  prop_df_ggplot <- melt(prop_df)
  
  # median for mRNA
  prop_df_mRNA <- df_final_mRNA[, , "0.8"]
  prop_mRNA_med <- apply(prop_df_mRNA, 2, function(x) median(x, na.rm = T))
  prop_mRNA_med <- stack(prop_mRNA_med)
  
  # compute the one-sided p-value
  if(main == "IBS"){
    alternative <- "less"
  }else{
    alternative <- "greater"
  }
  
  p_val_prop_test <- rep(NA, ncol(prop_df) - 1)
  for(i in 1:(ncol(prop_df)-1)){
    p_val_tmp <- wilcox.test(prop_df[,ncol(prop_df)], prop_df[,i], 
                             alternative = alternative)$p.value
    
    p_val_prop_test[i] <- my_stars.pval(p_val_tmp)
  }
  
  p_val_prop_test_mRNA <- rep(NA, ncol(prop_df_mRNA) - 1)
  for(i in 1:(ncol(prop_df_mRNA)-1)){
    p_val_tmp <- wilcox.test(prop_df_mRNA[,ncol(prop_df_mRNA)], prop_df_mRNA[,i], 
                             alternative = alternative)$p.value
    
    p_val_prop_test_mRNA[i] <- my_stars.pval(p_val_tmp)
  }
  
  prop_ggplot <- ggplot() +
    geom_boxplot(data = prop_df_ggplot, aes(x = factor(Var2), y = value),
                 fill = "darkslategray", position = position_dodge(0.75)) +
    geom_point(data = prop_mRNA_med, aes(x = factor(ind), y = values, 
                                         col = "red"), shape = 4, size = 1, stroke = 2) +
    geom_line(data = prop_mRNA_med, aes(x = factor(ind), y = values, group = 1), 
              col = "red", size = 1) +
    geom_text(aes(x = 1:length(p_val_prop_test_mRNA), y = rep(y_pval_mRNA, length(p_val_prop_test_mRNA)), 
                  label = p_val_prop_test_mRNA), col = "red", size = 5) +
    geom_text(aes(x = 1:length(p_val_prop_test), y = rep(y_pval_miRNA, length(p_val_prop_test)), 
                  label = p_val_prop_test), size = 5) +
    theme_Publication() + xlab("Fold reduction") + ylab(main) + ylim(ylim) +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 15))  +
    scale_x_discrete(labels= 1 / proportions) + ggtitle(cancer)
  
  
  return(prop_ggplot)
}


