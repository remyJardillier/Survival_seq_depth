
# Learn models ------------------------------------------------------------

alpha_vect <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

if(learn_new_models){
  
  for(cancer in cancers_final){
    
    print(paste("*** start learning for", cancer, "***"))
    
    # load the data ---
    source(file = "load_data/load_data_final.R")
    
    # logCPM data
    logCPM_data <- log.cpm(count)
    logCPM_data_std <- scale(logCPM_data)
    
    # Surv object for the Cox model
    y_cox <- Surv(time = clin$time, event = clin$status)
    
    # learn models ---
    dev_mat <- data.frame(matrix(ncol = 4, nrow = length(alpha_vect)))
    row.names(dev_mat) <- alpha_vect
    colnames(dev_mat) <- c("alpha", "low", "min", "up")
    dev_mat$alpha <- alpha_vect
    
    n_genes_slc = comp_time <- rep(NA, length(alpha_vect))
    names(n_genes_slc) = names(comp_time) <- alpha_vect
    
    for(i in 1:length(alpha_vect)){
      
      print(paste0("Start learning for alpha = ", alpha_vect[i]))
      
      comp_time[i] <- system.time(fit_tmp <- learn_model_EN(logCPM_data_std, clin, alpha_vect[i])) 
      
      id_min <- which.min(fit_tmp$cvm)
      dev_mat[i, "min"] <- fit_tmp$cvm[id_min]
      dev_mat[i, "low"] <- fit_tmp$cvlo[id_min]
      dev_mat[i, "up"] <- fit_tmp$cvup[id_min]
      
      n_genes_slc[i] <- length(genes_selected_func(fit_tmp, "lambda.min"))
    }
    
    save(dev_mat, comp_time, n_genes_slc,
         file = paste0("data_fit/", cancer, "/choose_alpha_EN.RData"))
    
  }
}



# Ggplot - all ------------------------------------------------------------------

if (!dir.exists("pdf_tmp"))
  dir.create("pdf_tmp", recursive = T)

for(cancer in cancers_all){
  
  load(paste0("data_fit/", cancer, "/choose_alpha_EN.RData"))
  
  # Déviance
  plot_dev <- ggplot() +
    geom_segment(data=dev_mat, aes(x = alpha, xend = alpha, y = low, yend = up), col = "darkgrey", size = 2) +
    geom_point(data=dev_mat, aes(x = alpha, y = min), col="red", size = 4) +
    theme_Publication() + xlab(expression(alpha)) + ylab("Déviance")
  
  # Nombre de gènes sélectionnés
  plot_genes <- ggplot(data = data.frame(n_genes = n_genes_slc, alpha = as.numeric(names(n_genes_slc))), 
         aes(x = alpha, y = n_genes)) +
    geom_point(col = "red", shape = 4, size = 2, stroke = 2) + theme_Publication() +
    xlab(expression(alpha)) + ylab("Nombre de gènes sélectionnés")
  
  # Temps de calcul
  plot_comp_time <- ggplot(data = data.frame(comp_time = comp_time, alpha = as.numeric(names(comp_time))), 
         aes(x = alpha, y = comp_time)) +
    geom_point(col = "red", shape = 4, size = 2, stroke = 2) + theme_Publication() +
    xlab(expression(alpha)) + ylab("Temps de calcul")
  
  plot_final <- ggarrange(plot_dev, plot_genes, plot_comp_time,
                          nrow = 3, ncol = 1, labels = c("A", "B", "C"))
  plot_final <- annotate_figure(plot_final, top = cancer)
  
  ggsave(plot_final, filename = paste0("pdf_tmp/", cancer, ".pdf"), 
         width = 297, height = 210, units = "mm")
}

pdf_names <- list.files("pdf_tmp", full.names = T)
pdf_combine(pdf_names, output = "pdf/alpha_EN.pdf")
unlink("pdf_tmp", recursive = T)

# Ggplot - one ------------------------------------------------------------------

cancer <- "KIRP"

load(paste0("data_fit/", cancer, "/choose_alpha_EN.RData"))

# Deviance
plot_dev <- ggplot() +
  geom_segment(data=dev_mat, aes(x = alpha, xend = alpha, y = low, yend = up), col = "darkgrey", size = 2) +
  geom_point(data=dev_mat, aes(x = alpha, y = min), col="red", size = 4) +
  theme_Publication() + xlab(expression(alpha)) + ylab("Deviance")
plot_dev

# Number of genes selected
plot_genes <- ggplot(data = data.frame(n_genes = n_genes_slc, alpha = as.numeric(names(n_genes_slc))), 
       aes(x = alpha, y = n_genes)) +
  geom_point(col = "red", shape = 4, size = 2, stroke = 2) + theme_Publication() +
  xlab(expression(alpha)) + ylab("Number of genes selected")
plot_genes

# computation time
plot_time <- ggplot(data = data.frame(comp_time = comp_time, alpha = as.numeric(names(comp_time))), 
       aes(x = alpha, y = comp_time)) +
  geom_point(col = "red", shape = 4, size = 2, stroke = 2) + theme_Publication() +
  xlab(expression(alpha)) + ylab("computation time")
plot_time

plot_final <- ggarrange(plot_dev, plot_genes, nrow = 1, ncol = 2,
                        labels = c("A", "B"))
plot_final

ggsave(plot_final, filename = paste0("pdf/alpha_EN_", cancer, ".pdf"))



# Choice lambda -----------------------------------------------------------

cancer <- "KIRC"

print(paste("*** start learning for", cancer, "***"))

# load the data ---
source(file = "load_data/load_data_final.R")

# logCPM data
logCPM_data <- log.cpm(count)
logCPM_data_std <- scale(logCPM_data)

# Surv object for the Cox model
y_cox <- Surv(time = clin$time, event = clin$status)

fit_EN <- learn_models(logCPM_data_std, clin, "EN")

pdf(file = "pdf/ch1_lambda_cv.pdf", width=6, height=4)
plot(fit_EN, xlim = c(-2.2,-0.5), ylim = c(12,14), ylab = "Déviance")
text(x = log(fit_EN$lambda.min), y=12, labels = expression(paste(lambda, ".min")), cex = 1)
text(x = log(fit_EN$lambda.1se), y=12, labels = expression(paste(lambda, ".1se")), cex = 1)
abline(h = fit_EN$cvup[which(fit_EN$lambda == fit_EN$lambda.min)], lty = "dashed", lwd = 2)
dev.off()
