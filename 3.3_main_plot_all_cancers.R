
# EN ----------------------------------------------------------------------

# create temporary folders
if (!dir.exists("pdf_tmp"))
  dir.create("pdf_tmp", recursive = T)

for(cancer in cancers_final){
  
  print(paste0("Start learning for: ", cancer))
  
  source(file = "load_data/load_data_final.R")
  
  # load the files ---
  load(file = paste0("data_fit/", cancer,  "/pred_degrade_seq_depth_miRNA_EN.RData"))
  load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_EN.RData"))
  
  # plot ---
  # C-index
  ggplot_C <- my_plot(C_df_final, C_df_final_mRNA, n_genes_df_final, count, count_mRNA, y_pval_mRNA = 1,
                        y_pval_miRNA = 0.97, y_n_pat = 0.46, main = "C-index", title = "C-index", 
                        legend_tick = c(0.5, 0.6,0.7,0.8,0.9,1),
                        ylim = c(0.45, 1), na.rm = T)
  
  # IBS
  ggplot_IBS <- my_plot(IBS_df_final, IBS_df_final_mRNA, n_genes_df_final, count, count_mRNA, y_pval_mRNA = 0.35, 
                        y_pval_miRNA = 0.32, y_n_pat = 0, main = "IBS", title = "IBS",
                        legend_tick = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),
                        ylim = c(0, 0.35), na.rm = T)
  
  test <- ggarrange(ggplot_C, ggplot_IBS,
                    nrow = 2) 
  test <- annotate_figure(test, fig.lab = cancer, fig.lab.pos = "top.left")
  
  ggsave(test, filename = paste0("pdf_tmp/", cancer, "_test.pdf"), width = 210, height = 297, 
         units = "mm", device = cairo_pdf)
  
}

# build pdf files, save all the files for all cancers, and remove temporary folders
pdf_names <- list.files("pdf_tmp", full.names = T)
pdf_combine(pdf_names, output = "pdf/all_cancers_grid_EN.pdf")
unlink("pdf_tmp", recursive = T)

print("Figure saved in: pdf/all_cancers_grid_EN.pdf")


# RF ----------------------------------------------------------------------

# create temporary folders
if (!dir.exists("pdf_tmp"))
  dir.create("pdf_tmp", recursive = T)

for(cancer in cancers_final){
  
  print(paste0("Start learning for: ", cancer))
  
  source(file = "load_data/load_data_final.R")
  
  # load the files ---
  load(file = paste0("data_fit/", cancer,  "/pred_degrade_seq_depth_miRNA_RF.RData"))
  load(file = paste0("data_fit/", cancer, "/pred_degrade_seq_depth_mRNA_RF.RData"))
  
  # plot ---
  # C-index
  ggplot_C <- my_plot(C_df_final, C_df_final_mRNA, n_genes_df_final, count, count_mRNA, y_pval_mRNA = 1,
                        y_pval_miRNA = 0.97, y_n_pat = 0.46, main = "C-index", title = "C-index", 
                        legend_tick = c(0.5, 0.6,0.7,0.8,0.9,1),
                        ylim = c(0.45, 1), na.rm = T)
  
  # IBS
  ggplot_IBS <- my_plot(IBS_df_final, IBS_df_final_mRNA, n_genes_df_final, count, count_mRNA, y_pval_mRNA = 0.35, 
                        y_pval_miRNA = 0.32, y_n_pat = 0, main = "IBS", title = "IBS",
                        legend_tick = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),
                        ylim = c(0, 0.35), na.rm = T)
  
  test <- ggarrange(ggplot_C, ggplot_IBS,
                    nrow = 2) 
  test <- annotate_figure(test, fig.lab = cancer, fig.lab.pos = "top.left")
  
  ggsave(test, filename = paste0("pdf_tmp/", cancer, "_test.pdf"), width = 210, height = 297, 
         units = "mm", device = cairo_pdf)
  
}

# build pdf files, save all the files for all cancers, and remove temporary folders
pdf_names <- list.files("pdf_tmp", full.names = T)
pdf_combine(pdf_names, output = "pdf/all_cancers_grid_RF.pdf")
unlink("pdf_tmp", recursive = T)

print("Figure saved in: pdf/all_cancers_grid_RF.pdf")