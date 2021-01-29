
# Load the packages needed ------------------------------------------------

source(file = "load_data/load_packages.R")

# load functions
invisible(lapply(list.files(path = "functions/", pattern = "[.]R$", 
                            recursive = TRUE, full.names = T), source))

# list of cancers
cancers_all <- list.files("../RData_cancer_seqDepth",  )
cancers_all <- tools::file_path_sans_ext(cancers_all)

# no miRNA data available for GBM
cancers_all <- cancers_all[! cancers_all %in% "GBM"]

# parameters ---
# number of repetitions (10 in the paper)
n_rep <- 1

# number of folds (5 in the paper)
K_folds <- 3

# parameter of the binomial law for subsampling
proportions <- c(10^-5, 0.0001, 0.001, 0.01, 0.1, 1)

# percentage of patients in the training set
perc_train <- seq(0.2, 0.8, by=0.2) # seq(0.1, 0.8, by=0.1) in the paper 

# minimum percentage of patients for which clinical data (age, gender, grade...) 
# is available
perc_min_clin <- 90

# Create folders ----------------------------------------------------------

# folders for the fit
for(cancer in cancers_all){
  if (!dir.exists(paste0("data_fit/", cancer)))
    dir.create(paste0("data_fit/", cancer), recursive = T)
}

# folders for the figures, pdf, and tables
for(cancer in cancers_all){
  if (!dir.exists("Figures"))
    dir.create("Figures", recursive = T)
  
  if (!dir.exists("tables"))
    dir.create("tables", recursive = T)
  
  if (!dir.exists("pdf"))
    dir.create("pdf", recursive = T)
}


# p-value - choose cancers ---------------------------------------------------------------

# learn models and choose cancers ---
# cancers retained: for wich the median C-index is above 0.6 for Cox or random forest
learn_new_models <- T # choose if new models have to be learned

signif_level <- 0.5 # 0.01 in the paper
source(file = "1.1_choose_cancer_models.R")

print(cancers_final)


# plot ---
source(file = "1.2_choose_cancer_plot.R")

# cancer characteristics
source(file = "1.3_ch_all_cancers.R")

# Compare seq depth of miRNA and mRNA ---
source(file = "1.4_seq_depth_all_cancers.R")


# Clinical variables ------------------------------------------------------

# clinical data available for the cancers choosed
source(file = "2.1_clin_available.R")

# Added value of miRNA data to clinical data - Cox-EN
learn_new_models <- T # choose if new models have to be learned
source(file = "2.2_clin_vs_clin_miRNA_EN.R")

# Added value of miRNA data to clinical data - RF
learn_new_models <- T # choose if new models have to be learned
source(file = "2.3_clin_vs_clin_miRNA_EN.R")


# Prediction obtained after different degradation parameters --------------

# Degradation of the data in terms of sequencing depth and number of patients ---

# Cox-EN !!! the following scrip can be time demanding !!!
learn_new_models <- T # choose if new models have to be learned
source(file = "3.1_degrade_data_EN.R")

# RF !!! the following scrip can be time demanding !!!
learn_new_models <- T # choose if new models have to be learned
source(file = "3.1_degrade_data_RF.R")

# main plot for all cancers
# !!! errors may occur due to the low number of iterations !!!
source(file = "3.3_main_plot_all_cancers.R")

# sum up of the results and little degradation for the cancers for which 
# it is impossible to degrade by a factor 2 or 5
learn_new_models <- T # choose if new models have to be learned (new degradation 
                      # factors: 2 and 5)

signif_level <- 0.5 # 0.05 in the article

source(file = "3.4_sum_up_deg_tables")

# comparison of prediction obtanined with miRNA-seq and mRNA-seq data
source(file = "3.5_comp_pred_miRNA_mRNA")


# Why prediction are degraded ? -------------------------------------------

# number of genes detected after dowsampling of the data
#   - the most expressed genes are detected
#   - the lower the sequencing depth, the lower the number of genes detected
all_cancers <- T # generate plot for all cancers ?
source(file = "4.1_n_genes_detected.R")

# Test if the degradation of the prediction is due to the number of genes detected: 
# - green : prediction obtained with the subsampled data (10^-4)
# - yellow : prediction obtained with the same genes, but (~200) but without downsampling
# - blue : prediction obtained with all the genes without downsampling
learn_new_models <- F
source(file = "4.2_why_C_degraded.R")

learn_new_models <- F
source(file = "4.3_why_C_degraded_RF.R")