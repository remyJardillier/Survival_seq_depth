
# load the median sequencing depth and the number of patients
load("data_fit/seq_depth.RData")
med_seq_depth
n_patients

signif_level <- 0.5 # 0.05 in the paper

# Sequencing depth --------------------------------------------------------

# degradation if the two metrics are degraded 
signif_level <- 0.5 # 0.05 in the paper

# miRNA ---
sum_up_patients_miRNA_EN <- sum_up_miRNA_func_seq("EN", "both", signif_level)
sum_up_patients_miRNA_RF <-sum_up_miRNA_func_seq("RF", "both", signif_level)

write.xlsx(sum_up_patients_miRNA_EN, file = "tables/patients_sum_up_miRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_patients_miRNA_RF, file = "tables/patients_sum_up_miRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/patients_sum_up_miRNA_EN.xlsx")
print("Table saved in: tables/patients_sum_up_miRNA_RF.xlsx")

# mRNA ---
sum_up_patients_mRNA_EN <-sum_up_mRNA_func_seq("EN", "both", signif_level)
sum_up_patients_mRNA_RF <-sum_up_mRNA_func_seq("RF", "both", signif_level)

write.xlsx(sum_up_patients_mRNA_EN, file = "tables/patients_sum_up_mRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_patients_mRNA_RF, file = "tables/patients_sum_up_mRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/patients_sum_up_mRNA_EN.xlsx")
print("Table saved in: tables/patients_sum_up_mRNA_RF.xlsx")

# Patients in the training dataset ----------------------------------------

# digredation if the two metrics are degraded 

signif_level <- 0.5 # 0.05 in the paper

# miRNA ---
sum_up_patients_miRNA_EN <- sum_up_miRNA_func_patients("EN", "both", signif_level)
sum_up_patients_miRNA_RF <-sum_up_miRNA_func_patients("RF", "both", signif_level)

write.xlsx(sum_up_patients_miRNA_EN, file = "tables/patients_sum_up_miRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_patients_miRNA_RF, file = "tables/patients_sum_up_miRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/patients_sum_up_miRNA_EN.xlsx")
print("Table saved in: tables/patients_sum_up_miRNA_RF.xlsx")

# mRNA ---
sum_up_patients_mRNA_EN <-sum_up_mRNA_func_patients("EN", "both", signif_level)
sum_up_patients_mRNA_RF <-sum_up_mRNA_func_patients("RF", "both", signif_level)

write.xlsx(sum_up_patients_mRNA_EN, file = "tables/patients_sum_up_mRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_patients_mRNA_RF, file = "tables/patients_sum_up_mRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/patients_sum_up_mRNA_EN.xlsx")
print("Table saved in: tables/patients_sum_up_mRNA_RF.xlsx")