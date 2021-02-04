
# load the median sequencing depth and the number of patients
load("data_fit/seq_depth.RData")
med_seq_depth
n_patients

# Sequencing depth --------------------------------------------------------

# miRNA ---
sum_up_seq_miRNA_EN <- sum_up_miRNA_func_seq("EN", "both", signif_level)
sum_up_seq_miRNA_RF <-sum_up_miRNA_func_seq("RF", "both", signif_level)

write.xlsx(sum_up_seq_miRNA_EN, file = "tables/sum_up_seq_miRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_seq_miRNA_RF, file = "tables/sum_up_seq_miRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/sum_up_seq_miRNA_EN.xlsx")
print("Table saved in: tables/sum_up_seq_miRNA_RF.xlsx")

# mRNA ---
sum_up_seq_mRNA_EN <-sum_up_mRNA_func_seq("EN", "both", signif_level)
sum_up_seq_mRNA_RF <-sum_up_mRNA_func_seq("RF", "both", signif_level)

write.xlsx(sum_up_seq_mRNA_EN, file = "tables/sum_up_seq_mRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_seq_mRNA_RF, file = "tables/sum_up_seq_mRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/sum_up_seq_mRNA_EN.xlsx")
print("Table saved in: tables/sum_up_seq_mRNA_RF.xlsx")

# Patients in the training dataset ----------------------------------------

# miRNA ---
sum_up_patients_miRNA_EN <- sum_up_miRNA_func_patients("EN", "both", signif_level)
sum_up_patients_miRNA_RF <-sum_up_miRNA_func_patients("RF", "both", signif_level)

write.xlsx(sum_up_patients_miRNA_EN, file = "tables/sum_up_patients_miRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_patients_miRNA_RF, file = "tables/sum_up_patients_miRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/patients_sum_up_miRNA_EN.xlsx")
print("Table saved in: tables/patients_sum_up_miRNA_RF.xlsx")

# mRNA ---
sum_up_patients_mRNA_EN <-sum_up_mRNA_func_patients("EN", "both", signif_level)
sum_up_patients_mRNA_RF <-sum_up_mRNA_func_patients("RF", "both", signif_level)

write.xlsx(sum_up_patients_mRNA_EN, file = "tables/sum_up_patients_mRNA_EN.xlsx", rowNames = T)
write.xlsx(sum_up_patients_mRNA_RF, file = "tables/sum_up_patients_mRNA_RF.xlsx", rowNames = T)

print("Table saved in: tables/patients_sum_up_mRNA_EN.xlsx")
print("Table saved in: tables/patients_sum_up_mRNA_RF.xlsx")
