#  --- Working directories ---
wd_cancer <- paste0("../data_cancer/", cancer, "/")

count <- t(read.csv(file = paste0(wd_cancer, cancer, ".miRseq_raw_counts.txt"),
                    stringsAsFactors = F, sep = "\t", header =F))
count[1:5, 1:5]

# colnames, rownames, NA values
colnames(count) <- count[1,]
count <- count[-1,]
row.names(count) <-  count[,1]
count <- count[,-1]
count <- as.matrix(count)
class(count) <- "numeric"

count[1:5, 1:5]

# keep only primary solid tumor samples 
# (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
freq_sample_types <- as.data.frame(table(substr(row.names(count), 14, 15)))
freq_sample_types[,2] <- paste0("(", freq_sample_types[,2], ")") 
print(paste("Sample types miRNA:", paste(paste(freq_sample_types[,1], freq_sample_types[,2]), collapse = ", ")))

if(cancer == "LAML"){
  id_prim_tumor <- as.numeric(substr(row.names(count), 14, 15)) == 3
}else{
  id_prim_tumor <- as.numeric(substr(row.names(count), 14, 15)) == 1
}

count <- count[id_prim_tumor, ]
row.names(count) <- substr(row.names(count),1,12)
dim(count)
count[1:5, 1:5]

