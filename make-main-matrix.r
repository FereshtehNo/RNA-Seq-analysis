library(dplyr)
CB_data = read.csv('GSE190496_Gene_counts_CBdata.csv')
colnames(CB_data) = c(CB_data[1, 1:5], colnames(CB_data)[6:13])
CB_data = CB_data[-1, ]
CB_data = CB_data[, -14]
NB_data = read.csv('GSE190496_Gene_counts_NBdata.csv')
merged_matrix = merge(NB_data, CB_data, by='COUNTS')
write.csv(merged_matrix, "merged_all_main_matrix__NB_CB.csv")


