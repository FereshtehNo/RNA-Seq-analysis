# Read in the two CSV files
counts <- read.csv("merged_all_main_matrix__NB_CB.csv", header=TRUE)
mart_export <- read.csv("gene_length-originaal.csv", header=TRUE)

# Merge the two dataframes based on the common column
merged_df <- merge(counts, mart_export, by.x="ENSEMBL_GENE_ID", by.y="Gene.stable.ID")

# Save the merged dataframe as a new CSV file
write.csv(merged_df, file="merged_matrix-main_with-gene_length.csv", row.names=FALSE)