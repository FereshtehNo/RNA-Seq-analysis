# Filter the data based on significance thresholds and retrieve FDR and gene length
up_regulated_genes <- count_data[which(count_data$FDR < 0.05 & count_data$LogFC > 1.5), c("ENSEMBL_GENE_ID", "FDR", "Gene.length")]
down_regulated_genes <- count_data[which(count_data$FDR < 0.05 & count_data$LogFC < -1.5), c("ENSEMBL_GENE_ID", "FDR", "Gene.length")]

# Save the gene lists to CSV files
write.csv(up_regulated_genes, file = "up_regulated_genes.csv", row.names = FALSE)
write.csv(down_regulated_genes, file = "down_regulated_genes.csv", row.names = FALSE)

# Print a message indicating that the files have been saved
cat("Up-regulated genes saved to up_regulated_genes.csv\n")
cat("Down-regulated genes saved to down_regulated_genes.csv\n")