# Load the edgeR library
library(edgeR)

# Read in the count data matrix from the CSV file
counts <- read.csv("merged_all_main_matrix__NB_CB.csv", header = TRUE, row.names = 1)

# Determine the number of genes in the count matrix
n_genes <- nrow(counts)

# Print the number of genes in the count matrix to the console
cat("Number of genes in the count matrix:", n_genes, "\n")
