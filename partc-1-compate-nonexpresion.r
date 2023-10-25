# Read the CSV file
counts <- read.csv("Merge-count-matrix_SRR17172481-SRR17172486.csv", header = TRUE, row.names = 1)

# Count the number of genes with zero counts in the control sample
zero_counts_control <- sum(counts[, "Counts.81"] == 0)

# Print the number of genes with zero counts in the control sample
cat("Number of genes with zero counts in the control sample:", zero_counts_control, "\n")

# Count the number of genes with zero counts in the COVID sample
zero_counts_covid <- sum(counts[, "Counts.86"] == 0)

# Print the number of genes with zero counts in the COVID sample
cat("Number of genes with zero counts in the COVID sample:", zero_counts_covid, "\n")