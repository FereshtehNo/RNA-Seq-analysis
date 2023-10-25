# Load the required library
library(edgeR)

# Read in the count matrix from the CSV file
counts <- read.csv("merged_all_main_matrix__NB_CB.csv", header = TRUE, row.names = 1)

# Replace any NA values with 0
counts[is.na(counts)] <- 0

# Write the modified count matrix to a new CSV file
write.csv(counts, file = "modify-merge-main-matrix.csv", row.names = TRUE)

# Print a message indicating that the file has been saved
cat("Modified count matrix saved to modify-merge-main-matrix.csv\n")

# Identify any non-numeric values in each column of the data frame
non_numeric_cols <- sapply(counts, function(x) any(grepl("[^0-9.]", x)))

# Exclude the non-numeric columns from the data frame
gene_expression_numeric <- counts[, !non_numeric_cols]

# Define the experimental design and create the DGEList object
group <- factor(c(rep("normal", 5), rep("covid", 8)))
dge <- DGEList(counts = as.matrix(gene_expression_numeric), group = group)

# Identify which columns correspond to normal samples and which ones correspond to COVID samples
normal_cols <- c("NB5", "NB4", "NB1", "NB2", "NB3")
covid_cols <- c("CB4", "CB5", "CB6", "CB7", "CB2", "CB3", "CB1", "CB8")

# Filter the genes that have low expression and/or low variability
keep <- rowSums(cpm(dge[, normal_cols]) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
p <- glmQLFTest(fit, coef = 2)
keep <- p$table$PValue < 0.05
dge <- dge[keep, ]

# Get the number of differentially expressed genes
n_genes <- sum(p$table$PValue < 0.05 & rowMeans(cpm(dge[, normal_cols])) > 1)
cat("Number of differentially expressed genes:", n_genes, "\n")

# Get the top differentially expressed genes
#tags <- topTags(fit, n = n_genes)
#de_genes <- rownames(tags$table)

# Print the top differentially expressed genes with their statistics
#cat("Top differentially expressed genes:\n")
#print(tags$table)