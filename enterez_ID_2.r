# Load data frames from CSV files
df1 <- read.csv("up_regulated_genes.csv")
df2 <- read.csv("down_regulated_genes.csv")

# Connect to BioMart database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping between ENSEMBL_GENE_ID and Entrez_ID
mapping <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), 
                 filters="ensembl_gene_id", 
                 values=unique(c(df1$ENSEMBL_GENE_ID, df2$ENSEMBL_GENE_ID)),
                 mart=ensembl)

# Merge mapping with cleaned data frames
df1 <- merge(df1, mapping, by.x="ENSEMBL_GENE_ID", by.y="ensembl_gene_id", all.x=TRUE)
df2 <- merge(df2, mapping, by.x="ENSEMBL_GENE_ID", by.y="ensembl_gene_id", all.x=TRUE)

# Rename Entrez_ID column
colnames(df1)[ncol(df1)] <- "Entrez_ID"
colnames(df2)[ncol(df2)] <- "Entrez_ID"

# Save updated data frames to CSV files
write.csv(df1, "clean_up_regulated_genes_entrez-2.csv", row.names=FALSE)
write.csv(df2, "clean_down_regulated_genes_entrez-2.csv", row.names=FALSE)