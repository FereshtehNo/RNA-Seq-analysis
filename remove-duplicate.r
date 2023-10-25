library(dplyr)

# Load data
df1 <- read.csv("up_regulated_genes.csv")
df2 <- read.csv("down_regulated_genes.csv")

# Check for duplicates in df1
df1_dup <- df1[duplicated(df1$Entrez_ID),]
if (nrow(df1_dup) > 0) {
  # Remove duplicates from df1
  df1 <- distinct(df1, Entrez_ID, .keep_all=TRUE)
  message(paste0("Removed ", nrow(df1_dup), " duplicate(s) from df1"))
}

# Check for duplicates in df2
df2_dup <- df2[duplicated(df2$Entrez_ID),]
if (nrow(df2_dup) > 0) {
  # Remove duplicates from df2
  df2 <- distinct(df2, Entrez_ID, .keep_all=TRUE)
  message(paste0("Removed ", nrow(df2_dup), " duplicate(s) from df2"))
}
# Check for duplicates in df1
if (any(duplicated(df1$Entrez_ID))) {
  message("df1 contains duplicates")
}

# Check for duplicates in df2
if (any(duplicated(df2$Entrez_ID))) {
  message("df2 contains duplicates")
}
# Save cleaned data frames to CSV files
write.csv(df1, "clean_up_regulated_genes.csv", row.names=FALSE)
write.csv(df2, "clean_down_regulated_genes.csv", row.names=FALSE)
