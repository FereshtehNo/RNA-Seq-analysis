# Replace "your_fastq_file.fastq" with the path to your fastq file
fastq_file <- file("SRR17172486-trimed.fastq", open = "r")

# Initialize variables for counting bases and reads
bases <- 0
reads <- 0

# Loop through the fastq file and count the number of bases and reads
while (length(header <- readLines(fastq_file, n = 1)) > 0) {
  sequence <- readLines(fastq_file, n = 1)
  readLines(fastq_file, n = 1)
  quality <- readLines(fastq_file, n = 1)
  bases <- bases + nchar(sequence)
  reads <- reads + 1
}

# Calculate the average read length
avg_read_length <- bases/reads

# Print the average read length
cat("Average read length:", avg_read_length, "\n")

# Close the fastq file
close(fastq_file)