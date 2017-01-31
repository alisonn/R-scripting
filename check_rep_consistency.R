# Checks consistency between the biological replicates (same genotype and treatment)
# for p-values, frequency of summer or winter alleles, readDepth consistency
# Output: 4 plots in pdf format and a .tsv with the merged replicates to avoid multiple file-handling downstream

# Input: rep1 and rep2 files - these are outputs from alisonn@github: python-scripting/find_parent_contribution.py
args = commandArgs(TRUE)
rep1_file = args[1]
rep2_file = args[2]
out_pdf = args[3]
out_file = args[4]

# Read in replicate files
rep1 <- read.table(rep1_file, stringsAsFactors = FALSE, header = TRUE)
rep2 <- read.table(rep2_file, stringsAsFactors = FALSE, header = TRUE)

# create a new column merging chr# and position to make matching easier
rep1$SNP_pos <- paste(rep1$contig, rep1$pos, sep = ":")
rep2$SNP_pos <- paste(rep2$contig, rep2$pos, sep = ":")

# create a summer allele frequency column for each (i.e. Parent = 1 = Summer)
rep1$f.summer <- ifelse(rep1$altParent == "1", rep1$altFreq, (1-rep1$altFreq))
rep2$f.summer <- ifelse(rep2$altParent == "1", rep2$altFreq, (1-rep2$altFreq))

# Merge the two data frames together based only on the new column (chr#(####))
reps_merged <- merge(x = rep1, y = rep2, by = "SNP_pos", all = TRUE)

# Filter by read depth
rd_cutoff = 100
reps_merged_filtered <- subset(reps_merged, reps_merged$readDepth.x >= rd_cutoff & reps_merged$readDepth.y >= rd_cutoff)

pdf(out_pdf, width = 10, height = 10)
par(mfrow = c(2,2))

# Plot log_10(p-value) of the replicates against one another to check for consistency
# Are the p-values of ASE at each particular site consistent within the biological replicates?
plot(log(reps_merged_filtered$p.value.x, base = 10) ~ log(reps_merged_filtered$p.value.y, base = 10), main = "SNP ASE P-value Replicate Consistency", ylab = "Replicate 1 Log10(SNP P-value)", xlab = "Replicate 2 Log10(SNP P-value)", ylim = c(-200, 0), xlim = c(-200, 0), cex = 0.6)

plot(log(reps_merged_filtered$p.value.x, base = 10) ~ log(reps_merged_filtered$p.value.y, base = 10), main = "SNP ASE P-value Replicate Consistency", ylab = "Replicate 1 Log10(SNP P-value)", xlab = "Replicate 2 Log10(SNP P-value)", ylim = c(-15, 0), xlim = c(-15, 0), cex = 0.6)

# Plot summer allele frequency of the replicates against one another to check for consistency
# Are the allele frequencies at each SNP position consistent within the biological replicatesl
plot(reps_merged_filtered$f.summer.x ~ reps_merged_filtered$f.summer.y, main = "Summer Allele Freq Replicate Consistency", ylab = "Replicate 1 Summer Allele Freq", xlab = "Replicate 2 Summer Allele Freq", cex = 0.6) 

# plot readdepth at each snp site
plot(reps_merged_filtered$readDepth.x ~ reps_merged_filtered$readDepth.y, ylim = c(100, 1000), xlim = c(100,1000), xlab = "Replicate 2 Read Depth", ylab = "Replicate 1 Read Depth", cex = 0.6, main = "Read Depth Consistency")
dev.off()

# Simplified dataframe for further downstream processing
out_dataframe = reps_merged[, c(1, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16, 9, 17)]

# write the data.frame
#x is rep1, y is rep2, keep all sites (no filters)
write.table(out_dataframe, file = out_file, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE) 
