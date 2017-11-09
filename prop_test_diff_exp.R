#!/bin/env R

# A test of proportions using differential expression data at two levels:
# 1. per gene basis and 2. per overall sample basis

# Load in the two differential expression files as df
df1 <- read.table(args[1], stringsAsFactors=F, header=T)
df2 <- read.table(args[2], stringsAsFactors=F, header=T)

out_file <- args[3]

# merge by column values (namely Fbgn###)
temp_df <- merge(df1, df2, by="gene_id")

#calculate totals for shared genes from cross3 and cross4 to prepare for prop.test()
temp_df$total.x <- temp_df$value_1.x + temp_df$value_2.x
temp_df$total.y <- temp_df$value_1.y + temp_df$value_2.y

#filter for totals > 0 because prop.test() needs n > 0
merged_df <- subset(temp_df, temp_df$total.x > 0 & temp_df$total.y > 0)

# cannot have gene_id in the data frame for prop.test so must append afterward and preserve indexing
# apply(df, by row, perform prop.test by row)
final_df <- merge_df[,c(8,21,28,29)]
final_df$prop_test_pval <- 1.00

# calculates if the proportions of upregulated genes per gene are equivalent in both crosses (merged biological replicates)
final_df$prop_test_pval <- apply(final_df, 1, function(row) prop.test(x=c(row[1], row[2]), n=c(row[3],row[4]))$p.value)
final_df$gene_id <- merge_df$gene_id
final_df$log2.fold_change.x <- merge_df$log2.fold_change..x
final_df$log2.fold_change.y <- merge_df$log2.fold_change..y

hist(final_df$prop_test_pval, breaks=50, main = "P-value Distribution of Prop.Test(): Up-regulated genes across Cross 3 & 4", ylab = "Number of Genes", xlab = "P-value from Prop.Test() to reject NULL")

## Write table with Cross-Level: genes consistent in transcripts
write.table(final_df, file=out_file, append=FALSE, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

## want to compare $log2.fold_change
merged_df2 <- subset(merged_df, merged_df$status.x == "OK" & merged_df$status.y == "OK")

total_nrow <- nrow(merged_df2)
success.x <- nrow(subset(merged_df2, merged_df2$log2.fold_change..x > 0))
success.y <- nrow(subset(merged_df2, merged_df2$log2.fold_change..y > 0))

# get statistic for either rejecting or not rejecting alt hypothesis that the two crosses are different
object <- prop.test(x=c(success.x, success.y), n=c(total_nrow,total_nrow))
## object$p.value

## want the indices of consistently expressed genes
indices <- which(final_df$prop_test_pval <= 0.05)
## get all original information of the consistently expressed genes
temp_df3 <- merge_df[indices,]

# subset those genes based on if sig. diffntl expression
sig_diff_exp <- subset(temp_df3, temp_df3$significant.x == 'yes' & temp_df3$significant.y == 'yes')

# calculate the pearson correlation (linear) and the spearman correlation (monotonic, always increasing or decreasing)
res <- cor.test(sig_diff_exp$log2.fold_change..y ~ sig_diff_exp$log2.fold_change..x, method = 'pearson')
res2 <- cor.test(sig_diff_exp$log2.fold_change..y ~ sig_diff_exp$log2.fold_change..x, method = 'spearman')

## grab summary statistics of pearson correlation
