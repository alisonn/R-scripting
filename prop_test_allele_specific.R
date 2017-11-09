#!/bin/env R
# Test of proportions using allele specific enrichment/expression at 3 levels
# 1 between replicates 2 between treatments 3 genic level

# Load the replicates, treatments 
# a is summer, b is winter
args=commandArgs(TRUE)
a_rep1 <- read.table(args[1], stringsAsFactors=F)
a_rep2 <- read.table(args[2], stringsAsFactors=F)
b_rep1 <- read.table(args[3], stringsAsFactors=F)
b_rep2 <- read.table(args[4], stringsAsFactors=F)

a_merge <- merge(a_rep1, a_rep2, by=c("V1", "V2"))
b_merge <- merge(b_rep1, b_rep2, by=c("V1", "V2"))

cov_cutoff <- 10
a_merge_filter <- subset(a_merge, as.numeric(a_merge$V6.x) >= cov_cutoff & as.numeric(a_merge$V6.y) >= cov_cutoff )
b_merge_filter <- subset(b_merge, as.numeric(b_merge$V6.x) >= cov_cutoff & as.numeric(b_merge$V6.y) >= cov_cutoff )

### PHASE 1: PROP.TEST() per SNP SITE
# generating the appropriate data table for prop.test
# columns are chr, pos, alt_x, alt_y, coverage_x, coverage_y, alt_parent (1=summer), and the last is the pval from prop.test()
col_keep <- c(1,2,5,10,6,11,12)
a_proptest <- a_merge_filter[,col_keep]
b_proptest <- b_merge_filter[,col_keep]
a_proptest$pval <- 0; b_proptest$pval <- 0

# conduct proptest between replicates from same cross and get only the p-values
a_proptest$pval <- apply(a_proptest, 1, function(row) prop.test(x=as.numeric(c(row[3], row[4])), n=as.numeric(c(row[5], row[6])))$p.value)
b_proptest$pval <- apply(b_proptest, 1, function(row) prop.test(x=as.numeric(c(row[3], row[4])), n=as.numeric(c(row[5], row[6])))$p.value)

# calculate the average by merging the treatment allelic counts together - this is the same as weighted average
a_proptest$alt_freq <- (as.numeric(a_proptest$V5.x) + as.numeric(a_proptest$V5.y))/(as.numeric(a_proptest$V6.x) + as.numeric(a_proptest$V6.y))
b_proptest$alt_freq <- (as.numeric(b_proptest$V5.x) + as.numeric(b_proptest$V5.y))/(as.numeric(b_proptest$V6.x) +as.numeric(b_proptest$V6.y))

final_result <- merge(a_proptest, b_proptest, by=c("V1","V2"))

## OUT_FILE reports the consistency between replicates of same treatments in the cross
col_names = c("chr", "pos", "alt_S1", "alt_S2", "cov_S1", "cov_S2", "alt_parent", "prop_pval_S", "alt_freq_S", "alt_W1", "alt_W2", "cov_W1", "cov_W2", "alt_parent", "prop_pval_W", "alt_freq_W")
write.table(final_result, file=args[5], append=FALSE, quote = FALSE, sep = "\t", col.names = col_names, row.names=FALSE)
