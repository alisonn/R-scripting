#!/usr/bin/env R

# CALIBRATE GENOME SIZES IN R
# Takes in a data frame output by CALCULATE_GENOME_SIZES.R and
# first calculates the daily standard estimate and uses that estimate to calculate another estimate for the strain
# then uses the two estimates for the strain to calculate an average estimate
# output: pdf plots of each estimate based off of each standard

args = commandArgs(TRUE)
data_file = args[1]
out_pdf = args[2]
out_file = args[3]

my_data = read.table(data_file, stringsAsFactors = FALSE, header = TRUE)

# first assign new columns to fill in values
my_data$ere.mean = 0 
my_data$y.size2 = 0
my_data$y.size3 = 0

# for each of the days, calculate the mean of d.ere (if applicable)
#  assign those rows with given day the mean value for their d.ere column
print("We are looping through the erecta sizes to find the daily mean")
day_min = range(my_data$day)[1]
day_max = range(my_data$day)[2]

my_data$ere.mean = my_data$ere.size 

for (i in day_min : day_max ) {
	currDereSize = mean(my_data[grep( toString(i), my_data$ere.size ), ]$ere.size )
	my_data[grep(toString(i), my_data$day), ]$ere.mean = currDereSize
}

print("We are calculating the 2 other genome sizes")
# calculate y.size2 (y size based off of d. erecta) and y.size3 (mean of y.size and y.size2)
#my_data$y.size2 = my_data$peak2 / my_data$peak1 * my_data$ere.mean
#my_data$y.size3 = (my_data$y.size + my_data$y.size2) * 1/2

print("We are filtering the initial data set to keep only Y-lines")
# filter out Standards[1,2] and yw which were incorporated into D. ere daily average
my_data_pre1 = my_data[!grepl("1495", my_data$strain), ]
my_data_pre2 = my_data_pre1[!grepl("Standards1", my_data_pre1$strain), ]
filtered_data = my_data_pre2[!grepl("Standards2", my_data_pre2$strain), ]

print("We are merging the within day replicates and making a final data set")

# need to merge the daily replicates into one line to avoid double-counting for boxplots
# get the replicates (code: STRAIN + "-" + "[0-9]") and use aggregate to find the average
without_day_rep_data = filtered_data[!grepl("*-[1,2]", filtered_data$strain), ]
within_day_reps = filtered_data[grepl("*-[0-9]", filtered_data$strain), ]
within_day_reps$mainStrain = substr(within_day_reps$strain, 1, nchar(within_day_reps$strain)-2)
avg_within_day = aggregate(within_day_reps, list(day = within_day_reps$day, strain = within_day_reps$mainStrain), FUN = mean)
head(avg_within_day);
head(without_day_rep_data);
keep_col = c("strain", "day", "y.size", "y.size2", "y.size3", "ere.mean", "peak1", "peak2", "peak3")

# warning error here is OK, it cannot merge rows with characters - okay
# final_data has all yw and standards gone and all within-day replicates merged together
final_data = rbind(without_day_rep_data[keep_col], avg_within_day[keep_col])
print("We are calculating the 2 other genome sizes")
# calculate y.size2 (y size based off of d. erecta) and y.size3 (mean of y.size and y.size2)
final_data$y.size2 = final_data$peak2 / final_data$peak1 * final_data$ere.mean
final_data$y.size3 = (final_data$y.size + final_data$y.size2) * 1/2
final_data$type1 = with(final_data, reorder(final_data$strain, final_data$y.size, median))
final_data$type2 = with(final_data, reorder(final_data$strain, final_data$y.size2, median))
final_data$type3 = with(final_data, reorder(final_data$strain, final_data$y.size3, median))

# this is a series of plots demonstrating by-day and by-strain variation
pdf(out_pdf, width = 14, height = 24) 
par(mfrow = c(4,1))

plot(as.numeric(final_data$y.size) ~ final_data$day, main = "Daily Batch Effects in Size Estimates", ylab = "Genome Size Estimate (Mb)", xlab = "Batch by Day", las = 2)
boxplot(final_data$y.size ~ final_data$type1, data = final_data, main = "Genome Size by D. virilis Standard", xlab = "Y-replacement Line", ylab = "Genome Size Estimate (Mb)", las = 2, col = 'lightpink')
boxplot(final_data$y.size2 ~ final_data$type2, data = final_data, main = "Genome Size by D. erecta Standard", xlab = "Y-replacement Line", ylab = "Genome Size Estimate (Mb)", las = 2, col = 'lightgreen')
boxplot(final_data$y.size3 ~ final_data$type3, data = final_data, main = "Genome Size Averaging", xlab = "Y-replacement Line", ylab = "Genome Size Estimate (Mb)", las = 2, col = 'lightblue')

dev.off()
# write to a separate output file (Y-line_all_estimates.tsv)
write.table(final_data, file = out_file, append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
