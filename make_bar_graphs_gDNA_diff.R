#!/bin/env R

### Makes Figure 3A bar graph showing differences in TE abundances, specifically LargeY-SmallY and MediumY-SmallY
## extra supplement figures are plotted below

library(ggplot2)
library(reshape2)
library(plyr)
library(viridis)

### inputs a bed file -- format: chr\tstart\tend\coverage
### next input is a size/length file -- format: chr\tlength
ssr2_file <- "DBCC035B6.betancourt2.sorted.bed"
kb10_file <- "DBCC035B7.betancourt2.sorted.bed"
dpse124_file <- "DBCC035B8.betancourt2.sorted.bed"
size_file <-"dpse-bar-graph.info.v2"
sizes <- read.table(size_file, header=F, stringsAsFactors=F)


## process files separately
process_file <- function(fileName,sizes) {
	data <- read.table(fileName, stringsAsFactors=F, header=F)	
	data <- join(data, sizes, by = 'V1', type = 'left', match = 'all')
	names(data) <- c('group', 'start', 'end', 'rawCov', 'fam', 'groupLen')
	data$wMean <- (data$end - data$start) / data$groupLen * data$rawCov ## gives coverage divided by proportion of the group
	agg <- aggregate(wMean ~ group + fam + groupLen, data, sum) ## sum all parts of the individual means

	## result is this which is good to go for subsetting and normalization
		## agg == c(group, groupLen, wMean ## the real wMean)
	muller <- c("Muller_A-AD", "Muller_B", "Muller_C", "Muller_E", "Muller_F", "Y_chrom", "Garbage")
	autosomes <- c("Muller_B", "Muller_C", "Muller_E", "Muller_F")
	normFactor <- mean((subset(agg, (agg$group %in% autosomes)))$wMean)
	print(normFactor)

	agg$wMean <- agg$wMean / normFactor * 2
		#print(subset(agg, agg$group %in% muller))
	result <- subset(agg, !(agg$group %in% muller))
	return(result)
}


ssr2 <- process_file(ssr2_file, sizes)
kb10 <- process_file(kb10_file, sizes)
dpse124 <- process_file(dpse124_file, sizes)

merge_data <- function(data1, data2, prefix) {
	result <- merge(data1, data2, by = "group")
	result$mean <- rowMeans(result[,c("wMean.x", "wMean.y")])
	result <- data.frame(result$group, result$fam.x, result$groupLen.x, result$wMean.x, result$wMean.y, result$mean)
	names(result) <- c("group", "fam", "length", paste(prefix, "1", sep = "_"), paste(prefix, "2", sep = "_"), paste(prefix, "mean", sep = "_"))
	return(result)
}

final <- merge(ssr2, kb10, by = "group")
final <- final[c("group", "fam.x", "groupLen.x", "wMean.x", "wMean.y")]
final <- merge(final, dpse124, by = "group")

final <- final[c("group", "fam", "groupLen", "wMean.x", "wMean.y", "wMean")]
names(final) <- c("TE", "Family", "Length", "SSR2_Count", "KB10_Count", "D124_Count")
final$ssr2_sum <- as.numeric(final$Length) * final$SSR2_Count 
final$kb10_sum <- as.numeric(final$Length) * final$KB10_Count
final$d124_sum <- as.numeric(final$Length) * final$D124_Count

final$kb10_plot <- final$kb10_sum - final$ssr2_sum
final$d124_plot <- final$d124_sum - final$ssr2_sum
toPlot <- data.frame(final$Family, final$kb10_plot, final$d124_plot)
names(toPlot) <- c("Family", "Medium_Y", "Large_Y")
toPlot <- aggregate(.~Family, toPlot, sum)
final.m <- melt(toPlot, id.vars="Family")

png("y_chromosome_diff.png", height = 400, width = 1000)
p <- ggplot(final.m, aes(reorder(Family, -value), value/(1e06))) + geom_bar(aes(fill=variable), position = "dodge", stat = "identity", colour='black') 
p <- p + xlab("Transposable element family") + ylab("Relative TE Content (Mb)") + theme_bw()
p <- p + theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size=15), axis.text.y=element_text(size=12.5, colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, size =9, colour = 'black'), legend.text=element_text(size=10)) 
p <- p + scale_fill_manual(values=c('#FDE725FF', '#287D8EFF'), name="Y-replacement Line")
print(p); dev.off()


## extra supplementary figures below 
## plotting absolute values

# val <- 25 ## example cutoff for values to plot
# toPlot <- data.frame(final$TE_Family, final$ssr2_sum, final$kb10_sum, final$d124_sum)
# toPlot <- data.frame(final$TE_Family, final$SSR2_Count, final$KB10_Count, final$D124_Count)
# names(toPlot) <- c("Family", "SSR2", "KB10", "D124")
# print(toPlot)
# toPlot <- aggregate(.~Family, toPlot, sum)
# toPlot <- subset(toPlot, toPlot$SSR2 > val & toPlot$KB10 > val & toPlot$D124 > val)

## toPlot.m <- melt(toPlot, id.vars="Family")
# p <- ggplot(toPlot.m , aes(reorder(Family, -value), value)) + geom_bar(aes(fill = variable), position = "dodge", stat = "identity") + xlab("TE family/group") + ylab("Estimated bp gain from TE") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size =15)) + scale_fill_manual(values = c('dodgerblue', 'dodgerblue3', 'dodgerblue4'))

#print(p)
#dev.off()
