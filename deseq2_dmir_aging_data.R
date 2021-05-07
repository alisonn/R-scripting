### cleaning up featurecounts output from newgtf with exon information
data <- read.table("D.mir_aging_RNA_tissues_featureCounts_genes.new.gff.counts", stringsAsFactors = F, header=T)
data$Chr <- lapply(strsplit(data$Chr,';'), function(l) l[[1]])
data$Start <- lapply(strsplit(data$Start, ';'), function(l) l[[1]])
data$End <- sapply(strsplit(data$End, ';'), tail, 1)
data$Strand <- lapply(strsplit(data$Strand, ';'), function(l) l[[1]])
data2 <- data.frame(lapply(data, as.character), stringsAsFactors = FALSE)
write.table(data2, file="D.mir_aging_RNA_tissues_featureCounts_genes.new.gff.reformat.counts", quote=FALSE, col.names=TRUE, row.names=FALSE, sep = "\t")

### Load libraries ###

library(gplots)
library(DESeq2)
library(EnhancedVolcano)

data <- read.table("D.mir_aging_RNA_samples_featureCounts_genes.new.gff.reformat.counts", header=T, stringsAsFactors = F ,row= 1)
data1 <- read.table("D.mir_aging_RNA_samples_featureCounts_betancourtTEs.counts", header=F, stringsAsFactors = F, row=1)
data1 <- read.table("D.mir_aging_RNA_samples_featureCounts_mirKevinRepeats.counts", header=F, stringsAsFactors = F, row=1)
colnames(data1) <- colnames(data)
data <- rbind(data, data1)
#data <- read.table("D.mir_aging_RNA_samples_featureCounts_genesAndTEs.counts", header=T, stringsAsFactors=F ,row= 1)

######################### MALE DATA ##################################
maleData <- (data[,12:17])
meta <- read.table("D.mir_aging_RNA_samples_metadata.tsv",stringsAsFactors = F, header = T)
colDataMale <- meta[meta$sex == 'Male', ]

dds_male <- DESeqDataSetFromMatrix(countData=maleData, colData=colDataMale, design = ~ age)

colData(dds_male)$age <- factor(colData(dds_male)$age, levels=c("Old", "Young"))  ## log2(fc) is log2(old/young)
dds_male <- DESeq(dds_male)
res_m <- results(dds_male)

normMaleCounts <- as.data.frame(counts(dds_male, normalized=TRUE))
maleTECounts <- normMaleCounts[!grepl("gene", row.names(normMaleCounts)),]
maleTECounts <- maleTECounts+1
row.names(maleTECounts) <- lapply(strsplit(row.names(maleTECounts), '.',fixed = TRUE), function(l) l[[1]][1])

## check replicate data in individual plots
plot(log(maleTECounts$MO1,base=2) ~ log(maleTECounts$MY1,base=2), pch=20); abline(0,1,col='red')
with(maleTECounts, text(log(MO1,base=2) ~ log(MY1,base=2), labels=row.names(maleTECounts), pos=4, cex=.5))
plot(log(maleTECounts$MO2,base=2) ~ log(maleTECounts$MY2,base=2),pch=20); abline(0,1,col='red')
with(maleTECounts, text(log(MO2,base=2) ~ log(MY2,base=2), labels=row.names(maleTECounts), pos=4, cex=1))
plot(log(maleTECounts$MO3,base=2) ~ log(maleTECounts$MY3,base=2),pch=20); abline(0,1,col='red')
with(maleTECounts, text(log(MO3,base=2) ~ log(MY3,base=2), labels=row.names(maleTECounts), pos=3, cex=1))

## check data again
l <- c(0,20)
plot(log(normMaleCounts$MO1,base=2) ~ log(normMaleCounts$MY1,base=2), ylim=l, xlim=l); abline(0,1,col='red')
plot(log(normMaleCounts$MO2,base=2) ~ log(normMaleCounts$MY2,base=2), ylim=l, xlim=l); abline(0,1,col='red')
plot(log(normMaleCounts$MO3,base=2) ~ log(normMaleCounts$MY3,base=2), ylim=l, xlim=l); abline(0,1,col='red')


### plotting volcano plot for Repeats Only in males
res_m_df <- as.data.frame(res_m)
res_m_df_repeats <- res_m_df[!grepl("gene", row.names(res_m_df)),]

EnhancedVolcano(res_m_df_repeats,
                lab = rownames(res_m_df_repeats),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8, 8),
                pCutoff = 0.05,
                FCcutoff = .58,
                pointSize = 3.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

### volcano Plot of genes only - male samples
res_m_df_genes <- res_m_df[grepl("gene", row.names(res_m_df)),]
EnhancedVolcano(res_m_df_genes,
                lab = rownames(res_m_df_genes),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8, 8),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

################ FEMALE DATA ######################
femaleData <- data[,6:11]
meta <- read.table("D.mir_aging_RNA_samples_metadata.tsv",stringsAsFactors = F, header = T)
colDataFemale <- meta[meta$sex == 'Female', ]

dds_female <- DESeqDataSetFromMatrix(countData=femaleData, colData=colDataFemale, design = ~ age)
colData(dds_female)$sex <- factor(colData(dds_female)$age, levels=c("Old", "Young")) 
dds_female <- DESeq(dds_female)
res_f <- results(dds_female)
normFemaleCounts <- as.data.frame(counts(dds_female, normalized=TRUE))
femaleTECounts <- normFemaleCounts[!grepl("gene", row.names(normFemaleCounts)),]
femaleTECounts <- femaleTECounts+1
row.names(femaleTECounts) <- lapply(strsplit(row.names(femaleTECounts), '.',fixed = TRUE), function(l) l[[1]][1])

### check female data
plot(log(femaleTECounts$FO1,base=2) ~ log(femaleTECounts$FY1,base=2), pch=20); abline(0,1,col='red')
with(femaleTECounts, text(log(FO1,base=2) ~ log(FY1,base=2), labels=row.names(femaleTECounts), pos=4, cex=1))
plot(log(femaleTECounts$FO2,base=2) ~ log(femaleTECounts$FY2,base=2),pch=20); abline(0,1,col='red')
with(femaleTECounts, text(log(FO2,base=2) ~ log(FY2,base=2), labels=row.names(femaleTECounts), pos=4, cex=1))
plot(log(femaleTECounts$FO3,base=2) ~ log(femaleTECounts$FY3,base=2),pch=20); abline(0,1,col='red')
with(femaleTECounts, text(log(FO3,base=2) ~ log(FY3,base=2), labels=row.names(femaleTECounts), pos=3, cex=1))

## check female data again
l <- c(0,20)
plot(log(normFemaleCounts$FO1,base=2) ~ log(normMaleCounts$FY1,base=2), ylim=l, xlim=l); abline(0,1,col='red')
plot(log(normFemaleCounts$FO2,base=2) ~ log(normFemaleCounts$FY2,base=2), ylim=l, xlim=l); abline(0,1,col='red')
plot(log(normFemaleCounts$FO3,base=2) ~ log(normFemaleCounts$FY3,base=2), ylim=l, xlim=l); abline(0,1,col='red')


### plotting volcano plot for female repeat expression
res_f_df <- as.data.frame(res_f)
res_f_df_repeats <- res_f_df[!grepl("gene", row.names(res_f_df)),]

EnhancedVolcano(res_f_df_repeats,
                lab = rownames(res_f_df_repeats),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8, 8),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


### Volcano plot of genes in females
res_f_df_genes <- res_f_df[grepl("gene", row.names(res_f_df)),]
EnhancedVolcano(res_f_df_genes,
                lab = rownames(res_f_df_genes),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8, 8),
                ylim=c(0,70),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


#### Boxplot REPEATS
maleTECounts$MOMean <- rowMeans(subset(maleTECounts, select = c("MO1", "MO2", "MO3"), na.rm = TRUE))
maleTECounts$MYMean <- rowMeans(subset(maleTECounts, select = c("MY1", "MY2", "MY3"), na.rm = TRUE))
femaleTECounts$FOMean <- rowMeans(subset(femaleTECounts, select = c("FO1", "FO2", "FO3"), na.rm = TRUE))
femaleTECounts$FYMean <- rowMeans(subset(femaleTECounts, select = c("FY1", "FY2", "FY3"), na.rm = TRUE))
femaleTECounts$FFC <- femaleTECounts$FOMean/femaleTECounts$FYMean 
maleTECounts$MFC <- maleTECounts$MOMean/maleTECounts$MYMean 
mergedTECounts <- cbind(maleTECounts, femaleTECounts)[c("FYMean", "FOMean", "MYMean", "MOMean", "FFC", "MFC")]

### pick which boxplot to graph: samples or fold-change
mergedTECounts.s <- subset(mergedTECounts, mergedTECounts$FYMean > 1 & mergedTECounts$FOMean > 1 & mergedTECounts$MYMean > 1 &  mergedTECounts$MOMean > 1)[,c("FYMean", "FOMean", "MYMean", "MOMean")]
mergedTECounts.s <- subset(mergedTECounts, mergedTECounts$FYMean > 0 | mergedTECounts$FOMean > 0 | mergedTECounts$MYMean > 0 |  mergedTECounts$MOMean > 0)[,c("FYMean", "FOMean", "MYMean", "MOMean")]


mergedTECounts.s <- subset(mergedTECounts, mergedTECounts$FYMean > 0 & mergedTECounts$FOMean > 0 & mergedTECounts$MYMean > 0 &  mergedTECounts$MOMean > 0)[,c("FFC", "MFC")]
mergedTECounts.s <- subset(mergedTECounts, mergedTECounts$FYMean > 0 | mergedTECounts$FOMean > 0 | mergedTECounts$MYMean > 0 |  mergedTECounts$MOMean > 0)[,c("FFC", "MFC")]

mergedTECounts.s.m <- melt(mergedTECounts.s)
## pick which ggplot for the graph - samples or foldchange
ggplot(mergedTECounts.s.m, aes(x=variable, y=log(value, base=2), fill=variable)) + geom_boxplot(notch=TRUE) + theme_bw()+ coord_cartesian(ylim=c(0,20)) + scale_fill_manual(values=c("red","pink", "blue", "cyan"))
ggplot(mergedTECounts.s.m, aes(x=variable, y=log(value, base=2), fill=variable)) + geom_boxplot(notch=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-2,2)) + scale_fill_manual(values=c("red", "blue"))



## Heatmap REPEATS
library(RColorBrewer)
myColors <- colorRampPalette(c("blue","white", "red"))(n=100)

### Heatmap repeats by log2(female/male) foldchange
ylog2 <- read.table("/Users/alisonn/Desktop/data/dmir_repeats/TE.DNA.FC.txt", header=F, stringsAsFactors = F)
ylog2 <- ylog2[order(ylog2["V2"], decreasing=T),]
ylog2$V1 <- lapply(strsplit(ylog2$V1, ':'), function(l) l[[1]])
## sort mergedTECounts by ylog2 
test <- mergedTECounts[order(match(row.names(mergedTECounts),ylog2[,1])),]
test <- log(test, base=2)
test <- subset(test, test$FYMean > 0 & test$FOMean > 0 & test$MYMean > 0 &  test$MOMean > 0)
png(filename = "test2.png", height=1000, width=400, units="px")
heatmap.2(as.matrix(test[,c(1:4)]), trace='none', Rowv = F, Colv = F, srtCol = 30, col=myColors,
          scale='row', lmat=rbind(c(4, 2), c(1, 3)), lhei=c(2, 8), lwid=c(4, 1), density.info='none')
dev.off()


### plotting scatterfplots with error bars by age (so separate plots for each sex)
mergedTECounts <- cbind(maleTECounts, femaleTECounts) ### this keeps rownames :)
#mergedTECounts.s <- subset(mergedTECounts, mergedTECounts$FYMean > 0 & mergedTECounts$FOMean > 0 & mergedTECounts$MYMean > 0 &  mergedTECounts$MOMean > 0)
mergedTECounts.s <- mergedTECounts 

## convert all -Infs (initially 0 expression) to 0's -- tell Doris this
# mergedTECounts.s.log2 <- (log(mergedTECounts.s, base=2))
# 
# ### TRYING THIS OUT 3/26/2020
# is.na(mergedTECounts.s.log2)<-sapply(mergedTECounts.s.log2, is.infinite)
# mergedTECounts.s.log2[is.na(mergedTECounts.s.log2)]<-0

## calculate ranges for each sample type 
mergedTECounts.s.log2$MO_SD <- apply(mergedTECounts.s.log2[,c("MO1", "MO2", "MO3")], 1,sd )
mergedTECounts.s.log2$MY_SD <- apply(mergedTECounts.s.log2[,c("MY1", "MY2", "MY3")], 1,sd )
mergedTECounts.s.log2$FO_SD <- apply(mergedTECounts.s.log2[,c("FO1", "FO2", "FO3")], 1,sd )
mergedTECounts.s.log2$FY_SD <- apply(mergedTECounts.s.log2[,c("FY1", "FY2", "FY3")], 1,sd )
row.names(mergedTECounts.s.log2) <- row.names(maleTECounts)

## identifying significantly differentially expressed 
#################~~~~~~~~~~~~~~~~~~~~~~~~####################3
############ write the EXPRESSION RESULTS files!!!!#######################
write.table(mergedTECounts.s.log2, "MSH22_aging_Satellites_DESeq2_expression_v2.tsv", quote = F, sep='\t', row.names = T, col.names = T)


row.names(res_m_df_repeats) <- lapply(strsplit(row.names(res_m_df_repeats), '.', fixed=TRUE), function(l) l[[1]])
yenrich <- read.table("../../dmir_repeats/yenriched_TEs.txt", header=T, row=1)
yenriched <- subset(yenrich, V2 == "1")$V1
res_m_df_repeats$c <- 0 ## reset the column
res_m_df_repeats$c[row.names(res_m_df_repeats) %in% yenriched] <- 1
res_m_df_repeats$c[abs(res_m_df_repeats$log2FoldChange) >= 0.58 & res_m_df_repeats$pvalue < 0.05 ] <- 2
res_m_df_repeats$c[row.names(res_m_df_repeats) %in% yenriched & abs(res_m_df_repeats$log2FoldChange) >= 0.58 & res_m_df_repeats$pvalue < 0.05 ] <- 3
####
mergedTECounts.s.log2$c <- res_m_df_repeats$c



####### BETWEEN SEX ANALYSIS ### 

maleData <- (data[,12:17])
meta <- read.table("D.mir_aging_RNA_samples_metadata.tsv",stringsAsFactors = F, header = T)
colDataMale <- meta[meta$sex == 'Male', ]

colDataYoung <- meta[meta$age == 'Young', ]
youngData <- data[,c("FY1", "FY2", "FY3", "MY1", "MY2", "MY3")]

dds_young <- DESeqDataSetFromMatrix(countData=youngData, colData=colDataYoung, design = ~ sex)
## just let deseq2 do its default
colData(dds_young)$age <- factor(colData(dds_young)$age, levels=c("Female", "Male"))  ## log2(fc) is log2(old/young)
dds_young <- DESeq(dds_young)
res_y <- results(dds_young)
normYoungCounts <- as.data.frame(counts(dds_young, normalized=TRUE))
youngTECounts <- normYoungCounts[!grepl("gene", row.names(normYoungCounts)),]
youngTECounts <- youngTECounts+1


youngData <- data[,c("FY1", "FY2", "FY3", "MY1", "MY2", "MY3")]

colDataYoung <- meta[meta$age == 'Young', ]

dds_young <- DESeqDataSetFromMatrix(countData=youngData, colData=colDataYoung, design = ~ sex)
colData(dds_young)$sex<- factor(colData(dds_young)$sex, levels=c("Female", "Male"))  ## log2(fc) is log2(Male/Female)
dds_young <- DESeq(dds_young)
res_y <- results(dds_young)
normYoungCounts <- as.data.frame(counts(dds_young, normalized=TRUE))
youngTECounts <- normYoungCounts[!grepl("gene", row.names(normYoungCounts)),]
youngTECounts <- youngTECounts+1


res_y_df <- as.data.frame(res_y)
res_y_df_repeats <- res_y_df[!grepl("gene", row.names(res_y_df)),]

### coloring significantly and differentially expressed
yenriched <- subset(yenrich, V2 == "1")$V1
row.names(res_y_df_repeats) <- lapply(strsplit(row.names(res_y_df_repeats), '.', fixed=TRUE), function(l) l[[1]])

res_y_df_repeats$c <- 0 ## reset the column
res_y_df_repeats$c[row.names(res_y_df_repeats) %in% yenriched] <- 1
res_y_df_repeats$c[abs(res_y_df_repeats$log2FoldChange) >= 0.58 & res_y_df_repeats$pvalue < 0.05 ] <- 2
res_y_df_repeats$c[row.names(res_y_df_repeats) %in% yenriched & abs(res_y_df_repeats$log2FoldChange) >= 0.58 & res_y_df_repeats$pvalue < 0.05 ] <- 3
mergedTECounts.s.log2$c <- res_y_df_repeats$c

p <- ggplot(mergedTECounts.s.log2, aes(x=FYMean, y=MYMean)) + geom_point(aes(color=factor(c)))+ theme_bw() + geom_abline(slope=1,intercept=0,col='grey')
p <- p + geom_errorbar(aes(ymax = MYMean+MY_SD, ymin=MYMean-MY_SD,alpha=0.5), color=ifelse(mergedTECounts.s.log2$c == '3', "dodgerblue2", ifelse(mergedTECounts.s.log2$c == '2', "orange", ifelse(mergedTECounts.s.log2$c == '1', "cyan3", "black"))))
p <- p + geom_errorbarh(aes(xmax=FYMean+FY_SD, xmin=FYMean-FY_SD, alpha=0.5),color=ifelse(mergedTECounts.s.log2$c == '3', "dodgerblue2", ifelse(mergedTECounts.s.log2$c == '2', "orange", ifelse(mergedTECounts.s.log2$c == '1', "cyan3", "black"))))
p <- p + coord_cartesian(ylim=c(-1,20), xlim=c(-1,20))
p <- p + scale_color_manual(values=c("black", "cyan3", "orange","dodgerblue2")) + ylab("Log2(Young Male Expression)") + xlab("Log2(Young Female Expression)")
print(p)

####
row.names(res_y_df_repeats) <- lapply(strsplit(row.names(res_y_df_repeats), '.', fixed=TRUE), function(l) l[[1]])
res_y_df_repeats$c <- ifelse(abs(res_y_df_repeats$log2FoldChange) >= 0.58 & res_y_df_repeats$pvalue < 0.05, 1, 0)
res_y_df_repeats$c <- ifelse(abs(res_y_df_repeats$log2FoldChange) >= 0.58 & res_y_df_repeats$pvalue < 0.05 &row.names(res_y_df_repeats) %in% yenriched , 1, 0)

res_y_df_repeats$c[is.na(res_y_df_repeats$c)] <- 0
mergedTECounts.s.log2$c <- res_y_df_repeats$c

p <- ggplot(mergedTECounts.s.log2, aes(x=FYMean, y=MYMean)) + geom_point(aes(color=factor(c)))+ theme_bw() + geom_abline(slope=1,intercept=0,col='grey')
p <- p + geom_errorbar(aes(ymax = MYMean+MY_SD, ymin=MYMean-MY_SD,alpha=0.5), color=ifelse(mergedTECounts.s.log2$c == '1', "orange", "black"))
p <- p + geom_errorbarh(aes(xmax=FYMean+FY_SD, xmin=FYMean-FY_SD, alpha=0.5), color=ifelse(mergedTECounts.s.log2$c == '1', "orange", "black"))
p <- p + coord_cartesian(ylim=c(-1,20), xlim=c(-1,20))
p <- p + scale_color_manual(values=c("black", "orange")) + ylab("Log2(Young Male Expression)") + xlab("Log2(Young Female Expression)")
print(p)



#######
### Revised Old Data

dds_old <- DESeqDataSetFromMatrix(countData=oldData, colData=colDataold, design = ~ sex)
## just let deseq2 do its default
colData(dds_old)$age <- factor(colData(dds_old)$age, levels=c("Female", "Male"))  ## log2(fc) is log2(old/young)
dds_old <- DESeq(dds_old)
res_o <- results(dds_old)
normOldCounts <- as.data.frame(counts(dds_old, normalized=TRUE))
oldTECounts <- normOldCounts[!grepl("gene", row.names(normOldCounts)),]
oldTECounts <- oldTECounts+1


###

oldData <- data[,c("FO1", "FO2", "FO3", "MO1", "MO2", "MO3")]
colDataold <- meta[meta$age == 'Old', ]

dds_old <- DESeqDataSetFromMatrix(countData=oldData, colData=colDataold, design = ~ sex)
colData(dds_old)$sex<- factor(colData(dds_old)$sex, levels=c("Female", "Male"))  ## log2(fc) is log2(Male/Female)
dds_old <- DESeq(dds_old)
res_o <- results(dds_old)
normoldCounts <- as.data.frame(counts(dds_old, normalized=TRUE))
oldTECounts <- normoldCounts[!grepl("gene", row.names(normoldCounts)),]
oldTECounts <- oldTECounts+1


res_o_df <- as.data.frame(res_o)
res_o_df_repeats <- res_o_df[!grepl("gene", row.names(res_o_df)),]


row.names(res_o_df_repeats) <- lapply(strsplit(row.names(res_o_df_repeats), '.', fixed=TRUE), function(l) l[[1]])
#### Recolor with 4 levels
yenriched <- subset(yenrich, V2 == "1")$V1
res_o_df_repeats$c <- 0 ## reset the column
res_o_df_repeats$c[row.names(res_o_df_repeats) %in% yenriched] <- 1
res_o_df_repeats$c[abs(res_o_df_repeats$log2FoldChange) >= 0.58 & res_o_df_repeats$pvalue < 0.05 ] <- 2
res_o_df_repeats$c[row.names(res_o_df_repeats) %in% yenriched & abs(res_o_df_repeats$log2FoldChange) >= 0.58 & res_o_df_repeats$pvalue < 0.05 ] <- 3
####
mergedTECounts.s.log2$c <- res_o_df_repeats$c

p <- ggplot(mergedTECounts.s.log2, aes(x=FOMean, y=MOMean)) + geom_point(aes(color=factor(c)))+ theme_bw() + geom_abline(slope=1,intercept=0,col='grey')
p <- p + geom_errorbar(aes(ymax = MOMean+MO_SD, ymin=MOMean-MO_SD,alpha=0.5),color=ifelse(mergedTECounts.s.log2$c == '3', "dodgerblue2", ifelse(mergedTECounts.s.log2$c == '2', "orange", ifelse(mergedTECounts.s.log2$c == '1', "cyan3", "black"))))
p <- p + geom_errorbarh(aes(xmax=FOMean+FO_SD, xmin=FOMean-FO_SD, alpha=0.5),color=ifelse(mergedTECounts.s.log2$c == '3', "dodgerblue2", ifelse(mergedTECounts.s.log2$c == '2', "orange", ifelse(mergedTECounts.s.log2$c == '1', "cyan3", "black"))))
p <- p + coord_cartesian(ylim=c(-1,20), xlim=c(-1,20))
p <- p + scale_color_manual(values=c("black", "cyan3", "orange","dodgerblue2")) + ylab("Log2(Old Male Expression)") + xlab("Log2(Old Female Expression)")
print(p)


## 2 levels 
res_o_df_repeats$c <- ifelse(abs(res_o_df_repeats$log2FoldChange) >= 0.58 & res_o_df_repeats$pvalue < 0.05, 1, 0)
res_o_df_repeats$c[is.na(res_o_df_repeats$c)] <- 0
mergedTECounts.s.log2$c <- res_o_df_repeats$c


p <- ggplot(mergedTECounts.s.log2, aes(x=FOMean, y=MOMean)) + geom_point(aes(color=factor(c)))+ theme_bw() + geom_abline(slope=1,intercept=0,col='grey')
p <- p + geom_errorbar(aes(ymax = MOMean+MO_SD, ymin=MOMean-MO_SD,alpha=0.5), color=ifelse(mergedTECounts.s.log2$c == '1', "orange", "black"))
p <- p + geom_errorbarh(aes(xmax=FOMean+FO_SD, xmin=FOMean-FO_SD, alpha=0.5), color=ifelse(mergedTECounts.s.log2$c == '1', "orange", "black"))
p <- p + coord_cartesian(ylim=c(-1,20), xlim=c(-1,20))
p <- p + scale_color_manual(values=c("black", "orange")) #+ ylab("Log2(Old Male Expression)") + xlab("Log2(Old Female Expression)")
print(p)

