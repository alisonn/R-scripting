#!/usr/bin/env R

# CALCULATE GENOME SIZES IN R 
# Uses a kernel density approximation to identify the peak fluorescence of cell populations containing different amounts of DNA. 
# @param: a file with gate channels repeated as often as they are observed in .fsc and separated by new lines
# @param: the script tsv_to_hist.py will generate a compatible input file
# @return: a master-file with the top peaks of fluorescence for the different cell populations

args = commandArgs(TRUE)
my_file = args[1]
out_file = args[2]
out_pdf = args[3]

# phase 1: read in the histogram data
my_data = read.table(my_file, stringsAsFactors = FALSE, header = FALSE)

# phase 2: generate the kernel density 
# default model is gaussian and bw = bandwidth (no difference b/w (400,450) BUT check in case)
my_density = density(my_data$V1, bw = 400)
maxima = my_density$x[which(diff(sign(diff(my_density$y))) == -2)]

# filter the maxima based on empirical cutoffs for gate channels 
prefilterMaxima = c(maxima[maxima <= 19000], maxima[maxima >= 28000])
filterMaxima = c(prefilterMaxima[prefilterMaxima <= 33000])
subMaxima = c(filterMaxima[filterMaxima >= 12500])

# phase 3: generate the plot
pdf(out_pdf, width = 6, height = 4)

# obtain the name of the y-line strain and the day batch
file_split = strsplit(my_file, "_")
strain = toString(strsplit(file_split[[1]][2], ".tsv")[[1]][1])
batch = toString(strsplit(file_split[[1]][1], "y-lines")[[1]][2])

hist(my_data$V1, breaks = 400, prob = TRUE, main = paste(strain,"day", batch, sep = "_"), xlab = "Gate Channel", ylab = "Density", cex.lab = 0.75, cex.main = 0.75, cex.axis = 0.75)
lines(my_density, col = 'lightblue', lwd = 2)
abline(v = subMaxima, col = 'red', lty = 2, lwd = 2)

dev.off()

# phase 4: calculate genome size 
# [1] from D. vir [2] from D. ere
erecta = subMaxima[1]/subMaxima[length(subMaxima)]*328
estVir = subMaxima[2]/subMaxima[3]*328

# phase 5: write the output file

toWrite = paste(strain, batch, subMaxima[1], subMaxima[2], subMaxima[3], estVir, erecta, sep = "\t")
# y-strain \t batch(day) \t maxima \t estimates
# just have a master sheet for all of the y-lines - to make it easier for anova testing
write(toWrite, out_file, append = TRUE, sep = "\t")
# close
