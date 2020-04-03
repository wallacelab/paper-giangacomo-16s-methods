#! /usr/bin/Rscript

# Plot how well Deblur retained reads from samples

library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Stats CSV file exported from deblur step of QIIME2")
parser$add_argument("-o", "--outfile", help="Output graphic file")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/1_AssignOtus/BlockingOligos_3_4/')
# args=parser$parse_args(c("-i","deblur-stats.csv", "-o", "99_tmp.png"))

# Libraries
library(ggplot2)

# Load data
cat("Loading data\n")
stats = read.csv(args$infile)

# Reformat to be nice for ggplot
mapped=data.frame(sample = stats$sample.id, count = stats$reads.hit.reference, set="Mapped reads")
unmapped=data.frame(sample = stats$sample.id, count = stats$reads.raw - stats$reads.hit.reference, set="Unmapped reads")
mydata=rbind(mapped, unmapped)
mydata$set = factor(mydata$set, levels=c("Unmapped reads", "Mapped reads"))

# Make graphic
mygraphic = ggplot(data=mydata, mapping=aes(x=sample, y=count, fill=set)) +
    geom_bar(stat="identity", position="stack")
ggsave(args$outfile, plot=mygraphic)
