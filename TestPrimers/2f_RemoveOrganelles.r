#! /usr/bin/Rscript

library(argparse)
library(phyloseq)
library(ggplot2)
library(gridExtra)
options(stringsAsFactors=F)

# Command-line arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="RDS file containing the phyloseq object to analyze (from step 2a)")
parser$add_argument("-o", "--outprefix", help="Prefix for all output files")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 10 Mohsen Final Data/2_Analysis/')
# args=parser$parse_args(c("-i", "2b_filtered_data.phyloseq.RDS", "-o", "99_tmp"))

cat("Filtering out organelle data\n")

# Load phyloseq data
mydata = readRDS(args$infile)
taxonomy = as.data.frame(tax_table(mydata))
metadata = sample_data(mydata)

# Remove all organellar sequences
is_chloroplast = taxonomy$Order == 'Chloroplast'        # Chloroplasts are an order
is_mitochondria = taxonomy$Family == 'Mitochondria'     # Mitochondria are a family
to_remove = is_chloroplast | is_mitochondria
filtered = prune_taxa(mydata, taxa = !to_remove)

# Add columns to metadata for pre- and post-filtering
metadata$sample=rownames(metadata)
metadata$pre_filter=sample_sums(mydata)
metadata$post_filter=sample_sums(filtered)

# Make graphics of how read depth correlates before and after
cor_by_sample = ggplot(data=metadata, mapping=aes(x=pre_filter, y=post_filter, color=sample.type)) +
    geom_point(size=10, alpha=0.5) +
    xlab("Depth before filtering") + 
    ylab("Depth after filtering") + 
    ggtitle("Sample depth changes after removing organelles")
cor_by_treatment = ggplot(data=metadata, mapping=aes(x=pre_filter, y=post_filter, color=treatment)) +
    geom_point(size=10, alpha=0.5) +
    xlab("Depth before filtering") + 
    ylab("Depth after filtering") + 
    ggtitle("Sample depth changes after removing organelles")

# Export graphics
plots=list(cor_by_sample, cor_by_treatment)
png(paste(args$outprefix, ".changes.png", sep=""), width=15, height=5, units='in', res=300)
    grid.arrange(grobs=plots, ncol=length(plots))
dev.off()

# Write out filtered data
saveRDS(filtered, file=paste(args$outprefix, ".RDS", sep=""))
