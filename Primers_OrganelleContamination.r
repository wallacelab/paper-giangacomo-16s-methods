#! /usr/bin/Rscript

# Plot graphics determining how well each primer set prevented chloroplast contamination

# Libraries
library(argparse)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tidyr)
options(stringsAsFactors=F)

# Command-line arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Rdata file containing the phyloseq object to analyze (from step 2a)")
parser$add_argument("-o", "--outprefix", help="Prefix for all output files")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2020 03 Consolidated Pipeline/')
# args=parser$parse_args(c("-i", "TestPrimers/2_Analysis/2b_filtered_data.phyloseq.RDS", "-o", "99_tmp"))


# Load data
cat("Loading phyloseq data\n")
mydata = readRDS(args$infile) 

# Pull data back out of Phyloseq object so can plot exactly what we want
otus = otu_table(mydata)
taxonomy = as.data.frame(tax_table(mydata))
key = sample_data(mydata)

# Sanity checks to make sure things line up properly
if ( ! identical(rownames(otus), rownames(taxonomy))){
    stop("OTU names do not line up with taxonomy")
}
if ( ! identical(colnames(otus), rownames(key))){
    stop("Sample names do not line up with sample key")
}

# Calculate fraction of chloroplast and mitochondria contamination in each
cat("Calculating fraction organnelar DNA\n")
is_chloroplast = taxonomy$Order == 'Chloroplast'        # Chloroplasts are an order
is_mitochondria = taxonomy$Family == 'Mitochondria'     # Mitochondria are a family
chloro_counts = colSums(otus[is_chloroplast,])
mito_counts = colSums(otus[is_mitochondria,])
total_counts = colSums(otus)
other_counts = total_counts - chloro_counts - mito_counts

# Make data frame for plotting the above
tally = data.frame(sample=colnames(otus), targets=other_counts/total_counts, sample.type = key$sample.type, treatment = key$treatment) 
tally$organs = 1 - tally$targets

# Change levels of samples and treatments for publication
tally$treatment = sapply(as.character(tally$treatment), switch,
                                 BlockingOligos_v3v4="BO_3/4", 
                                 BlockingOligos_v5v7="BO_5/7",
                                 BlockingOligos_v5v7_noLinkers="BO_5/7",
                                 Discriminating = "Disc.",
                                 PNAs= "PNA", Universal="Univ.",
                                 NA) # NA catches anything that didn't match
tally$sample.type = sapply(as.character(tally$sample.type), switch,
                              "leaf-maize"="Maize Leaf", 
                              "leaf-soybean"="Soybean Leaf", 
                              "defined-community"="Defined Community", 
                              "soil-clay"="Soil 1",
                              "soil-flowerbed"="Soil 2",
                              NA) # NA catches anything that didn't match

# Split by sample type
plant_targets = droplevels(subset(tally, tally$sample.type %in% c("Maize Leaf", "Soybean Leaf")))
nonplant_targets = droplevels(subset(tally, tally$sample.type %in% c("Soil 1", "Soil 2", "Defined Community")))


# Plot actual graphic
plot_violins = function(mydata, mygroup=""){
  set.seed(1) # To keep jitters the same each time
  # Make plot
  myplot = ggplot(mydata, mapping=aes(x=treatment, y=organs)) + 
    geom_violin(scale="width", fill='lightblue', size=0.25) +                                          # Violin plots
    geom_jitter(mapping=aes(color=sample.type), alpha=0.5) +                                                       # Scatter points of actual data (jittered to be easier to see)
    labs(x="Primer Set", y="Fraction Organellar Reads") +
    theme_bw() + 
    theme(
      axis.title = element_text(size=11, face="bold"),
      axis.text = element_text(size=8, face="bold"), 
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), 
      legend.position = c(0.25, 0.82), # Legend in upper-left of plot
      legend.title = element_text(size=6, face='bold'),
      legend.text = element_text(size=5, face='bold'), 
      legend.margin = margin(3,3,3,3), 
      legend.box.background = element_rect(color='black'),
      legend.spacing.y = unit(0, "pt"),
      legend.key.size = unit(10, "pt") # adjusts the (hidden) background behind legend points; most important element for legend spacing
    ) +
    labs(color="Sample Type")
  # Save PNG and SVG
  ggsave(myplot, file= paste(args$outprefix, mygroup, "png", sep="."), width=2.5, height = 2.6, units="in", dpi=300)
  ggsave(myplot, file= paste(args$outprefix, mygroup, "svg", sep="."), width=2.5, height = 2.6, units="in", dpi=300)
}
plot_violins(plant_targets, "plant")
plot_violins(nonplant_targets, "nonplant")
