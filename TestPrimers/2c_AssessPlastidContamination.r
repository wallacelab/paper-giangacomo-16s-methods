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
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/2_Analysis/')
# args=parser$parse_args(c("-i", "2b_filtered_data.phyloseq.RDS", "-o", "2c_tmp"))



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

# Make a data frame for plotting the above
contamination = data.frame(sample=colnames(otus), chloro = chloro_counts, mito = mito_counts, other=other_counts, sample.type = key$sample.type, treatment = key$treatment)

# Tidy it up for plotting
counts = gather(contamination, source, amount, -sample, -sample.type, -treatment)
counts = arrange(counts, treatment, sample)

# Filter out the non-plant samples
plantdata = subset(counts, counts$sample.type %in% c("leaf-maize", "leaf-soybean"))
plantdata = droplevels(plantdata)

# Plot column plot of chloroplast, mitochondria, and other sequences
cat("Creating plots with output prefix",args$outprefix,"\n")
ggplot(plantdata, mapping=aes(x=sample, y=amount, fill=source)) + 
    geom_col(position="fill") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +   # Rotate text 90 degrees
    scale_y_continuous(expand = c(0,0)) +                   # Remove gray space above/below graphs
    facet_grid(. ~ sample.type, scales="free_x")  +         # Two panels, one for each sample type 
    scale_fill_manual("legend", values=c("chloro" = 'green', 'mito' = 'red', 'other' = 'blue'))
ggsave(paste(args$outprefix, ".column_plot.png", sep=""))

# Plot violin plot of fraction of non-organelle reads in each sample
got_targets = data.frame(sample=colnames(otus), targets=other_counts/total_counts, sample.type = key$sample.type, treatment = key$treatment) 
got_targets$organs = 1 - got_targets$targets
plant_targets = subset(got_targets, got_targets$sample.type %in% c("leaf-maize", "leaf-soybean"))
ggplot(plant_targets, mapping=aes(x=treatment, y=organs)) + 
  geom_violin(scale="width", fill='brown1') +                                          # Violin plots
  geom_jitter() +                                                       # Scatter points of actual data (jittered to be easier to see)
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))  +   # Rotate text 90 degrees
  labs(x="Primer Set", y="Fraction Organellar Reads", title="Fraction Organnellar Reads in Plant Samples")
ggsave(paste(args$outprefix, ".violin_plot.png", sep=""), width=78, height = 78, units="mm")

# Same plot, but for non-plant targets (=sanity check)
other_targets = subset(got_targets, got_targets$sample.type %in% c("defined-community", "soil-clay", "soil-flowerbed"))
ggplot(other_targets, mapping=aes(x=treatment, y=organs)) + 
  geom_violin(scale="width", fill='brown1') +                                          # Violin plots
  geom_jitter() +                                                       # Scatter points of actual data (jittered to be easier to see)
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))  +   # Rotate text 90 degrees
  labs(x="Primer Set", y="Fraction Organellar Reads", title="Fraction Organnellar Reads in Non-Plant Samples")
ggsave(paste(args$outprefix, ".violin_plot.non_plant_samples.png", sep=""))

# Write out plant organelle data table for raw numbers, including a quick summary
plant_targets_summary = split(plant_targets, droplevels(plant_targets$treatment))
plant_targets_summary = lapply(plant_targets_summary, function(t){
  data.frame(treatment=unique(t$treatment), min_organelle=min(t$organs), 
             max_organelle=max(t$organs))
})
plant_targets_summary = do.call(rbind, plant_targets_summary)
write.csv(plant_targets, file=paste(args$outprefix, ".organelle_fraction.csv", sep=""))
write.csv(plant_targets_summary, file=paste(args$outprefix, ".organelle_fraction.summary.csv", sep=""))



# ##########
# Violin plot for publication
# ##########
pretty_plants = plant_targets

# Change levels of treatment for publication
pretty_plants$treatment = sapply(as.character(pretty_plants$treatment), switch,
                                 BlockingOligos_v3v4="BO_3/4", 
                                 BlockingOligos_v5v7="BO_5/7",
                                 BlockingOligos_v5v7_noLinkers="BO_5/7",
                                 Discriminating = "Disc.",
                                 PNAs= "PNA", Universal="Univ.",
                                 NA) # NA catches anything that didn't match
# Plot actual graphic
set.seed(1) # To keep jitters the same each time
ggplot(pretty_plants, mapping=aes(x=treatment, y=organs)) + 
  geom_violin(scale="width", fill='lightblue', size=0.25) +                                          # Violin plots
  geom_jitter() +                                                       # Scatter points of actual data (jittered to be easier to see)
  labs(x="Primer Set", y="Fraction Organellar Reads") +
  theme_bw() + 
  theme(
    axis.title = element_text(size=11, face="bold"),
    axis.text = element_text(size=8, face="bold")
  )
ggsave(paste(args$outprefix, ".violin_plot.publication.png", sep=""), 
       width=78, height = 78, units="mm", dpi=300)

