#! /usr/bin/Rscript

# Analyze the effect of different methods on alpha diversity

library(phyloseq)
library(argparse)
library(ggplot2)
library(gridExtra)
options(stringsAsFactors=F)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="RDS file with saved phyloseq object for analysis")
parser$add_argument("-o", "--outprefix", help="Output file prefix for graphics")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2020 03 Consolidated Pipeline/')
# args=parser$parse_args(c("-i","TestExtraction/2_Analysis/2f_otu_table.no_organelles.RDS", "-o",'99_tmp'))


# Load data
cat("Loading data for alpha diversity analysis\n")
source("StandardizeLabels.r")
mydata = standardize_labels(readRDS(args$infile), type='extraction')


# ############
# Use facet_grid() to construct figure; since plants and soils are so different in scale, plots created separately and then manually joined together
# ############

# Custom labels for the facets
metrics = c("Observed OTUs", "Shannon Diversity")
names(metrics) = c("Observed", "Shannon")

# Split off data for plants and soils (set data to missing rather than subsetting so that image components stay the same size)
# plantdata = soildata = mydata
is_soil = sample_data(mydata)$sample.type=="Soil"
is_plant = !is_soil
plantdata = prune_samples(is_plant, mydata)
soildata = prune_samples(is_soil, mydata)
#otu_table(plantdata)[,!is_plant] = 0
#otu_table(soildata)[,!is_soil] = 0

# Make plots
richness_plot = function(plotdata, tag="", ...){
    myplot = plot_richness(plotdata, measures=c("Observed", "Shannon"), x="treatment", color="treatment") +
        geom_point(size=2.5, stroke=1) +
#         facet_wrap(variable ~ sample.type, scales="free", labeller=labeller(variable = metrics), nrow=2) +
        facet_grid(variable ~ sample.type, scales="free", labeller=labeller(variable = metrics)) +
        theme_bw() + 
        theme(text=element_text(face="bold"), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        xlab("Samples") +
        ylab("Alpha Diversity") +
        labs(color="Extraction Method")
    
    # Save
    ggsave(myplot, filename = paste(args$outprefix, tag, "png", sep="."), ...)
    ggsave(myplot, filename = paste(args$outprefix, tag, "svg", sep="."), ...)
}

richness_plot(plantdata, "plants", width=5, height=3)
richness_plot(soildata, "soil", width=3.345, height=3)
