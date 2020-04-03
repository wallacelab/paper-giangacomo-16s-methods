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
mydata = readRDS(args$infile)

# Create a new combined metadata column of sample type + treatment
metadata = sample_data(mydata)
metadata$sample.type = sapply(as.character(metadata$sample.type), switch,
                              "Leaf-Arabidopsis"="Arabidopsis", 
                              "Leaf-Corn"="Maize", 
                              "Leaf-Soybean"="Soybean", 
                              "Soil-Soil SS1"="Soil",
                                NA) # NA catches anything that didn't match
metadata$sample.type = factor(metadata$sample.type, levels=c("Arabidopsis", "Maize", "Soybean", "Soil"))
sample_data(mydata) = metadata


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
# richness_plot(mydata, "all")
    


    
    
    
    
    
    
    
    
    
    
# # Make a plot of how much each treatment captures all possible OTUs
# #   Stacked barplot with "shared" and "unique", split by sample type
# 
# # Collapse 
# cat("Plotting shared & unique OTUs\n")
# merged = merge_samples(mydata, group="sample_type_and_treatment")
# 
# # Split by sample type (don't want to compare soil versus leaf)
# mergecounts = as.data.frame(otu_table(merged))
# key=strsplit(rownames(mergecounts), split="~") # Split back into treatment and sample type
# key = as.data.frame(do.call(rbind, key))
# names(key) = c("sample.type", "treatment")
# 
# # Split by treatment
# splitcounts = split(mergecounts, key$sample.type)
# 
# # Convert to presence/absence
# presence = lapply(splitcounts, function(x){
#   x[x>1] = 1  # Take any count > 1 and turn to 1
#   return(x)
# })
# 
# # Get shared/unique OTUs in each
# shared = lapply(presence, function(mypresence){
#   
#   # Determine which OTUs are absent from this group, unique to one treatment, or shared among them
#   num_present = colSums(mypresence)
#   is_absent = num_present==0 # Not used, but good for sanity-checking
#   is_unique = num_present==1
#   is_shared = num_present>=2
#   
#   # Calculate shared & unique OTUs for each treatment
#   unique_counts = rowSums(mypresence[,is_unique, drop=FALSE])
#   shared_counts = rowSums(mypresence[,is_shared, drop=FALSE])
#   
#   # Create a final data frame 
#   tallied = strsplit(row.names(mypresence), split="~")
#   tallied = as.data.frame(do.call(rbind, tallied))
#   names(tallied) = c('sample.type', 'treatment')
#   myshared = data.frame(tallied, otu_count=shared_counts, count_type='shared')
#   myunique = data.frame(tallied, otu_count=unique_counts, count_type='unique')
#   
#   return(rbind(myshared, myunique))
# })
# shared=do.call(rbind, shared)
# shared$count_type = factor(shared$count_type, levels= c('unique', 'shared'))   # Put in the order I want for plotting
# 
# # Plot stacked barplot of shared & unique sequences
# ggplot(shared, mapping=aes(x=treatment, y=otu_count, fill=treatment, alpha=factor(count_type))) +
#   geom_bar(stat='identity', position='stack') + 
#   facet_wrap( ~ sample.type, scales='free') + # Facet into one plot per sample type
#   theme_bw() + 
#   theme(strip.background  = element_blank(),
#         panel.grid.major.x=element_blank(),
#         panel.grid.minor.x=element_blank()) + 
#   scale_alpha_manual(values=c(0.6, 1))  # Fix alpha values so aren't too transparent
# ggsave(paste(args$outprefix, ".unique_otus.png", sep=""))
# 
# 
# # Plot of fraction of all OTUs each method captures
# fraction_total = lapply(presence, function(mypresence){
#   
#   # Since presence is 0/1 matrix, easy to get fraction of total
#   mypresence = subset(mypresence, select = colSums(mypresence) != 0)    # Remove any OTUs not present in any of these samples
#   fraction_present = rowSums(mypresence) / ncol(mypresence) # Get fraction of total remaining OTUs in each method
#   
#   # Create a final data frame 
#   tallied = strsplit(row.names(mypresence), split="~")
#   tallied = as.data.frame(do.call(rbind, tallied))
#   names(tallied) = c('sample.type', 'treatment')
#   mycounts = data.frame(tallied, fraction_total=fraction_present)
#     
#   return(mycounts)
# })
# fraction_total=do.call(rbind, fraction_total)
# 
# # Plot stacked barplot of shared & unique sequences
# ggplot(fraction_total, mapping=aes(x=treatment, y=fraction_total, fill=treatment)) +
#   geom_bar(stat='identity') + 
#   facet_wrap( ~ sample.type, scales='free') + # Facet into one plot per sample type
#   theme_bw() + 
#   theme(strip.background  = element_blank(),
#         panel.grid.major.x=element_blank(),
#         panel.grid.minor.x=element_blank())
# ggsave(paste(args$outprefix, ".fraction_total_otus.png", sep=""))
