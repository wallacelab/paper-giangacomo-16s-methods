#! /usr/bin/Rscript

# Plot graphics showing how well (or not) each primer set kept the community as determined by the Universal (Earth Microbiome Project) primers)

# Libraries
library(argparse)
library(ggplot2)
library(phyloseq)
options(stringsAsFactors=F)

# Command-line arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="RDS file containing the phyloseq object to analyze (from step 2a)")
parser$add_argument("-o", "--outprefix", help="Prefix for all output files")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/2_Analysis/')
# args=parser$parse_args(c("-i", "2b_filtered_data.phyloseq.Rdata", "-o", "99_tmp"))

# Load data
cat("Loading phyloseq data\n")
mydata=readRDS(args$infile)


# Split data into plant and non-plant samples
is_plant = get_variable(mydata, "sample.type") %in% c("leaf-maize", "leaf-soybean")
non_plant = get_variable(mydata, "sample.type") %in% c("defined-community", "soil-clay", "soil-flowerbed")
plants = prune_samples(is_plant, mydata)
nonplants = prune_samples(non_plant, mydata) 

# Make my own function to collapse by levels, since taxa_glom() keeps giving integer overflow errors
collapse_taxa=function(phylo, level){
    # Check that the requested phylogenetic level is actually in the data
    if(! level %in% rank_names(phylo)){
        warning("Illegal taxonomic rank requested; returning unchanged object")
        return(phylo)
    }
    
    # Filter out OTUs with 0 reads
    has_reads = which(taxa_sums(phylo)>0)
    phylo = prune_taxa(names(has_reads), phylo)
    
    # Iteratively merge groups
    tempdata = phylo
    mytaxa = unique(as.character(tax_table(tempdata)[,level]))
    for(taxon in sort(mytaxa)){
        # cat("Collapsing", taxon,"\n")
        classifications = as.character(tax_table(tempdata)[,level])
        tomerge = which(classifications == taxon)
        tempdata = merge_taxa(tempdata, tomerge, archetype=1)
    }
    return(tempdata)
}


# Look over some high-level taxonomy to display results
for(level in c("Domain", "Phylum", "Class", "Order")){
    
    cat("\tPlotting distortion barplots at level",level,"\n")
    
    # Plant community data
    subplants = collapse_taxa(plants, level=level)
    subplants = transform_sample_counts(subplants, fun=function(x){x/sum(x)})   # Convert to relative abundance
    p1 = plot_bar(subplants, fill=level, title=paste("Plant community distortion at",level,"level")) + facet_grid(.~sample.type, scales="free_x")
    
    # Nonplant community data
    subnonplants = collapse_taxa(nonplants, level=level)
    subnonplants = transform_sample_counts(subnonplants, fun=function(x){x/sum(x)})   # Convert to relative abundance
    p2 = plot_bar(subnonplants, fill=level, title=paste("Non-plant community distortion at",level,"level")) + facet_grid(.~sample.type, scales="free_x")
    
    ggsave(paste(args$outprefix, level, "plant", "png", sep="."), plot=p1, width=12, height=8)
    ggsave(paste(args$outprefix, level, "nonplant", "png", sep="."), plot=p2, width=12, height=8)
}
