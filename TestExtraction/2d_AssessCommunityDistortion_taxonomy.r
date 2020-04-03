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
    
    # Community data
    subdata = collapse_taxa(mydata, level=level)
    subdata = transform_sample_counts(subdata, fun=function(x){x/sum(x)})   # Convert to relative abundance
    p1 = plot_bar(subdata, fill=level, title=paste("Community distortion at",level,"level")) + facet_grid(.~sample.type, scales="free_x")
    
    ggsave(paste(args$outprefix, level, "png", sep="."), plot=p1, width=12, height=8)
}
