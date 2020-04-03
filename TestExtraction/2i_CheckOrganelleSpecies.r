#! /usr/bin/Rscript

# Use organelle sequences to confirm the identity of each species

library(phyloseq)
library(argparse)
library(ggplot2)
library(gridExtra)
options(stringsAsFactors=F)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="RDS file with saved phyloseq object for analysis")
parser$add_argument("-n", "--top-otus", type='integer', default=10, help="Merge all but the top N OTUs for output display")
parser$add_argument("-o", "--outfile", help="Graphical output file")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 10 Mohsen Final Data/2_Analysis_clean/')
# args=parser$parse_args(c("-i","2b_filtered_data.phyloseq.RDS", "-o",'99_tmp.png'))

# Load data
mydata=readRDS(args$infile)
taxonomy = as.data.frame(tax_table(mydata))

# Pull out names of OTUs classified as organelles
mito_names = rownames(subset(taxonomy, taxonomy$Family=='Mitochondria'))
chloro_names = rownames(subset(taxonomy, taxonomy$Order=='Chloroplast'))

# Helper function to take data and organelle names and return a graphic of their abundance
plot_organelles = function(mytable, taxa_names, top_n=10, title=""){

    # Subset out the target organelles
    targets = prune_taxa(taxa_names, mytable)
    
    # If have too many to display, merge all but the top N
    if(length(targets) > top_n){
        counts = sort(rowSums(otu_table(targets)))
        to_merge = names(counts)[1:(length(counts)-top_n)]
        merged = merge_taxa(targets, to_merge)
        taxa_names(merged)[taxa_names(merged) %in% to_merge] = "all_else"
    }else{
        merged = targets
    }
    
    
    # Normal barplot
    plot1 = plot_bar(merged, fill='Genus', title=title) +
        facet_wrap(~sample.type, scales="free", nrow=1) +
        theme(axis.text.x=element_text(vjust=0.5))
        
    # Have to hack it to do a normalized barplot; code is taken from plot_bar above
    melted = psmelt(merged)
    plot2  = ggplot(melted, aes_string(x = "Sample", y = "Abundance", fill = 'Genus')) + 
        geom_bar(stat = "identity", position = "fill", color = "black") + 
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
        ggtitle(title) +
        facet_wrap(~sample.type, scales="free", nrow=1) 
    
    
    return(list(plot1, plot2))
}


mitoplot = plot_organelles(mydata, mito_names, top_n = args$top_otus, title = "Mitochondria")
chloroplot = plot_organelles(mydata, chloro_names, top_n = args$top_otus, title = "Chloroplast")

# Plot output
png(args$outfile, width=20, height=10, units='in', res=300)
    grid.arrange(grobs=c(mitoplot, chloroplot), nrow=2)
dev.off()
