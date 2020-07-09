#! /usr/bin/Rscript

# Identify which bacterial taxa are discrminated against by each group

library(argparse)
library(DESeq2)
library(ggplot2)
library(phyloseq)
options(stringsAsFactors=F)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="RDS file with saved phyloseq object for analysis")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-a", "--alpha", type="double", default=0.05, help="Significance cutoff for differentially expressed taxa")
parser$add_argument("-r", "--reference", default="PowerSoil", help="Reference treatment to contrast against")
parser$add_argument("-l", "--levels", nargs="*", default="Phylum", help="Space-separated list of taxonomic levels to collapse and test at")
parser$add_argument("-z", "--fix-zeros", default=FALSE, action='store_true', help="Adds 1 to all OTU counts to prevent DEseq from ignoring any with 0s")
parser$add_argument("-t", "--type", choices=c("extraction", "amplification"), default="extraction", help="Which experiment set this analysis belongs to")
parser$add_argument("-s", "--sample-type", nargs="*", help="Subset data to just these sample types")
parser$add_argument("-m", "--mean-fits", nargs="*", help="Use the 'mean' option to fit dispersion estimates for these taxonomix levels (default is parametric fit)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2020 03 Consolidated Pipeline/')
# args=parser$parse_args(c("-i","TestPrimers/2_Analysis/2f_otu_table.no_organelles.RDS", "-s", "Soil 1", "Soil 2", "-t", "amplification", 
#     "-r", "Universal", "-l", "Domain", "Phylum", "-m", "Domain", "-o",'99_tmp'))

# Load data
cat("Loading data from",args$infile,"\n")
source("StandardizeLabels.r")
mydata = standardize_labels(readRDS(args$infile), type=args$type)

# Subset to required sample types
if(!is.null(args$sample_type)){
    cat("\tSubsetting to just samples of type",args$sample_type,"\n")
    tokeep = data.frame(sample_data(mydata))$sample.type %in% args$sample_type
    mydata = subset_samples(mydata, tokeep)
}

# Collapse at different levels
cat("\tCollapsing at taxonomic ranks",args$levels,"\n")
collapsed = lapply(args$levels, function(mylevel){
    grouped = tax_glom(mydata, taxrank=mylevel)
    if(args$fix_zeros && any(otu_table(grouped)==0) ){
        otu_table(grouped) = otu_table(grouped)+1
    }
    return(grouped)
})
names(collapsed) = args$levels

# Conver to DEseq2 and analyze
cat("Running DEseq2 analysis")
de.data = lapply(collapsed, FUN = phyloseq_to_deseq2, design= ~ sample.type + treatment)
de.analysis = lapply(names(de.data), function(level){
    myfit="parametric"
    if(level %in% args$mean_fits){
        myfit="mean"
    }
    cat("\tFitting",myfit,"error dispersion for level",level,"\n")
    DESeq(de.data[[level]], test="Wald", fitType="mean")
})

# Compile results comparing against the supplied reference
treatments = levels(sample_data(mydata)$treatment)
treatments = treatments[treatments != args$reference]
results = lapply(1:length(de.analysis), function(level){
    source_data = collapsed[[level]]
    anal=de.analysis[[level]]
    results.list = list()
    for(i in 1:length(treatments)){
        # Extract results
        mytreat = treatments[i]
        myresults = results(anal, contrast=c('treatment', args$reference, mytreat), cooksCutoff = FALSE)
        mytaxa = tax_table(source_data)[rownames(myresults),]
        myresults = data.frame(reference=args$reference, contrast=mytreat, cbind(myresults, mytaxa))
        
        # Add to list
        results.list[[i]] = myresults
        names(results.list)[i] = mytreat
    }
    do.call(rbind, results.list)
})
# Add taxonomic level to data frame and combine
for(i in 1:length(args$levels)){
    results[[i]] = data.frame(collapse_by=args$levels[i], results[[i]])
}
results = do.call(rbind, results)
rownames(results) = 1:nrow(results) # For simplicity

# Subset to significant ones
significant = subset(results, results$padj <= args$alpha)

# Write out
outtext = paste(args$outprefix, ".txt", sep="")
cat("Writing significant results to", outtext,'\n')
write.table(significant, file=outtext, sep='\t', quote=F, row.names=F, col.names=T)

# ###########
# Make a hierarchical graphic to show how much is distorted
# ###########

cat("Plotting hierarchical rectangle plot\n")
taxa = as.data.frame(tax_table(mydata))

# Helper function to make plotting data
get_plot_data = function(equalize=FALSE, contrast=NULL){
    
    # Get counts of each taxon; 'equalize' shows all smallest levels equally instead based on read counts
    taxa_counts = taxa_sums(mydata)
    if(equalize){
        taxa_counts = rep(1, length(taxa_counts))
    }
    
    taxa_levels = names(taxa)[colSums(is.na(taxa))<nrow(taxa)] # Only levels not already collapsed
    rect_data = lapply(taxa_levels, function(level){
        index = which(names(taxa)==level)
        tax_strings = apply(taxa[,1:index,drop=F], FUN=paste, MARGIN=1, collapse=" ")
        
        # Aggregate into total
        sums = aggregate(taxa_counts, by=list(tax_strings), FUN=sum)
        names(sums) = c("taxon", "count")
        cumulative = cumsum(sums$count)
        
        # Create parameters for graphing as rectangles
        sums$left = index
        sums$right = index+1
        sums$bottom = c(0, cumulative[-length(cumulative)])
        sums$top = cumulative
        
        # Return compiled data frame
        return(sums)
    })
    rect_data = do.call(rbind, rect_data)
    
    # Determine which are significant
    sig_taxa = apply(significant[,names(taxa)], FUN=paste, MARGIN=1, collapse=" ")
    sig_taxa = sub(sig_taxa, pattern=" NA.+", repl="") # Removing NAs from collapsed levels
    increasing = sig_taxa[significant$log2FoldChange < 0] # log2 change < 0 means contrast (denominator) is bigger
    decreasing = sig_taxa[significant$log2FoldChange > 0] # log2 change > 0 means contrast (denominator) is smaller
    if(!is.null(contrast)){
        sig_taxa = sig_taxa[significant$contrast == contrast]
    }
    rect_data$taxon = sub(rect_data$taxon, pattern=" NA.+", repl="") 
    rect_data$significant = rect_data$taxon %in% sig_taxa
    rect_data$contrast = contrast
    
    # Get directionality
    rect_data$direction="none"
    rect_data$direction[rect_data$significant & rect_data$taxon %in% increasing] = "increasing"
    rect_data$direction[rect_data$significant & rect_data$taxon %in% decreasing] = "decreasing"
    rect_data$direction=factor(rect_data$direction, levels=c("none", "increasing", "decreasing"))
    
    return(rect_data)
}

# Helper function to get number and percent significant
get_percent_data = function(plotdata){
    levels = split(plotdata, plotdata$left) # Split by taxonomic level (coded by where rectangle's left edge starts)
    taxa_counts = lapply(levels, function(l){
        sum(l$significant)
    })
    taxa_fractions =  lapply(levels, function(l){
        fraction = sum(l$count[l$significant]) / sum(l$count)
    })
    
    # Combine
    taxa_counts = do.call(rbind, taxa_counts)
    taxa_fractions = do.call(rbind, taxa_fractions)
    percent_data = data.frame(level = as.numeric(rownames(taxa_counts)), taxa_counts, taxa_fractions, contrast=unique(plotdata$contrast))
    if(length(unique(percent_data$contrast)) != 1){
        warning("More than one contrast when calculating percent for plot label")
    }
        
    return(percent_data)
}

# Helper function to plot the hierarchical rectangles
plot_rects = function(rects, percents, title="", reference="reference"){
    # Plot values for later reference
    xmin=min(rects$left)
    xmax=max(rects$right)
    ymin=min(rects$bottom)
    ymax = max(rects$top)
    
    # Data for # and % OTUs
    percents$percent = paste(round(percents$taxa_fractions * 100), "%", sep="")
    
    # Visual borders around plots
    borders = data.frame(left=xmin, right=xmax, bottom=ymin, top=ymax, contrast = unique(rects$contrast))
    
    # Taxa names
    tax_names = subset(rects, (left==2) & (contrast == sort(unique(contrast))[1]) )
    tax_names = tax_names[order(tax_names$count, decreasing=T),][1:10,] # Only label top 10 phyla
    tax_names$taxon = sub(tax_names$taxon, pattern="Bacteria ", repl="")
    tax_names$y = (tax_names$top + tax_names$bottom) / 2
    
    # Actual plot
    ggplot() + 
        theme_void() + 
        geom_rect(rects, mapping=aes(xmin=left, xmax=right, ymin=bottom, ymax=top, fill=direction), size=0.1, color='black') +
        geom_rect(borders, mapping=aes(xmin=left, xmax=right, ymin=bottom, ymax=top), fill=NA, size=0.5, color='black') +  # Border around plots
        scale_x_continuous(breaks=1:ncol(taxa) + 0.5, labels=names(taxa)) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y = element_blank()) +
        facet_grid(~contrast) + 
        scale_fill_manual(values=c('white', 'royalblue', 'red3'),    # Color scale
                          labels=c('No change', paste('More than', reference),paste('Less than',reference))) +
        theme(legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=6),
              legend.key.size = unit(0.5, "cm")) +     # Legend adjustments
        theme(strip.text=element_text(size=10)) + # Adjust treatment labels on top
        geom_text(percents, mapping=aes(x=level, label=taxa_counts), y=ymax * -0.025, nudge_x=0.5, size=2) +  # Count labels
        geom_text(percents, mapping=aes(x=level, label=percent), y=ymax * -0.055, nudge_x=0.5, size=1.5) +    # Percentage labels
        ylim(ymax * -0.03, ymax) +
        geom_text(tax_names, mapping=aes(label=taxon, x=left, y=y), hjust="right", vjust="center", nudge_x = -1.5, size=2.5) + # Phylum labels
        theme(plot.margin=unit(c(0,0,0,3), "cm")) + # Adjust margins
        coord_cartesian(clip = 'off')  # Required so phylum labels don't get cut off
        
        
}

# Plot
rawdata = lapply(treatments, function(t){get_plot_data(equalize=FALSE, contrast=t)})
percent_data = lapply(rawdata, get_percent_data)
rawplot = plot_rects(do.call(rbind, rawdata), do.call(rbind, percent_data), reference=args$reference)

# Save 
ggsave(rawplot , file=paste(args$outprefix, ".png", sep=""), width=5, height=3, dpi=300)
ggsave(rawplot , file=paste(args$outprefix, ".svg", sep=""), width=5, height=3, dpi=300)
