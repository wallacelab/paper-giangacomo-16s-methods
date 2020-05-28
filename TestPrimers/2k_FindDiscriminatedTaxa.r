#! /usr/bin/Rscript

# Identify which bacterial taxa are discrminated against by each group

library(ape)
library(phyloseq)
library(argparse)
library(DESeq2)
library(ggplot2)
options(stringsAsFactors=F)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="RDS file with saved phyloseq object for analysis")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-a", "--alpha", type="double", default=0.05, help="Significance cutoff for differentially expressed taxa")
parser$add_argument("-r", "--reference", default="PNAs", help="Reference treatment to contrast against")
parser$add_argument("-l", "--levels", nargs="*", default="Phylum", help="Space-separated list of taxonomic levels to collapse and test at")
parser$add_argument("-z", "--fix-zeros", default=FALSE, action='store_true', help="Adds 1 to all OTU counts to prevent DEseq from ignoring any with 0s")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2020 03 Consolidated Pipeline/TestPrimers/2_Analysis/')
# args=parser$parse_args(c("-i","2f_otu_table.no_organelles.RDS", "-o",'99_tmp', '-l', 'Phylum', 'Genus'))


# Load data
cat("Loading data from",args$infile,"\n")
mydata=readRDS(args$infile)

# Collapse at different levels
cat("\tCollapsing at taxonomic ranks",args$levels,"\n")
collapsed = lapply(args$levels, function(mylevel){
    grouped = tax_glom(mydata, taxrank=mylevel)
    if(args$fix_zeros && any(otu_table(grouped)==0) ){
        otu_table(grouped) = otu_table(grouped)+1
    }
    return(grouped)
})

# Conver to DEseq2 and analyze
cat("Running DEseq2 analysis")
de.data = lapply(collapsed, FUN = phyloseq_to_deseq2, design= ~ sample.type + treatment)
de.analysis = lapply(de.data, FUN=DESeq, test="Wald", fitType="parametric")

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
    if(!is.null(contrast)){
        sig_taxa = sig_taxa[significant$contrast == contrast]
    }
    rect_data$taxon = sub(rect_data$taxon, pattern=" NA.+", repl="") 
    rect_data$significant = rect_data$taxon %in% sig_taxa
    rect_data$contrast = contrast
    
    return(rect_data)
}

# Helper function to plot the hierarchical rectangles
plot_rects = function(rects, title=""){
    ggplot(rects, mapping=aes(xmin=left, xmax=right, ymin=bottom, ymax=top, fill=significant)) + 
        geom_rect(size=0.01, color='black') +
        scale_x_continuous(breaks=1:ncol(taxa) + 0.5, labels=names(taxa)) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y = element_blank()) +
        facet_grid(~contrast) + 
        ggtitle(title)
}

# Plot
rawdata = lapply(treatments, function(t){get_plot_data(equalize=FALSE, contrast=t)})
equaldata = lapply(treatments, function(t){get_plot_data(equalize=TRUE, contrast=t)})
rawplot = plot_rects(do.call(rbind, rawdata), title="Raw counts")
equalplot = plot_rects(do.call(rbind, equaldata), title="Equalized counts")

# Save 
ggsave(rawplot, file=paste(args$outprefix, ".raw_counts.png", sep=""), width=5, height=3, dpi=300)
ggsave(equalplot, file=paste(args$outprefix, ".equal_counts.png", sep=""), width=5, height=3, dpi=300)
