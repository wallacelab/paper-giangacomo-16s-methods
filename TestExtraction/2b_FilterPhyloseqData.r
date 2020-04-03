#! /usr/bin/Rscript

# Fitler the Phyloseq object based on OTU and sample statistics.

# Libraries
library(argparse)
library(phyloseq)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Phyloseq object saved from 2a_ConvertQiimeToPhyloseqObject.r")
parser$add_argument("-s", "--sample-depth", type='integer', default=1, help="Minimum number of reads to keep samples")
parser$add_argument("-d", "--otu-depth", type='integer', default=1, help="Minimum number of samples in which an OTU has to appear to be kept. (OTUs must pass this AND the prevalence criterion)")
parser$add_argument("-p", "--otu-prevalence", type='integer', default=1, help="Minimum number of reads for an OTU for it to be keptsamples in which an OTU has to appear to be kept. (OTUs must pass this AND the depth criterion)")
parser$add_argument("-o", "--outfile", help="New output phyloseq object")
parser$add_argument("--outgraphics", help="PNG file to write graphics to")
parser$add_argument("--group-by", help="Taxonomic rank to group results by before applying filters (optional)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/2_Analysis/')
# args=parser$parse_args(c("-i","2a_combined_data.phyloseq.Rdata", "-s", "1000", "-d", "3", "-p", "2", "-o", '99_tmp.Rdata'))


# Load data
cat("Filtering Phyloseq data object from",args$infile,"\n")
mydata = readRDS(args$infile)   # Loads phyloseq data into variable 'mydata'

# Group by phylogenetic rank if requested
if(!is.null(args$group_by)){
    cat("\tGrouping data by", args$group_by,"\n")
    mydata = tax_glom(mydata, taxrank=args$group_by)
}


# ###################
# Filtering
# ###################

# Helper function to display percentages
get_percent=function(numerator, denominator){
    percent = numerator/denominator * 100
    percent = round(percent, 2)
    percent = paste(percent, "%", sep="")
    return(percent)
}

# Determine which samples to keep
samples = data.frame(id = sample_names(mydata), counts = sample_sums(mydata))
samples_tokeep = samples$id[samples$counts >= args$sample_depth]
cat("\tMinimum sample depth",args$sample_depth,"leaves",length(samples_tokeep),"out of",nrow(samples),"samples (", get_percent(length(samples_tokeep),nrow(samples)),")\n")

# Filter based on samples
filtered = prune_samples(samples_tokeep, mydata)


# Determine which OTUs to keep based on depth
taxa = data.frame(id = taxa_names(filtered), counts = taxa_sums(filtered))
depth_tokeep = taxa$id[taxa$counts >= args$otu_depth]
cat("\tMinimum OTU depth",args$otu_depth,"leaves",length(depth_tokeep),"out of",nrow(taxa),"OTUs (", get_percent(length(depth_tokeep),nrow(taxa)),")\n")

# Determine which OTUs to keep based on prevalence
presence = transform_sample_counts(filtered, fun = function(x){ifelse(x>=1, yes=1, no=0)})    # Converts to a 1/0 matrix of presence/absence
prevalence = data.frame(id = taxa_names(presence), counts = taxa_sums(presence))
prevalence_tokeep = prevalence$id[prevalence$counts >= args$otu_prevalence]
cat("\tMinimum OTU prevalence",args$otu_prevalence,"leaves",length(prevalence_tokeep),"out of",nrow(prevalence),"OTUs (", get_percent(length(prevalence_tokeep),nrow(prevalence)),")\n")

# Find which taxa to keep based on which passed both depth and prevalence criteria (=intersection)
taxa_tokeep = intersect(depth_tokeep, prevalence_tokeep)
cat("\tCombining the OTU filters leaves",length(taxa_tokeep),"out of",nrow(taxa),"OTUs (", get_percent(length(taxa_tokeep),nrow(taxa)),")\n")

# Fitler based on taxa
filtered = prune_taxa(taxa_tokeep, filtered)

cat("\tThe final filtered table retains", get_percent(sum(otu_table(filtered)), sum(otu_table(mydata))), "of all reads\n")



# ###################
# Plot summary graphs
# ###################

png(args$outgraphics, width=1000, height=1500)
par(mfrow=c(3,2), cex=1.2)

# Sample cutoffs
plot(sort(sample_sums(mydata)), main="Sample Depths (pre-filter)", xlab="", ylab="Counts")
abline(h=args$sample_depth, lty='dashed')
plot(sort(sample_sums(filtered)), main="Sample Depths (post-filter)", xlab="", ylab="Counts")

# OTU depths
plot(sort(taxa_sums(mydata)), main="Taxa Depths (pre-filter)", xlab="", ylab="Log Counts", log='y')
abline(h=args$otu_depth, lty='dashed')
plot(sort(taxa_sums(filtered)), main="Taxa Depths (post-filter)", xlab="", ylab="Log Counts", log='y')

# OTU prevalence
raw_presence = transform_sample_counts(mydata, fun = function(x){ifelse(x>=1, yes=1, no=0)})    # Converts to a 1/0 matrix of presence/absence
filtered_presence = transform_sample_counts(filtered, fun = function(x){ifelse(x>=1, yes=1, no=0)})    # Converts to a 1/0 matrix of presence/absence
hist(sort(taxa_sums(raw_presence)), main="Taxa Prevalence (pre-filter)", xlab="# OTUs", ylab="Counts", breaks=50)
abline(v=args$otu_prevalence, lty='dashed')
hist(sort(taxa_sums(filtered_presence)), main="Taxa Prevalence (post-filter)", xlab="# OTUs", ylab="Counts", breaks=50)

dev.off()



# ###################
# Export filtered data
# ###################

# Rename for consistency with original
mydata = filtered

# Need to reroot the phylogenetic tree because there's a chance the root was filtered out (which causes problems with Unifrac stuff downstream)
# Following function is from joey711 on Phyloseq Issue # 597 (https://github.com/joey711/phyloseq/issues/597)
pick_new_outgroup <- function(tree.unrooted){
    require("magrittr")
    require("data.table")
    require("ape") # ape::Ntip
    # tablify parts of tree that we need.
    treeDT <- 
      cbind(
        data.table(tree.unrooted$edge),
        data.table(length = tree.unrooted$edge.length)
      )[1:Ntip(tree.unrooted)] %>% 
      cbind(data.table(id = tree.unrooted$tip.label))
    # Take the longest terminal branch as outgroup
    new.outgroup <- treeDT[which.max(length)]$id
    return(new.outgroup)
  }
old_tree = phy_tree(filtered)
new_root = pick_new_outgroup(old_tree)
new_tree = ape::root(old_tree, outgroup=new_root, resolve.root=TRUE)
phy_tree(filtered) = new_tree

# Write out
cat("Writing filtered data to",args$outfile,"\n")
saveRDS(mydata, file=args$outfile)
