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
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 10 Mohsen Final Data/2_Analysis_clean/')
# args=parser$parse_args(c("-i","2b_filtered_data.phyloseq.RDS", "-o",'99_tmp.png'))


# Load data
cat("Loading data for alpha diversity analysis\n")
mydata = readRDS(args$infile)

# Create a new combined metadata column of sample type + treatment
metadata = sample_data(mydata)
metadata$sample_type_and_treatment = paste(metadata$sample.type, metadata$treatment, sep="~")
sample_data(mydata) = metadata

# Plot
cat("Plotting alpha diversity\n")
plots = lapply(unique(metadata$sample.type), function(my_sample){
    tokeep = rownames(metadata)[metadata$sample.type==my_sample]
    subdata = prune_samples(mydata, samples=tokeep) # subset_samples would be clearer but kept having an error
    plot_richness(subdata, measures=c("Observed", "Chao1", "Shannon"), x="treatment", 
              color="treatment") +
        ggtitle(my_sample)
})
png(paste(args$outprefix, ".alpha_diversity.png", sep=""), width=10, height=5*length(plots), units='in', res=300)
    grid.arrange(grobs=plots, nrow=length(plots), ncol=1)
dev.off()



# Make a plot of how much each treatment captures all possible OTUs
#   Stacked barplot with "shared" and "unique", split by sample type

# Collapse 
cat("Plotting shared & unique OTUs\n")
merged = merge_samples(mydata, group="sample_type_and_treatment")

# Split by sample type (don't want to compare soil versus leaf)
mergecounts = as.data.frame(otu_table(merged))
key=strsplit(rownames(mergecounts), split="~") # Split back into treatment and sample type
key = as.data.frame(do.call(rbind, key))
names(key) = c("sample.type", "treatment")

# Split by treatment
splitcounts = split(mergecounts, key$sample.type)

# Convert to presence/absence
presence = lapply(splitcounts, function(x){
  x[x>1] = 1  # Take any count > 1 and turn to 1
  return(x)
})

# Get shared/unique OTUs in each
shared = lapply(presence, function(mypresence){
  
  # Determine which OTUs are absent from this group, unique to one treatment, or shared among them
  num_present = colSums(mypresence)
  is_absent = num_present==0 # Not used, but good for sanity-checking
  is_unique = num_present==1
  is_shared = num_present>=2
  
  # Calculate shared & unique OTUs for each treatment
  unique_counts = rowSums(mypresence[,is_unique, drop=FALSE])
  shared_counts = rowSums(mypresence[,is_shared, drop=FALSE])
  
  # Create a final data frame 
  tallied = strsplit(row.names(mypresence), split="~")
  tallied = as.data.frame(do.call(rbind, tallied))
  names(tallied) = c('sample.type', 'treatment')
  myshared = data.frame(tallied, otu_count=shared_counts, count_type='shared')
  myunique = data.frame(tallied, otu_count=unique_counts, count_type='unique')
  
  return(rbind(myshared, myunique))
})
shared=do.call(rbind, shared)
shared$count_type = factor(shared$count_type, levels= c('unique', 'shared'))   # Put in the order I want for plotting

# Plot stacked barplot of shared & unique sequences
ggplot(shared, mapping=aes(x=treatment, y=otu_count, fill=treatment, alpha=factor(count_type))) +
  geom_bar(stat='identity', position='stack') + 
  facet_wrap( ~ sample.type, scales='free') + # Facet into one plot per sample type
  theme_bw() + 
  theme(strip.background  = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank()) + 
  scale_alpha_manual(values=c(0.6, 1))  # Fix alpha values so aren't too transparent
ggsave(paste(args$outprefix, ".unique_otus.png", sep=""))


# Plot of fraction of all OTUs each method captures
fraction_total = lapply(presence, function(mypresence){
  
  # Since presence is 0/1 matrix, easy to get fraction of total
  mypresence = subset(mypresence, select = colSums(mypresence) != 0)    # Remove any OTUs not present in any of these samples
  fraction_present = rowSums(mypresence) / ncol(mypresence) # Get fraction of total remaining OTUs in each method
  
  # Create a final data frame 
  tallied = strsplit(row.names(mypresence), split="~")
  tallied = as.data.frame(do.call(rbind, tallied))
  names(tallied) = c('sample.type', 'treatment')
  mycounts = data.frame(tallied, fraction_total=fraction_present)
    
  return(mycounts)
})
fraction_total=do.call(rbind, fraction_total)

# Plot stacked barplot of shared & unique sequences
ggplot(fraction_total, mapping=aes(x=treatment, y=fraction_total, fill=treatment)) +
  geom_bar(stat='identity') + 
  facet_wrap( ~ sample.type, scales='free') + # Facet into one plot per sample type
  theme_bw() + 
  theme(strip.background  = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())
ggsave(paste(args$outprefix, ".fraction_total_otus.png", sep=""))
