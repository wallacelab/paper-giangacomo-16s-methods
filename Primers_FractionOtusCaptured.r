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
parser$add_argument("--group-by", help="Taxonomic rank to group results by (optional)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2020 03 Consolidated Pipeline/')
# args=parser$parse_args(c("-i","TestPrimers/2_Analysis/2f_otu_table.no_organelles.RDS", "-o",'99_tmp'))


 # Load data
cat("Loading data for fraction of OTUs captured\n")
mydata = readRDS(args$infile)
mydata = prune_samples(mydata, samples= ! sample_data(mydata)$sample.type %in% c("water")) # Filter down to just desired sequences

# Create a new combined metadata column of sample type + treatment
metadata = sample_data(mydata)                            
metadata$sample.type = sapply(as.character(metadata$sample.type), switch,
                              "leaf-maize"="Maize Leaf", 
                              "leaf-soybean"="Soybean Leaf", 
                              "defined-community"="Defined Community", 
                              "soil-clay"="Soil 1",
                              "soil-flowerbed"="Soil 2",
                              NA) # NA catches anything that didn't match
                              
metadata$treatment = sapply(as.character(metadata$treatment), switch,
                         BlockingOligos_v3v4="BO_3/4", 
                         BlockingOligos_v5v7="BO_5/7",
                         BlockingOligos_v5v7_noLinkers="BO_5/7",
                         Discriminating = "Discriminating",
                         PNAs= "PNA", Universal="Universal",
                         NA) # NA catches anything that didn't match

metadata$sample_type_and_treatment = paste(metadata$sample.type, metadata$treatment, sep="~")
sample_data(mydata) = metadata

# Group by phylogenetic rank if requested
if(!is.null(args$group_by)){
    cat("\tGrouping data by", args$group_by,"\n")
    mydata = tax_glom(mydata, taxrank=args$group_by)
}


# ###############
# Shared and Unique OTUs
# ###############

# Split by sample type (don't want to compare soil versus leaf)
merged = merge_samples(mydata, group="sample_type_and_treatment")
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
  total_otus = sum(is_unique | is_shared) # Total possible OTUs in this sample set

  # Create a final data frame
  tallied = strsplit(row.names(mypresence), split="~")
  tallied = as.data.frame(do.call(rbind, tallied))
  names(tallied) = c('sample.type', 'treatment')
  myshared = data.frame(tallied, otu_count=shared_counts, total_possible = total_otus, count_type='Shared')
  myunique = data.frame(tallied, otu_count=unique_counts, total_possible = total_otus, count_type='Unique')

  return(rbind(myshared, myunique))
})
shared=do.call(rbind, shared)
shared$count_type = factor(shared$count_type, levels= c('Unique', 'Shared'))   # Put in the order I want for plotting

# Get value of max possible OTUs in each category
max_possible = unique(shared[,c("sample.type", "total_possible")])
max_possible$label = paste("Total unique =", max_possible$total_possible)
max_possible$label_position = max_possible$total_possible * 1.01

# Determine Y axis label
ylabel="# OTUs"
if(!is.null(args$group_by)){
    ylabel = paste(args$group_by, "count")
}

# Control sample order; (much more complicated than it really should be)
sample.order = c("Defined Community", "Soil 1", "Soil 2", "Maize Leaf", "Soybean Leaf")
shared$sample.type = as.numeric(factor(shared$sample.type, levels=sample.order))
max_possible$sample.type = as.numeric(factor(max_possible$sample.type, levels=sample.order)) 
label_key = setNames(sample.order, 1:length(sample.order))
label_key[label_key=="Defined Community"] = "Defined\nCommunity"

# Plot stacked barplot of shared & unique sequences
myplot = ggplot(shared, mapping=aes(x=treatment, y=otu_count)) +
  geom_bar(stat='identity', position='stack', mapping=aes(fill=treatment, alpha=factor(count_type))) +
  #facet_wrap( ~ sample.type, scales='free') + 
  facet_wrap( ~ sample.type, scales='free', nrow=1, labeller=labeller(sample.type=label_key)) + 
  theme_bw() +
  theme(strip.background  = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        text=element_text(face="bold"), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_alpha_manual(values=c(0.6, 1))  + # Fix alpha values so aren't too transparent
  labs(alpha="OTU category", fill="Primer Set") +
  ylab(ylabel) +
  xlab(NULL) + 
  # Add horizontal line of total possible OTUs in each sample type
  geom_hline(data=max_possible, mapping=aes(yintercept=total_possible), linetype='dashed') +
  geom_text(data=max_possible, mapping=aes(x=0.5, y=label_position, label=label),
            hjust="left", vjust="bottom", size=2.5, fontface='italic')
  
ggsave(myplot, filename = paste(args$outprefix, ".png", sep=""), width=8, height=3)
ggsave(myplot, filename = paste(args$outprefix, ".svg", sep=""), width=8, height=3)
