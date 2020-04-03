#! /usr/bin/Rscript

library(argparse)
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(rbiom)
options(stringsAsFactors=F)

# Command-line arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="RDS file containing the phyloseq object to analyze (from step 2a)")
parser$add_argument("-r", "--rarefaction", type="integer", help="Text file of Weighted UniFrac distances")
parser$add_argument("--force-rarefaction-level", type="logical", default=FALSE, help="Use the specified rarefaction level even if it is lower than the lowest sample depth. (Default is to use the lowest sample depth if it's higher than the value given by --rarefaction.)")
parser$add_argument("-o", "--outprefix", help="Prefix for all output files")
args=parser$parse_args()
# setwd('/home/jgwall/Documents/Papers/16sMethodsDevelopment_Cecelia_Mohsen/Pipeline/')
# args=parser$parse_args(c("-i", "TestExtraction/2_Analysis/2f_otu_table.no_organelles.RDS", "-o", "99_tmp", "-r", "2000" ))

cat("Assessing community distortion with beta diversity metrics\n")
set.seed(1)

# Load phyloseq data
mydata = readRDS(args$infile)
metadata = sample_data(mydata)
mytable=otu_table(mydata)
mytree=phy_tree(mydata)

# Check if the specified rarefaction is lower than the smallest sample depth and change if it is (and user didn't overrule this behavior)
if(!args$force_rarefaction_level && args$rarefaction < min(sample_sums(mydata))){
    cat("\tNote: Specified rarefaction is level is less than the minimum sample depth, so minimum sample depth will be used instead.\n")
    args$rarefaction = min(sample_sums(mydata))
}

# Rarefy matrix with rbiom
cat("Calculating distance metrics\n")
rarefied = rbiom::rarefy(mytable, depth=args$rarefaction)  # Specify rbiom:: to be absolutely certain we don't use vegan's function of the same name
cat("\tRemoved", ncol(mytable) - ncol(rarefied), "samples for having fewer than", args$rarefaction,"total reads for rarefaction\n")

# Adjust metadata to reflect fewer samples (potentially)
metadata$sample = rownames(metadata)
metadata = subset(metadata, metadata$sample %in% colnames(rarefied))

# Calculate UniFrac distances with rbiom
weighted = unifrac(rarefied, tree=mytree, weighted=TRUE)
unweighted = unifrac(rarefied, tree=mytree, weighted=FALSE)

# Calculate Bray-Curtis distance matrix with vegan
bray = as.matrix(vegdist(t(rarefied), method='bray'))

# Combined distances into a single list item
distances=list("Weighted UniFrac"=as.matrix(weighted), "Unweighted UniFrac"=as.matrix(unweighted), "Bray-Curtis"=bray)

# Standardize treatment s Sample names
metadata$treatment_std = sapply(as.character(metadata$treatment), switch,
                                ExtracNAmp="ExtractNAmp", 
                                MoBioPowerSoil="PowerSoil",
                                QiagenDNeasyPlant="DNeasyPlant",
                                ZymoEasyDna = "EasyDNA",
                                KazuBuffer = "KazuBuffer",
                                NA) # NA catches anything that didn't match
metadata$sample.type = sapply(as.character(metadata$sample.type), switch,
                              "Leaf-Arabidopsis"="Arabidopsis Leaf", 
                              "Leaf-Corn"="Maize Leaf", 
                              "Leaf-Soybean"="Soybean Leaf", 
                              "Soil-Soil SS1"="Soil",
                                NA) # NA catches anything that didn't match


# Helper function to subset and do MDS each time, making a list of output results
subMDS = function(mydistances, metadata, mysamples){
  
  myresults = lapply(mydistances, function(mydist){
    mydist = mydist[mysamples, mysamples]
    myMDS = cmdscale(mydist, eig=T)
    
    # Combine into a single data frame
    mymeta = metadata[mysamples,]
    mymeta$PC1 = myMDS$points[,1]
    mymeta$PC2 = myMDS$points[,2]
    
    # Variance per PC; saving as a data frame column for convenience
    pc_variance = myMDS$eig / sum(myMDS$eig)
    mymeta$PC1_var = pc_variance[1]
    mymeta$PC2_var = pc_variance[2]
    return(mymeta)
  })
  
  return(myresults)
  
}


# #############
# Main text figure - Construct with grid.arrange()
# #############

# Helper function for plotting MDS plots; returns a single ggplot item
plot.mds = function(mydata, type="", metric="", legend.title=NULL, ...){
  myplot = ggplot(data=mydata, mapping=aes(x=PC1, y=PC2, ...), ...) +
      xlab(paste("PC1 (", round(mydata$PC1_var[1]*100, 1), "%)", sep="")) +
      ylab(paste("PC2 (", round(mydata$PC2_var[2]*100, 1), "%)", sep="")) +
      geom_point(size=4, alpha=0.65) +
      #ggtitle(paste(type, metric, sep=" - ")) +
      ggtitle(type) + 
      theme_bw() +
      theme(aspect.ratio=1, plot.title = element_text(size=10, face="bold"), 
            axis.title = element_text(size=10, face="bold"), 
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.title=element_text(size=10, face="bold"))
  if(! is.null(legend.title)){
      myplot = myplot + labs(color=legend.title)
  }
  return(myplot)
}

# Get MDS plots of everything
alldata  = subMDS(distances, metadata, metadata$sample)
all_weighted = plot.mds(alldata[['Weighted UniFrac']], color=sample.type, type="All", metric="Weighted UniFrac", legend.title="Sample Type")
all_bray = plot.mds(alldata[['Bray-Curtis']], color=sample.type, type="All", metric="Bray-Curtis", legend.title="Sample Type")
all_plots = list(all_weighted, all_bray)
# top_plots = list(all_weighted) # Don't need Bray-Curtis

# Sample-specific plots
metric="Weighted UniFrac"
sample_types = c("Soil","Arabidopsis Leaf","Maize Leaf","Soybean Leaf")
sample_plots = lapply(sample_types, function(mytype){
  samples = rownames(metadata)[metadata$sample.type == mytype ]
  targets = subMDS(distances, metadata, samples)
  plot.mds(targets[[metric]], color=treatment_std, type=mytype, metric=metric, legend.title="Extraction Method") +
    scale_color_brewer(palette = "Dark2") # Change color scale
})

# # Combine
# all_plots = c(top_plots, sample_plots)

# Output Everything together
png(paste(args$outprefix, "all", "png", sep="."), width=8, height=2, units='in', res=300)
  grid.arrange(grobs=all_plots, ncol=2)
dev.off()
svg(paste(args$outprefix, "all", "svg", sep="."), width=8, height=2)
  grid.arrange(grobs=all_plots, ncol=2)
dev.off()




# # Output by-sample plots graphics
png(paste(args$outprefix, "by_sample", "png", sep="."), width=16, height=2, units='in', res=300)
  grid.arrange(grobs=sample_plots, nrow=1)
dev.off()

svg(paste(args$outprefix, "by_sample", "svg", sep="."), width=16, height=2)
  grid.arrange(grobs=sample_plots, nrow=1)
dev.off()
