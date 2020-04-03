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
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 10 Mohsen Final Data/2_Analysis/')
# args=parser$parse_args(c("-i", "2b_filtered_data.phyloseq.RDS", "-o", "99_tmp", "-r", "1000" ))

cat("Assessing community distortion with beta diversity metrics\n")

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
distances=list(weighted=as.matrix(weighted), unweighted=as.matrix(unweighted), bray=bray)

# Standardize treatment names
metadata$treatment_std = sapply(as.character(metadata$treatment), switch,
                                ExtracNAmp="ExtractNAmp", 
                                MoBioPowerSoil="PowerSoil",
                                QiagenDNeasyPlant="DNeasyPlant",
                                ZymoEasyDna = "EasyDNA",
                                KazuBuffer = "KazuBuffer",
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

# Helper function for plotting and saving MDS plots
plot.mds = function(mydistances, outfile, type="", ...){
  
  plots = lapply(names(mydistances), function(metric){
    mydata=mydistances[[metric]]
    p = qplot(data=mydata, x=PC1, y=PC2, size=10, alpha=0.5, main=paste(type, metric, sep=" - "), ...) +
        xlab(paste("PC1 (", round(mydata$PC1_var[1]*100, 1), "%)", sep="")) + 
        ylab(paste("PC2 (", round(mydata$PC2_var[2]*100, 1), "%)", sep="")) +
        geom_text(mapping=aes(label=rownames(mydata)), size=2, color='black')
    return(p)
    })
    
  png(outfile, width=20, height=5, units='in', res=300)
    grid.arrange(grobs=plots, ncol=length(plots))
  dev.off()
}



targets  = subMDS(distances, metadata, metadata$sample) # Should take all samples

# Plot by sample type
plot.mds(targets, paste(args$outprefix, ".all_by_sample_type.png", sep=""), color=sample.type)
plot.mds(targets, paste(args$outprefix, ".all_by_treatment.png", sep=""), color=treatment_std)


# Loop over sample types and plot individualls
for(mytype in unique(metadata$sample.type)){
  samples = rownames(metadata)[metadata$sample.type == mytype ]
  targets = subMDS(distances, metadata, samples)
  plot.mds(targets, paste(args$outprefix, ".", mytype, ".png", sep=""), color=treatment_std, type=mytype)
}


