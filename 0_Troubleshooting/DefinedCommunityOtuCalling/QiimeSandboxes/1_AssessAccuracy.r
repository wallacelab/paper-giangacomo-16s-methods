#! /usr/bin/Rscript

# Convert QIIMe output to Phyloseq and determine how well defined community calling went

library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-b", "--biom", help="Text-formatted BIOM file of microbiome observations from QIIME")
parser$add_argument("-x", "--taxonomy", help="Taxonomy key for sequences in the BIOM file")
parser$add_argument("-k", "--keyfile", help="QIIME-formatted keyfile of sample data")
parser$add_argument("-t", "--targets", help="File of target OTUs to go for (with Genus and OTU columns)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/0_Troubleshooting/DefinedCommunityOtuCalling/QiimeSandboxes/1_AssignOtus/')
# args=parser$parse_args(c("-b","clustered-seqs.0.99.biom.txt", "-x", "/home/jgwall/Projects/0_RawData/Silva_132_release/majority_taxonomy_7_levels.99.txt", "-o",'99_tmp.png', "-t", "../0_atcc_species.99.txt", "-k", "../defined_community_test_keyfile.tsv"))

# Libraries
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(psadd)

# Load data
cat("Loading data\n")
biom = read.delim(args$biom, skip=1, check.names=F, row.names=1)
taxonomy = read.delim(args$taxonomy, row.names=1, header=FALSE)
targets = read.delim(args$targets)
key = read.delim(args$keyfile, row.names=1)[-1,]

cat("Reformatting data\n")

# Reformat taxonomy from SILVA to phyloseq-compatible
taxlevels = strsplit(taxonomy$V2, split=";")
taxlevels = lapply(taxlevels, function(x){
    x = sub(x, pattern="^D_.__", repl="")
    if(length(x) < 7){
        x=c(x, rep("unknown", times=7-length(x)))   # Pad out any where classification stops above Species. (from QIIME assignment)
    }
    names(x) = c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species")
    return(x)
})
taxlevels = do.call(rbind, taxlevels)
rownames(taxlevels) = rownames(taxonomy)

# Convert data to formats as required by phyloseq
biom = otu_table(as.matrix(biom), taxa_are_rows=TRUE)
taxonomy = tax_table(as.matrix(taxlevels))
key = sample_data(type.convert(key))    # type.convert changes things to factors, integer, etc. instead of the default "character" 

# Merge into a single object
mydata = phyloseq(biom, taxonomy, key)



# ##############
# OTU count report
# ##############

cat("Outputting OTU counts\n")
otu_counts = taxa_sums(mydata)
otu_summary = data.frame(as.data.frame(tax_table(mydata)), otu=names(otu_counts), count=otu_counts)
otu_summary = otu_summary[order(otu_summary$count, decreasing=TRUE),]
write.table(otu_summary, file=paste(args$outprefix, ".otu_summary.txt", sep=""), row.names=F, col.names=T, sep='\t', quote=F)



# ###############
# Summary Graphics Calculations
# ###############

mytaxa = as.data.frame(tax_table(mydata))
otu_counts = taxa_sums(mydata)
master = data.frame(counts=otu_counts, mytaxa, otu=names(otu_counts))

# Create a data table summing things by genus
by_genus = master
by_genus$Genus[! by_genus$Genus %in% targets$Genus ] = "Other"
by_genus$Class[! by_genus$Genus %in% targets$Genus ] = "Other"
by_genus = aggregate(counts ~ Class + Genus, data=by_genus, FUN=sum)
by_genus$Level = "Genus"

# Same, but by OTU/species
by_otu = master
speciesmatch = match(by_otu$otu, targets$OTU)
by_otu$Species = targets$Species[speciesmatch]
by_otu$Species[is.na(by_otu$Species)] = "Other"
by_otu$Class[by_otu$Species == "Other"] = "Other"
by_otu = aggregate(counts ~ Class + Species, data=by_otu, FUN=sum)
by_otu$Level = "OTU"


# ##############
# Stacked barplot of how many reads assigned to correct taxon by Genus and by exact OTU
# ##############

cat("Creating plots\n")

# Combine into one frame for plotting stacked plots
combined = rbind(by_genus[,c("Class", "counts", "Level")], by_otu[,c("Class", "counts", "Level")])
combined$Class = factor(combined$Class)
combined$Class = relevel(combined$Class, ref="Other")

# Calculate fraction mapped for each
totals = split(combined, combined$Level)
totals = lapply(totals, function(x){
    sum(x$counts[x$Class != "Other"])/sum(x$counts) # Fraction non-Other counts
})
totals=do.call(rbind, totals)
totals = data.frame(Level=rownames(totals), fraction=totals[,1])
totals$fraction = paste(round(totals$fraction * 100, digits=1), "%", sep="")

# Make color palette so "Other" is grayed out
colors = scales::hue_pal()(length(levels(combined$Class)))
colors[1] = "darkgray"
names(colors) = levels(combined$Class)

# Actual plot
stacked_plot = ggplot(combined, mapping=aes(x=Level, y=counts)) + 
   geom_bar(position='fill', stat='identity', color="black", mapping=aes(fill=Class)) +
   scale_fill_manual(values=colors) +
   geom_text(data=totals, mapping=aes(x=Level, y=1.1, label=fraction), size=8) +
   ggtitle("Fraction reads clustered to target")

   
# ##############
# Staggered barplot of how much in each OTU
# ##############

# Redo levels so plotted in preferred order (largest to smallest, roughly)
by_otu$Species = sub(by_otu$Species, pattern=" ATCC.+", repl="")    # Remove ATCC number for easier plotting
no_other = subset(by_otu, by_otu$Species !="Other")
by_otu$Species = factor(by_otu$Species, levels = c("Other", no_other$Species[order(no_other$counts, decreasing=T)]))

# Class to a factor for color plotting
by_otu$Class = factor(by_otu$Class)
by_otu$Class = relevel(by_otu$Class, ref="Other")

# Make color palette so "Other" is grayed out
colors = scales::hue_pal()(length(levels(by_otu$Class)))
colors[1] = "darkgray"
names(colors) = levels(by_otu$Class)

# Height of bar representing 5% of the data
goal = sum(otu_counts) / 20

# Actual plot
side_plot = ggplot(by_otu, mapping=aes(x=Species, y=counts, fill=Class)) +
   geom_bar(stat='identity') +
   scale_fill_manual(values=colors) +
   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
   geom_hline(yintercept=goal, color="black", linetype="dashed") +
   annotate(geom="text", x=nrow(by_otu), y=goal, label="5% of total", vjust=-0.5, hjust=1, fontface="italic") +
   ggtitle("Spread of reads among defined community OTUs")
   
   
   
# Write out plots
png(paste(args$outprefix, ".png", sep=""), width=15, height=5, units='in', res=300)
  grid.arrange(grobs=list(stacked_plot, side_plot), ncol=2, widths=c(1,2))
dev.off()



# ##############
# Krona plot of everything; throws an error b/c not an interactive session, so kept last
# ##############

cat("Plotting KRONA plots\n")

# Krona plot of actual data
plot_krona(mydata, output=paste(args$outprefix, ".krona", sep=""), variable="treatment")
