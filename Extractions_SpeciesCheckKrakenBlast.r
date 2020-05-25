#! /usr/bin/Rscript

library(argparse)
library(ggplot2)
library(gridExtra)
library(phyloseq)
options(stringsAsFactors=F)

# Command-line arguments
parser=ArgumentParser()
parser$add_argument("-b", "--blast-results", nargs="*", help="List of BLAST output files")
parser$add_argument("-k", "--kraken-results", nargs="*", help="List of Kraken output files")
parser$add_argument("-r", "--rds-file", help="RDS file with phyloseq object containing the samples to limit display to")
parser$add_argument("-m", "--min-cutoff", type="double", default=0.01, help="Minimum fraction of total BLAST scores or Kraken hits for a taxon to be displayed")
parser$add_argument("-l", "--taxon-level", default="F", help="Default taxonomic level to show for kraken results (1- or 2-letter abbreviation)")
parser$add_argument("-t", "--taxonomy", help="Taxonomy key for the reference sequneces")
parser$add_argument("--keyfile", help="Keyfile to use for metadata")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2020 03 Consolidated Pipeline/')
# args=parser$parse_args(c("-b", "0_Troubleshooting/CheckSpeciesByBlast/1_M99.blast.txt", "0_Troubleshooting/CheckSpeciesByBlast/1_M100.blast.txt",
#                          "-k", "0_Troubleshooting/KrakenCheckExtractions/0a_kraken_report.M79.txt", "0_Troubleshooting/KrakenCheckExtractions/0a_kraken_report.M80.txt",
#                          "--keyfile", "TestExtraction/16s_extractions_keyfile.tsv", "-t", "~/Projects/0_RawData/Silva_132_release/majority_taxonomy_7_levels.99.txt",
#                          "-o", "99_tmp", "-r", "TestExtraction/2_Analysis/2f_otu_table.no_organelles.RDS"))


cat("Plotting species checks for BLAST and Kraken results\n")

# Load sample keyfile
sample_key = read.delim(args$keyfile, comment.char = "#", row.names=1)

# Pre-format taxonomy
taxonomy = read.delim(args$taxonomy, row.names=1, header=FALSE, col.names=c("id", "tax_string"))
tax_levels = do.call(rbind, strsplit(taxonomy$tax_string, ";"))
colnames(tax_levels) = c("kingdom","phylum","class","order","family","genus","species")
tax_levels = sub(tax_levels, pattern="^D_.__", repl="")    
taxonomy = cbind(taxonomy, tax_levels)
taxonomy$is_organ = taxonomy$family=="Mitochondria" | taxonomy$order=="Chloroplast"

# Identify target samples and make display labels
targets = readRDS(args$rds_file)
target_samples = sample_names(targets)
target_data = sample_data(targets)
target_breaks = paste(target_data$sample.type, target_data$treatment, rownames(target_data), sep="|")
target_labels = paste(target_data$sample.type, target_data$treatment, sep= " - ")
target_labels= sub(target_labels, pattern="^[^-]+-", repl="")

###############
# BLAST Results
###############

cat("Processing BLAST results\n")
blast_columns=c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data = lapply(args$blast_results, read.delim, col.names=blast_columns)

# Process down to hits per reference sequnece (multiple best hits are split evenly)
tally_hits=function(results){
     results = split(results, results$qaccver) # Split by query
     tallied = lapply(results, function(r){
         r = subset(r, r$bitscore == max(r$bitscore)) # Only take best hits
         r$value = 1 / nrow(r)
         return(r)
     })
     tallied = do.call(rbind, tallied)

     # Collapse by reference sequence (I'm sure there's a better way to do this)
     counts = split(tallied, tallied$saccver)
     counts = lapply(counts, function(c){
         return(data.frame(subject = c$saccver[1], score=sum(c$value)))
     })
     counts=do.call(rbind, counts)
     counts$fraction_total = counts$score / sum(counts$score)
     return(counts)
}
compiled = lapply(data, tally_hits)

# Join together
names(compiled) = sub(args$blast_results, pattern=".+1_(.+)\\..+txt$", repl="\\1")
for (n in names(compiled)){
    compiled[[n]]$sample = n
}
alldata = do.call(rbind, compiled)

# Add in taxonomy keys
alldata$species = taxonomy[alldata$subject, "species"]
alldata$is_organ = taxonomy[alldata$subject, "is_organ"]
alldata = subset(alldata, alldata$is_organ)

# Remove hits to species present in too few instances (default is no more than 1% in any sample)
totals = aggregate(alldata$fraction_total, by=list(alldata$sample, alldata$species), FUN=sum)
names(totals) = c('sample', 'species', 'total')
totals$fraction = totals$total / sum(totals$total)
to_keep = unique(totals$species[totals$total >= args$min_cutoff])
plotdata = subset(alldata, alldata$species %in% to_keep)

# Simplify and agregate for plotting
plotdata = aggregate(plotdata$fraction_total, by=list(plotdata$sample, plotdata$species), FUN=sum)
names(plotdata) = c('sample','species', 'fraction_total')

# Add in sample metadata
sample_data = sample_key[plotdata$sample, c("sample.type", "treatment")]
plotdata$sample_data = paste(sample_data$sample.type, sample_data$treatment, plotdata$sample, sep="|")
plotdata$sample_data = factor(plotdata$sample_data, levels=sort(target_breaks))

# Limit to just target samples
plotdata = subset(plotdata, plotdata$sample %in% target_samples)


###############
# Kraken Results
###############

cat("Processing Kraken results\n")
kraken_cols = c("cum_total", 'cum_hits', 'exact_hits', 'tax_level', 'tax_id', 'tax_name')
kraken = lapply(args$kraken_results, read.delim, header=F, col.names=kraken_cols, strip.white=TRUE)


# Join together
names(kraken) = sub(args$kraken_results, pattern=".+_report\\.(.+)\\.txt$", repl="\\1")
for (n in names(kraken)){
    kraken[[n]]$sample = n
}
allkraken = do.call(rbind, kraken)

# Put fraction total on same scale as BLAST (0-1 instead of 0-100%)
allkraken$fraction_total = allkraken$cum_total/100

# Filter to just the desired taxonomic level
plotkraken = subset(allkraken, allkraken$tax_level == args$taxon_level)

# Remove hits to species present in too few instances (default is no more than 1% in any sample)
to_keep = unique(plotkraken$tax_name[plotkraken$fraction_total >= args$min_cutoff])
plotkraken = subset(plotkraken, plotkraken$tax_name %in% to_keep)

# Add in sample metadata
sample_data = sample_key[plotkraken$sample, c("sample.type", "treatment")]
plotkraken$sample_data = paste(sample_data$sample.type, sample_data$treatment, plotkraken$sample, sep="|")
plotkraken$sample_data = factor(plotkraken$sample_data, levels=sort(target_breaks))

# Limit to just target samples
plotkraken = subset(plotkraken, plotkraken$sample %in% target_samples)


###############
# Output combined graphic
###############

cat("Outputting graphics to", args$outprefix, "\n")
mywidth=8
myheight=6

# Join datasets together
plotdata$set="blast"
plotdata$taxon=plotdata$species
plotkraken$set="kraken"
plotkraken$taxon=plotkraken$tax_name
columns = c("sample", "taxon", "fraction_total", "sample_data", "set")
masterplot = rbind(plotdata[,columns], plotkraken[,columns])

# Make graphic
output = ggplot(masterplot, mapping=aes(x=taxon, y=sample_data, fill=fraction_total)) +
    geom_tile() +
    facet_grid(~ set, scales="free_x") + 
    scale_fill_gradient(low="white", high="darkblue") +
    theme_classic() +
    theme( axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    theme(axis.text = element_text(size=6)) +
    labs(x="Taxon Matches", y="Sample", fill="Fraction Total") +
    scale_y_discrete(breaks=target_breaks, labels=target_labels, drop=FALSE)

ggsave(output, filename=paste(args$outprefix, "png", sep="."), width=mywidth, height=myheight, units='in', dpi=300)
ggsave(output, filename=paste(args$outprefix, "svg", sep="."), width=mywidth, height=myheight, units='in', dpi=300)

# Save R state to a file (used to rapidly load processed results to fiddle with figure formatting)
save.image(file=paste(args$outprefix, "Rdata", sep="."))
