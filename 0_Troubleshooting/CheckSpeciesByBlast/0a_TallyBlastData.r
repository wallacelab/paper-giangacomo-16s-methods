#! /usr/bin/Rscript

# Tally up BLAST results to get an idea of species identity

# Libraries
library(argparse)
library(ggplot2)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs="*", help="List of BLAST output files")
parser$add_argument("-o", "--outfile", help="Consolidated output file")
parser$add_argument("-m", "--min-total-score", default=1, type='double', help="do not display sequences with fewer than this many total hits in the heatmap")
parser$add_argument("--outgraphics", help="PNG file to write graphics to")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2020 03 Consolidated Pipeline/0_Troubleshooting/CheckSpeciesByBlast/')
# args=parser$parse_args(c("-i","1_M100.blast.txt", "1_M101.blast.txt", "1_M102.blast.txt", "-o", "99_tmp.txt", "--outgraphics", "99_tmp.png"))

# Load data
cat("Compiling BLAST results from", length(args$infiles), "input files\n")
blast_columns=c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data = lapply(args$infiles, read.delim, col.names=blast_columns)

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

    return(counts)
}
compiled = lapply(data, tally_hits)

# Join together
names(compiled) = sub(args$infiles, pattern="1_(.+)\\..+txt$", repl="\\1")
for (n in names(compiled)){
    compiled[[n]]$sample = n
}
alldata = do.call(rbind, compiled)

# Remove too low-scoring hits
totals = aggregate(alldata$score, by=list(alldata$subject), FUN=sum)
names(totals) = c('subject', 'total')
to_keep = totals$subject[totals$total >= args$min_total_score]
alldata = subset(alldata, alldata$subject %in% to_keep)
cat("\t", length(to_keep), "reference seqs out of", nrow(totals), "have min total score of",args$min_total_score,"\n")

# Plot
cat("Graphing results to", args$outgraphics, "\n")
myplot = ggplot(alldata, mapping=aes(x=subject, y=sample, fill=score)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="blue")
ggsave(myplot, file=args$outgraphics)
