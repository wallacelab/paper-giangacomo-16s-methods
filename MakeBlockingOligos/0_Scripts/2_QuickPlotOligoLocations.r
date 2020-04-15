#! /usr/bin/Rscript
 
library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Compiled BLAST hits")
parser$add_argument("-o", "--outfile")
args=parser$parse_args()
# setwd("/home/jgwall/Projects/MicrobiomeMethodsDevelopment/2017 09 Blocking oligo design/2_BlockingOligos/")
# args=parser$parse_args(c("-i","2b_combined_blast.txt", "-o","99_tmp.png"))

cat("Plotting crude locations of oligos and sequences\n")
# Load data
data=read.delim(args$infile)
seqs = split(data, data$sseqid)

# Plotter function
plot_region=function(myseq){
  title=unique(myseq$sseqid[1])
  yvals=1:nrow(myseq)
  xlim=range(c(myseq$sstart, myseq$send))
  ylim=range(yvals)
  
  
  # Plot blank
  plot(NA, NA, xlim=xlim, ylim=ylim, main=title, xlab="", ylab="")
  
  # Add lines
  segments(x0=myseq$sstart, x1=myseq$send, y0=yvals, col=myseq$qseqid, lwd=6)
  
  # Add text
  text(x=max(xlim), y=yvals, labels=as.character(myseq$qseqid), pos=4)
}

# Actually plot and save
nrow=length(seqs)
png(args$outfile, width=800, height=300 * nrow)
  par(mfrow=c(nrow, 1), mar=c(4,4,2,10), cex=1.2, xpd=NA)
  tmp=lapply(seqs, plot_region)
dev.off()