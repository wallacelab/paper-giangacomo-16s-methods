#! /usr/bin/Rscript

# Convert OTU counts (transformed or otherwise) into BLUPs

library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-b", "--blast", nargs="*", help="List of BLAST results matching ASVs to Silva OTUs")
parser$add_argument("-t", "--tree", help="NEWICK-formatted phylogenetic tree of sequences")
parser$add_argument("-x", "--taxonomy", help="Taxonomy key for sequences in the phylogenetic tree")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/1_AssignOtus/')
# args=parser$parse_args(c("-b","PNAs/1g_rep-seqs.blast_results.txt", "Universal/1g_rep-seqs.blast_results.txt", "-t", "/home/jgwall/Projects/0_RawData/Silva_132_release/silva_99_otus.tre", 
#   "-x", "/home/jgwall/Projects/0_RawData/Silva_132_release/majority_taxonomy_7_levels.99.txt", "-o",'99_tmp'))

# Libraries
library(ape)

# Function to add a new tip to a phylogenetic tree; extracted from the phytree package at https://github.com/liamrevell/phytools/blob/master/R/utilities.R [7 Feb 2020, commit 68ae5cf]
# written by Liam J. Revell 2012, 2013, 2014, 2015
bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL,position=0,interactive=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	use.edge.length<-if(is.null(tree$edge.length)) FALSE else TRUE
	if(use.edge.length==FALSE) tree<-compute.brlen(tree)
	if(interactive==TRUE){
		plotTree(tree,...)
		cat(paste("Click where you would like to bind the tip \"",tip.label,"\"\n",sep=""))
		flush.console()
		obj<-get.treepos(message=FALSE)
		where<-obj$where
		position<-obj$pos
	} else if(is.null(where)) where<-Ntip(tree)+1
	if(where<=Ntip(tree)&&position==0){
		pp<-1e-12
		if(tree$edge.length[which(tree$edge[,2]==where)]<=1e-12){
			tree$edge.length[which(tree$edge[,2]==where)]<-2e-12
			ff<-TRUE
		} else ff<-FALSE
	} else pp<-position
	if(is.null(edge.length)&&is.ultrametric(tree)){
		H<-nodeHeights(tree)
		if(where==(Ntip(tree)+1)) edge.length<-max(H)
		else edge.length<-max(H)-H[tree$edge[,2]==where,2]+position
	}
	tip<-list(edge=matrix(c(2,1),1,2),
		tip.label=tip.label,
		edge.length=edge.length,
		Nnode=1)
		class(tip)<-"phylo"
	obj<-bind.tree(tree,tip,where=where,position=pp)
	if(where<=Ntip(tree)&&position==0){
		nn<-obj$edge[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label)),1]
		obj$edge.length[which(obj$edge[,2]==nn)]<-obj$edge.length[which(obj$edge[,2]==nn)]+1e-12
		obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label))]<-0
		obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tree$tip.label[where]))]<-0
	}
	root.time<-if(!is.null(obj$root.time)) obj$root.time else NULL
	#obj<-untangle(obj,"read.tree")    # Takes a huge amount of time; tree will be written out later, so just removed
	if(!is.null(root.time)) obj$root.time<-root.time
	if(interactive) plotTree(obj,...)
	if(!use.edge.length) obj$edge.length<-NULL
	obj
}

# Load data
cat("Loading data from", length(args$blast), "input files\n")
blast = lapply(args$blast, read.delim, col.names=c("query", "subject", "percent_identity","length", "mismatch", "gapopen", "quert_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"))
taxonomy = read.delim(args$taxonomy, row.names=1, header=FALSE, col.names=c("OTU", "tax_string"))
tree = read.tree(args$tree)

# Check for same ASVs across datasets (can happen with same primers)
cat("Checking for common ASVs across datasets\n")
unique_queries = lapply(blast, function(b){unique(b$query)})
query_counts = table(unlist(unique_queries))
duplicated = data.frame(query=names(query_counts), counts=as.numeric(query_counts))
duplicated = subset(duplicated, duplicated$counts > 1)
cat("\t",nrow(duplicated),"ASV identifiers are duplicated; may want to manually check to make sure are the same\n")
write.table(duplicated, file=paste(args$outprefix, ".duplicated_asvs_from_blast.txt", sep=""), quote=F, row.names=F, col.names=T, sep='\t')
blast=do.call(rbind, blast)

# Reduce to just the best BLAST hits for each sequence (determined by bitscore)
cat("Finding best BLAST hits for",length(unique(blast$query)),"ASVs\n")
hits = split(blast, blast$query)
best_hits = lapply(hits, function(myhits){
    myhits = subset(myhits, myhits$bitscore==max(myhits$bitscore))# Only take maximum bitscore
})

# Compare against taxonomy and make sure everything is the same at the genus level
cat("Sanity-checking with taxonomy to be sure genera are mostly the same\n")
taxonomy$genus = sub(taxonomy$tax_string, pattern=".+D_5__(.+);.+", repl="\\1")
add_genus = lapply(best_hits, function(myhits){
    tax_index = match(myhits$subject, rownames(taxonomy))
    myhits$genus = taxonomy$genus[tax_index]
    myhits$tax_string = taxonomy$tax_string[tax_index]
    return(myhits)
})
check_genus = sapply(add_genus, function(x){
    return(length(unique(x$genus))==1) # Returns TRUE if only 1 genus in dataset, FALSE otherwise    
})
cat("\t",sum(check_genus),"ASVs have perfect matches at genus level;", sum(!check_genus),"have mismatches.\n")
cat("\t\t(Mismatches are usually due to ambiguous taxa; check the raw output to confirm\n")

# Reduce to the first BLAST hit (under assumption that anything with the same, best bitscore is nearly equivalent)
cat("Taking the first BLAST hit as a representative OTU\n")
final_hits = lapply(best_hits, function(x){
    x[1,]
})
final_hits=do.call(rbind, final_hits)

# Add ASVs to the taxonomic tree as children of their matched OTU
cat("Adding ASVs to phylogenetic tree as children of the selected OTUs (Note: this will take a while)\n")
# # # # # bind.tip<-function(mytree,tip.label,edge.length=NULL,where=NULL){ # Function from Liam Revell at http://blog.phytools.org/2012/11/adding-single-tip-to-tree.html
# # # # #     if(is.null(where)) where<-length(mytree$tip)+1
# # # # #     tip<-list(edge=matrix(c(2,1), nrow=1, ncol=2),
# # # # #               tip.label=tip.label,
# # # # #               edge.length=edge.length,
# # # # #               Nnode=1)
# # # # #     class(tip)<-"phylo"
# # # # #     obj<-bind.tree(mytree,tip,where=where, position=1e-12)  # Need to have a non-zero position or just overwrites existing node
# # # # #     return(obj)
# # # # # }
newtree=tree # Make a copy just in case; helps with debugging 
for(i in 1:nrow(final_hits)){
    if(i %% 100 == 0){
        cat("\tProcessed",i,"ASVs into the tree\n")
    }
    target_tip = which(newtree$tip.label==final_hits$subject[i])
    newtree = bind.tip(newtree, tip.label=final_hits$query[i], edge.length=0, where=target_tip)
}

# Write out results
cat("Writing out results\n")
write.tree(newtree, file=paste(args$outprefix, ".tre", sep=""))
genus_table=do.call(rbind, add_genus)[,c("query", "subject", "bitscore", "genus", "tax_string")]
write.table(genus_table, file=paste(args$outprefix, ".genera.txt", sep=""), quote=F, row.names=F, col.names=T, sep='\t')

