#! /usr/bin/R

setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/0_Troubleshooting/troubleshoot_phyloseq_unifrac/qiime2_moving_pictures_rarefied/')

# Install needed packages
install.packages(c('devtools', 'picante', 'phyloseq', 'rbiom', 'ape'))
devtools::install_github("cmmr/rbiom")

# Load QIIME output files
counts=as.matrix(read.delim("feature-table.biom.txt", header=T, row.names=1, skip=1))
tree=ape::read.tree("tree.nwk")


# ########
# Generate UniFrac distances of everything
# ########

# QIIME
qiime_weighted = as.dist(read.delim("weighted_unifrac.tsv", row.names=1))
qiime_unweighted = as.dist(read.delim("unweighted_unifrac.tsv", row.names=1))

# Phyloseq
library(phyloseq)
phyloseq_data = phyloseq(otu_table(counts, taxa_are_rows=TRUE), tree)
phyloseq_weighted = phyloseq::UniFrac(phyloseq_data, weighted=TRUE)
phyloseq_unweighted = phyloseq::UniFrac(phyloseq_data, weighted=FALSE)

# rbiom
rbiom_weighted = rbiom::unifrac(counts, weighted=TRUE, tree=tree)
rbiom_unweighted = rbiom::unifrac(counts, weighted=FALSE, tree=tree)

# GUniFrac
gunifracs = GUniFrac::GUniFrac(t(counts), tree=tree, alpha=1)$unifracs
guni_weighted = gunifracs[, , "d_1"]
guni_unweighted = gunifracs[, , "d_UW"]

# picante (unweighted unifrac only, and definitely the slowest)
picante_unweighted = picante::unifrac(t(counts), tree)


# ###########
# Compare methods
# ###########

compare_dists = function(mydists){
    comparison = matrix(NA, nrow=length(mydists), ncol=length(mydists), dimnames=list(names(mydists), names(mydists)))
    for(i in 1:length(mydists)){
        for(j in i:length(mydists)){
            comparison[i,j] = cor(mydists[[i]], mydists[[j]])
        }
    }
    return(comparison)
}

weighted = list(qiime=qiime_weighted, phyloseq = phyloseq_weighted, rbiom = rbiom_weighted, gunifrac = as.dist(guni_weighted))
unweighted = list(qiime=qiime_unweighted, phyloseq = phyloseq_unweighted, rbiom = rbiom_unweighted, gunifrac = as.dist(guni_unweighted), picante=picante_unweighted)
compare_dists(weighted)
compare_dists(unweighted)
