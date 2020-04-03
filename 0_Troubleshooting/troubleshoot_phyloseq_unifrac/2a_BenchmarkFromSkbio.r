#! /usr/bin/R

# # Install needed packages
# install.packages(c('devtools', 'picante', 'phyloseq', 'rbiom', 'ape'))
# devtools::install_github("cmmr/rbiom")

setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/0_Troubleshooting/troubleshoot_phyloseq_unifrac/')

# Load QIIME output files
counts = matrix(c(1, 3, 0, 1, 0, 0, 2, 0, 4, 4, 0, 0, 6, 2, 1, 0, 0, 1, 1, 1,
             5, 3, 5, 0, 0, 0, 0, 0, 3, 5), nrow=6, ncol=5, byrow=TRUE,
             dimnames=list(LETTERS[1:6], paste("OTU",1:5,sep="")))
counts=t(counts)    # Put into normal OTUs-as-rows format
tree=ape::read.tree(text="(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;")

benchmark_weighted=as.dist(read.csv('skbio_benchmark.weighted.csv', row.names=1))
benchmark_unweighted=as.dist(read.csv('skbio_benchmark.unweighted.csv', row.names=1))

# ########
# Generate UniFrac distances of everything
# ########

# # QIIME
# qiime_weighted = as.dist(read.delim("weighted_unifrac.tsv", row.names=1))
# qiime_unweighted = as.dist(read.delim("unweighted_unifrac.tsv", row.names=1))

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

weighted = list(benchmark=benchmark_weighted, phyloseq = phyloseq_weighted, rbiom = rbiom_weighted, gunifrac = as.dist(guni_weighted))
unweighted = list(benchmark=benchmark_unweighted, phyloseq = phyloseq_unweighted, rbiom = rbiom_unweighted, gunifrac = as.dist(guni_unweighted), picante=picante_unweighted)
compare_dists(weighted)
compare_dists(unweighted)
