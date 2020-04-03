#! /usr/bin/Rscript

# Benchmark UniFrac calculations in Phyloseq against other options

# ##########
# Install needed libraries
# ##########

install.packages(c('devtools', 'picante', 'phyloseq', 'rbiom', 'ape'))
devtools::install_github("cmmr/rbiom")


# ############
# Benchmark dataset (very simple)
# ############

# Make OTU table
u_counts = c(1, 0, 0, 4, 1, 2, 3, 0)
v_counts = c(0, 1, 1, 6, 0, 1, 0, 0)
counts=as.matrix(data.frame(u=u_counts, v=v_counts, row.names=paste("OTU", 1:8, sep="")))
tree = ape::read.tree(text="(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,(OTU5:0.5,((OTU6:0.33,OTU7:0.62):0.5,OTU8:0.5):0.5):0.5):1.25):0.0)root;")


# ########
# Generate UniFrac distances of everything
# ########

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

# Weighted values
rbiom_weighted
guni_weighted

# Unweighted values
rbiom_unweighted
guni_unweighted
picante_unweighted





# OLD DATA

# #############
# Generate artificial data of 1000 OTUs at total depth 10k across 100 samples
# #############

set.seed(1)
nsamples = 100
notus = 1000
depth = 10000

# Labels
sample_names = paste("sample",1:nsamples,sep="_")
otu_names = paste("otu",1:notus,sep="_")

# Make biom table
counts = matrix(0, nrow=notus, ncol=nsamples, dimnames=list(otu_names, sample_names))
for(i in 1:nsamples){
    reads = sample(1:nrow(counts), size=depth, replace=TRUE, prob=runif(notus))    # Gives a different sampling probability for each sample
    reads = table(reads)
    index = as.numeric(names(reads))
    counts[index,i] = reads
}

# Make tree
library(ape)
tree = rtree(n=notus, rooted=TRUE, tip.label=otu_names)


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

weighted = list(phyloseq = phyloseq_weighted, rbiom = rbiom_weighted, gunifrac = as.dist(guni_weighted))
unweighted = list(phyloseq = phyloseq_unweighted, rbiom = rbiom_unweighted, gunifrac = as.dist(guni_unweighted), picante=picante_unweighted)
compare_dists(weighted)
compare_dists(unweighted)
