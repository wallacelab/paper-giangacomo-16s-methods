#! /usr/bin/R

# # Install needed packages
# install.packages(c('devtools', 'picante', 'phyloseq', 'rbiom', 'ape'))
# devtools::install_github("cmmr/rbiom")

library(phyloseq)

# Sample Global Patterns subset data from the phyloseq unit tests
treeFile = system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")
GP500File = system.file("extdata", "GP_otu_table_rand_short.txt.gz", package = "phyloseq")
GP500 = import_qiime(GP500File, treefilename = treeFile)

# Now import the results with read.table()
gp500_uuf = read.csv(system.file("extdata", "gp500-uuf.csv", package = "phyloseq"), header = FALSE, fill = TRUE)
gp500_wuf = read.csv(system.file("extdata", "gp500-wuf.csv", package = "phyloseq"), header = FALSE, fill = TRUE)
gp500_wufu = read.csv(system.file("extdata", "gp500-wufu.csv", package = "phyloseq"), header = FALSE, fill = TRUE)
# Add the sample names
colnames(gp500_uuf) <- rownames(gp500_uuf) <- colnames(gp500_wuf) <- rownames(gp500_wuf) <- colnames(gp500_wufu) <- rownames(gp500_wufu) <- sample_names(GP500)
# Coerce to Distance Matrices for comparison `"dist"` class
gp500_wufu <- as.dist(gp500_wufu)
gp500_wuf <- as.dist(gp500_wuf)
gp500_uuf <- as.dist(gp500_uuf)


# ########
# Generate UniFrac distances of everything
# ########

counts = otu_table(GP500)
tree = phy_tree(GP500)


# Phyloseq
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

weighted = list(gp500_wufu=gp500_wufu, gp500_wuf=gp500_wuf, phyloseq = phyloseq_weighted, rbiom = rbiom_weighted, gunifrac = as.dist(guni_weighted))
unweighted = list(gp500_uuf=gp500_uuf, phyloseq = phyloseq_unweighted, rbiom = rbiom_unweighted, gunifrac = as.dist(guni_unweighted), picante=picante_unweighted)
compare_dists(weighted)
compare_dists(unweighted)
