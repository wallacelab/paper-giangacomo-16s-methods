# Compare QIIME- and phyloseq-filtered OTU tables to see how compare

library(phyloseq)

load("C:/Users/jgwall/Desktop/2019 07 Cecelia Final Data/2_Analysis/2b_filtered_data.phyloseq.Rdata")
phylo=as.data.frame(otu_table(mydata))
qiime=read.delim("C:/Users/jgwall/Desktop/000_UpdateNAPB_DC/CeceliaMicrobiome/1_AssignOTUs_additions/feature-table.biom.txt", 
                 skip=1, check.names=F, row.names=1)


# Compare samples
all_samples = union(names(phylo), names(qiime))
shared_samples = intersect(names(phylo), names(qiime))
different_samples = setdiff(all_samples, shared_samples)
cat("QIIME and Phyloseq filtering have",length(different_samples),"different samples:",different_samples,"\n")
identical(names(phylo), names(qiime))

# Compare OTUs
all_otus = union(rownames(phylo), rownames(qiime))
shared_otus = intersect(rownames(phylo), rownames(qiime))
different_otus = setdiff(all_otus, shared_otus)
cat("QIIME and Phyloseq filtering have",length(different_otus),"different otus:",different_otus,"\n")
identical(rownames(phylo), rownames(qiime))

# Reorder QIIME to match phyloseq
qiime2 = qiime[rownames(phylo), names(phylo)]
identical(names(phylo), names(qiime2))
identical(rownames(phylo), rownames(qiime2))

# Compare sample and OTU depth
identical(colSums(phylo), colSums(qiime2))
identical(rowSums(phylo), rowSums(qiime2))


# Okay, looks like the filtering worked exactly the same, so I can integrate QIIME filtering 
# into the pipeline and just pick it up from there.