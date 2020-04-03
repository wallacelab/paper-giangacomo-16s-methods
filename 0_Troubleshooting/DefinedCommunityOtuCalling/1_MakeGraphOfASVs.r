
# Distance file formed from command 'clustalo -i tmp.fa --distmat-out=tmp.dists --full --force'

library(igraph)
setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/0_Troubleshooting/DefinedCommunityOtuCalling/QiimeSandboxes/1_AssignOtus/')

# Read in data
dists=read.table('derep-seqs.dists.txt', row.names=1, skip=1)
graph = graph_from_adjacency_matrix(1-as.matrix(dists), weighted=T, mode='upper') # 1-dists because graphs formed by adjacency/similarities, not distances
# write_graph(graph, "gephi.gml", format="gml")

# Subset for quick checks
graph = graph_from_adjacency_matrix(1-as.matrix(dists)[1:1000,1:1000], weighted=T, mode='upper') # 1-dists because graphs formed by adjacency/similarities, not distances


# Prune edges below a certain value, and also all the self-edges
to_trim = (E(graph)$weight < 0.9) | (E(graph)$weight >= 1) # Everything below 0.9 or self-circles at 1
trimmed = delete_edges(graph, E(graph)[to_trim])


# Bring in taxonomy data for analysis
library(phyloseq)
phyloseq = readRDS('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/2_Analysis/2a_combined_data.phyloseq.RDS')

# Plot results
V(trimmed)$label = NA # Remove node labels
plot(trimmed)
