# Look at defined community
library(phyloseq)
library(psadd)
options(stringsAsFactors=F)

# Input
setwd("/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/0_Troubleshooting/DefinedCommunityOtuCalling/")
infile="../../2_Analysis/2b_filtered_data.phyloseq.RDS"

# Load data
cat("Loading data for fraction of OTUs captured\n")
mydata = readRDS(infile)
defined = prune_samples(mydata, samples=sample_data(mydata)$sample.type == "defined-community")
defined = prune_taxa(defined, taxa=taxa_sums(defined)>0)

# Plot 
otus = as.data.frame(otu_table(defined))
counts = rowSums(otus)
sort(counts, decreasing=T)[1:10]
plot(sort(counts, decreasing=T)[1:50])

# Krona plot of actual data
plot_krona(defined, output="defined_comm.filtered", variable="treatment")

# ############
# Theoretical defined community
# ############

# Build up theoretical values for each
zymo_species = c("P_aeruginosa", "E_coli","S_enterica", "L_fermentum","E_faecalis","S_aureus","L_monocytogenes","B_subtilis")
zymo_counts = c(4.2,10.1,10.4,18.4,9.9,15.5,14.1,17.4)
zymo_table = data.frame(theoretical=zymo_counts, row.names=zymo_species)  # OTU table
zymo_data=data.frame(treatment="theoretical", row.names="theoretical")  # Sample data
zymo_taxonomy = read.delim('zymo_taxonomy.tsv', row.names=1)

# Convert to necessary reqs
zymo_counts = otu_table(as.matrix(zymo_table), taxa_are_rows=TRUE)
zymo_data = sample_data(zymo_data)
zymo_taxonomy = tax_table(as.matrix(zymo_taxonomy))
zymo = phyloseq(zymo_counts, zymo_data, zymo_taxonomy)

# krona plot
plot_krona(zymo, output="defined_comm.theoretical", variable="treatment")


# ###############
# Try to figure out what happened
# ###############

mytaxons = as.data.frame(tax_table(defined))
# Looking for specific abundant taxa that shouldn't be there

acineto = prune_taxa(defined, taxa=mytaxons$Species=="Acinetobacter baumannii") 
# Most abundant is LN611354.1.1200. Taking the seq, it's >97% similar to a P. aeruginosa strain

bifido = prune_taxa(defined, taxa=mytaxons$Genus=="Bifidobacterium") 
# Most abudnant are JN020355.1.1242 and EU767376.1.1347


# ###############
# Turns out is actually ATCC community. Figure out which OTUs correspond to those
# ###############

universal = prune_samples(defined, samples = sample_data(defined)$treatment=="Universal")
univ_otus = otu_table(universal)
univ_taxa = as.data.frame(tax_table(universal))
atcc = data.frame(species=univ_taxa$Species, genus=univ_taxa$Genus, counts=rowSums(univ_otus))
atcc=atcc[order(atcc$counts, decreasing=T),]
head(atcc, 20)
plot(atcc$counts[1:100])
