#! /usr/bin/Rscript

# Convert QIIME derived files into a single phyloseq object

library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-b", "--biom", nargs="*", help="List of text-formatted BIOM file of microbiome observations from QIIME")
parser$add_argument("-t", "--tree", help="NEWICK-formatted phylogenetic tree of sequences (should have ASVs added)")
parser$add_argument("-x", "--taxonomy", nargs="*", help="List of taxonomy keys for sequences in the BIOM file")
parser$add_argument("-k", "--keyfile", help="QIIME-formatted keyfile of sample data")
parser$add_argument("-o", "--outprefix", help="File name prefix for all output files")
parser$add_argument("-s", "--samples-to-remove", nargs="*", help="Space-separated list of bad samples to be filtered out")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 10 Mohsen Final Data//1_AssignOtus/')
# args=parser$parse_args(c("-b","Plants/deblur-seqs.biom.txt", "Soils/deblur-seqs.biom.txt", "-t", "../2_Analysis/2_combine_phylogenies.tre", "-x", "Plants/taxonomy_formatted.tsv", "Soils/taxonomy_formatted.tsv",
#  "-k", "../16s_extractions_keyfile.tsv", "-o",'99_tmp.RDS', "-s", "M135"))

# Libraries
library(phyloseq)

# Load data
cat("Loading data\n")
biom = lapply(args$biom, read.delim, skip=1, check.names=F, row.names=1)
taxonomy = lapply(args$taxonomy, read.delim, row.names=1, header=FALSE)
key = read.delim(args$keyfile, row.names=1)[-1,]    # Remove the second line, which is just comments on the qiime data types
tree = read_tree(args$tree)

cat("Reformatting data\n")

# Reformat taxonomy from SILVA to phyloseq-compatible
taxlevels = list()
for(i in 1:length(taxonomy)){
    taxlevels[[i]] = strsplit(taxonomy[[i]]$V2, split=";")
    taxlevels[[i]] = lapply(taxlevels[[i]], function(x){
        x = sub(x, pattern="^D_.__", repl="")
        if(length(x) < 7){  # Add in any missing taxonomy levels as "unknown"
            x[(length(x)+1):7] = "unknown"
        }
        names(x) = c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species")
        return(x)
    })
    taxlevels[[i]] = do.call(rbind, taxlevels[[i]])
    rownames(taxlevels[[i]]) = rownames(taxonomy[[i]])
}

# Merge individual taxonomies into a single one
asvs = lapply(taxlevels, function(x){
    x=as.data.frame(x)
    x$asv=row.names(x)
    return(x)
})
asvs = as.data.frame(do.call(rbind, asvs))
asvs.unique = unique(asvs)
taxonomy = tax_table(as.matrix(asvs.unique))

# Identify if any ASVs still present in more than one copy (happens if entire taxonomy string is not identical)
tmp_taxonomy = as.data.frame(taxonomy)
asv_counts = table(tmp_taxonomy$asv)
duplicated = asv_counts[asv_counts >1]
if(length(duplicated)>0){
    dup_out = paste(args$outprefix, ".duplicated_asvs.txt", sep="")
    cat("\tWARNING:",length(duplicated),"ASVs are present in >1 copy in taxonomy; only the first will be used.\n")
    cat("\t\tThis usually happens if the taxonomy line is not 100% identical; check",dup_out,"for specific list.\n")
    dup_taxa = subset(tmp_taxonomy, tmp_taxonomy$asv %in% names(duplicated))
    write.table(dup_taxa, file=dup_out, sep="\t", quote=F, row.names=T, col.names=T)
}


# Convert biom tables and merge together
biom = lapply(biom, function(b){
    otu_table(as.matrix(b), taxa_are_rows=TRUE)
})
# Merge biom tables
biom.master = biom[[1]]
for(i in 2:length(biom)){
    biom.master = merge_phyloseq(biom.master, biom[[i]])    # Can only merge two at a time
}

# Convert sample key
key = sample_data(type.convert(key))    # type.convert changes things to factors, integer, etc. instead of the default "character" 

# Merge into a single phyloseq object
mydata = phyloseq(biom.master, taxonomy, tree, key) 

# Confirm all samples and ASVs made the merge
all_samples = unique(unlist(lapply(biom, colnames)))
all_asvs = unique(unlist(lapply(biom, rownames)))
missing_samples = setdiff(all_samples, sample_names(mydata))
missing_asvs = setdiff(all_asvs, taxa_names(mydata))
if(length(missing_samples) > 0){
    warning(length(missing_samples),"samples were excluded from the final BIOM file.\n")
}
if(length(missing_asvs) > 0){
    warning(length(missing_asvs)," ASVs were excluded from the final BIOM file.\n")
}
write(missing_samples, file=paste(args$outprefix, ".missing_samples.txt", sep=""))
write(missing_asvs, file=paste(args$outprefix, ".missing_asvs.txt", sep=""))


# Remove bad samples
cat("Removing", length(args$samples_to_remove), "bad samples\n")
mysamples = sample_names(mydata)
tokeep = mysamples[! mysamples %in% args$samples_to_remove]
mydata = prune_samples(mydata, samples=tokeep)
removed = setdiff(mysamples, sample_names(mydata))
cat("\tFound",length(removed),"samples to remove:",sort(removed),"\n")


# Write out the phyloseq object (so don't have to do this again)
out_rds = paste(args$outprefix, ".RDS", sep="")
cat("Writing phyloseq object to",out_rds, "\n")
cat("\tFinal phyloseq object has",nsamples(mydata),"samples and",ntaxa(mydata),"OTUs\n")
saveRDS(mydata, file=out_rds)


