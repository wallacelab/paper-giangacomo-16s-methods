#! /bin/bash

# Perform quality control and filtering of 16s reads

# Arguments
datadir=$1      # Directory of raw sequence reads
workdir=$2      # Working directory to output results
keyfile=$3      # QIIME-formatted keyfile
forward_primers=$4  # Space-separated list of forward primers to remove from reads
reverse_primers=$5  # Space-separated list of reverse primers to remove from reads
rev_trim=$6         # How many bases to trim off the 3' end of reverse reads
min_length=$7    # Minimum length of sequence to keep after joining
references="$8"    # Classifer artifact to use for adding taxonomy to reads
otu_identity=$9  # Fraction identity for OTU cluster (so, 0.99 = 99%)
num_cores=${10}  # number of CPU cores to use in cutadapt
#silva_taxonomy=${11} # Taxonomy file for use with reference OTUs
#defined_community=${12}

# Set up subdirectories
trimdir=$workdir/1_trimmed_seqs
filtdir=$workdir/1a_filtered_seqs

if [ ! -e $trimdir ]; then mkdir $trimdir; fi
if [ ! -e $filtdir ]; then mkdir $filtdir; fi


# #######################
# # Step 1 - Basic sequence cleaning
# #######################
# 
# # Trim the last X bases from the reverse reads, since this seems to be the source of most read-joining errors
# cp $datadir/*_R1_*.fastq.gz $trimdir        # Just copy over the forward reads
# for rev in $datadir/*_R2_*.fastq.gz; do     # Trim reverse reads 
#     filename=`basename $rev`
#     cutadapt --cut -$rev_trim --cores $num_cores -o "$trimdir/$filename" "$rev"
# done
# 
# # Filter out reads with no data
# for fwd in $trimdir/*_R1_*.fastq.gz; do     # Trim reverse reads 
#     rev=${fwd/_R1_/_R2_}
#     fwd_out=`basename $fwd`
#     rev_out=`basename $rev`
#     cutadapt --cores $num_cores -m 1 -o "$filtdir/$fwd_out" -p "$filtdir/$rev_out" "$fwd" "$rev"    # "-m 1" means "at least 1 bp long", so filter out failed sequences
# done
# 
# 
# # # Import sequences into QIIME (assumes Casava 1.8 paired-end demultiplexed data, which is what the original data was in; see https://docs.qiime2.org/2019.7/tutorials/importing/#importing-seqs )
# echo "Importing FASTQ sequences into QIIME artifact"
# qiime tools import \
#   --type 'SampleData[PairedEndSequencesWithQuality]' \
#   --input-path $filtdir \
#   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#   --output-path $workdir/1_demux-paired-end.qza
#   
# 
# # # Use cutadapt to remove the primer sequences from forward and reverse reads
# echo "Removing primers with cutadapt using $num_cores processing cores"
# qiime cutadapt trim-paired \
#     --i-demultiplexed-sequences $workdir/1_demux-paired-end.qza \
#     --o-trimmed-sequences $workdir/1a_cutadapt-paired-end.qza \
#     --p-front-f $forward_primers \
#     --p-front-r $reverse_primers \
#     --p-cores $num_cores
# 
#     
# # Join paired ends with vsearch
# echo "Joining paired reads"
# qiime vsearch join-pairs \
#   --i-demultiplexed-seqs $workdir/1a_cutadapt-paired-end.qza \
#   --o-joined-sequences $workdir/1b_joined-paired-end.qza \
#   --p-minmergelen $min_length \
#   --verbose 
#   
# # Summary/visualization of joined pairs
# echo "Summarizing joined pairs"
# qiime demux summarize \
#   --i-data $workdir/1b_joined-paired-end.qza \
#   --o-visualization $workdir/1b_joined-paired-end.summary.qzv
#   
# # Quality filter joined pairs
# echo "Quality filtering joined pairs"
# qiime quality-filter q-score-joined \
#   --i-demux $workdir/1b_joined-paired-end.qza \
#   --o-filtered-sequences $workdir/1c_quality-filtered.qza \
#   --o-filter-stats $workdir/1c_quality-filtered.stats.qza


#######################
# Step 2 - Clustering and OTU calling (not doing individual sequence variants with deblur or the like b/c of different amplicon lengths for the different primers)
#######################
 
# # Dereplicate sequences (=collapse ones with identical sequences)
# echo "Dereplicating reads with vsearch"
# qiime vsearch dereplicate-sequences \
#   --i-sequences $workdir/1c_quality-filtered.qza  \
#   --o-dereplicated-table $workdir/1d_derep-table.qza \
#   --o-dereplicated-sequences $workdir/1d_derep-seqs.qza


if [ $references == "NONE" ]; then
    echo "### De Novo OTU clustering ###"
    qiime vsearch cluster-features-de-novo \
      --i-table $workdir/1d_derep-table.qza \
      --i-sequences $workdir/1d_derep-seqs.qza \
      --p-perc-identity $otu_identity \
      --p-threads $num_cores \
      --o-clustered-table $workdir/1e_clustered-table.$otu_identity.qza \
      --o-clustered-sequences $workdir/1e_clustered-seqs.$otu_identity.qza \

else
    echo "### Reference-based OTU clustering ###"
    # Import the reference sequences to qiime
    echo "Importing reference sequences"
    qiime tools import \
      --input-path "$references" \
      --output-path $workdir/1e_reference-seqs.qza  \
      --type 'FeatureData[Sequence]'
    
    
    # Closed-reference OTU clustering with vsearch
    echo "Clustering OTUs with closed-reference pipeline using vsearch at $otu_identity identity"
    qiime vsearch cluster-features-closed-reference \
      --i-table $workdir/1d_derep-table.qza \
      --i-sequences $workdir/1d_derep-seqs.qza \
      --i-reference-sequences $workdir/1e_reference-seqs.qza  \
      --p-perc-identity $otu_identity \
      --p-threads $num_cores \
      --o-clustered-table $workdir/1e_clustered-table.$otu_identity.qza \
      --o-clustered-sequences $workdir/1e_clustered-seqs.$otu_identity.qza \
      --o-unmatched-sequences $workdir/1e_unmatched-seqs.$otu_identity.qza  
fi
    





# Summarize the resulting feature table  
echo "Summarizing feature table"
qiime feature-table summarize \
  --i-table $workdir/1e_clustered-table.$otu_identity.qza \
  --o-visualization $workdir/1e_clustered-table.$otu_identity.visualization.qzv \
  --m-sample-metadata-file $keyfile
   
  
###################
# Step 3 - Extract data for other programs to access
###################

qiime tools export --input-path $workdir/1d_derep-table.qza --output-path $workdir
qiime tools export --input-path $workdir/1d_derep-seqs.qza --output-path $workdir
mv $workdir/feature-table.biom $workdir/derep-seqs.biom
mv $workdir/dna-sequences.fasta $workdir/derep-seqs.fasta

# Extract needed data to put everything into Phyloseq (=feature table and phylogenetic tree)
qiime tools export --input-path $workdir/1e_clustered-table.$otu_identity.qza --output-path $workdir
qiime tools export --input-path $workdir/1e_clustered-seqs.$otu_identity.qza --output-path $workdir 
mv $workdir/feature-table.biom $workdir/clustered-seqs.$otu_identity.biom
mv $workdir/dna-sequences.fasta $workdir/clustered-seqs.$otu_identity.fasta

# Convert biom table to text (less likely to have import issues into phyloseq)
biom convert -i $workdir/clustered-seqs.$otu_identity.biom -o $workdir/clustered-seqs.$otu_identity.biom.txt --to-tsv



###################
# DEPRECATED Step 4 - Alignment network - 
###################

# echo "Calculating pairwise distances with ClustalOmega"
# clustalo -i $workdir/clustered-seqs.fasta --distmat-out=$workdir/clustered-seqs.dists.txt --full --force --threads $num_cores > /dev/null
# clustalo -i $workdir/derep-seqs.fasta --distmat-out=$workdir/derep-seqs.dists.txt --full --force --threads $num_cores > /dev/null

# # Making a Gephi graph in R
# library(igraph)
# dists=read.table('derep-seqs.short.dists.txt', row.names=1, skip=1)
# graph = graph_from_adjacency_matrix(dists, weighted=T, mode='upper')
# write_graph(graph, "gephi.gml", format="gml")

# TODO: Still trying to figure out how to separate out clusters well here. Need to cut/remove edges that are too weak?
