#! /usr/bin/env bash 

# This is the master script to analyze the data for Giangacomo et al 2019, on comparing methods for extracting and amplifying microbiomes from plant-associated samples
# This script was run in Linux Mint 18.3 with KDE 64-bit, although the majority of actual computation is meant to take place within a Conda environment

############
# SETUP
############

# List of subdirectories
datadir=0_data                      # Folder for raw sequencing files and keyfile to go in

# Core files and options
keyfile=defined_community_test_keyfile.tsv    # QIIME-formatted keyfile of samples and metadata
SILVA=/home/jgwall/Projects/0_RawData/Silva_132_release # Base directory to SILVA data
silva_seqs=$SILVA/silva_132_99_16S.fna # Path to QIIME-compatible SILVA 16s RNA library at 99% identity
silva_taxonomy=$SILVA/majority_taxonomy_7_levels.99.txt # Path to the corresponding taxonomy table for the above
silva_tree=$SILVA/silva_99_otus.tre # SILVA phylogenetic tree
ncores=7    # Number of processor cores to use
conda_qiime2=qiime2-2019-7      # Conda environment with QIIME 2019.7
conda_phyloseq=phyloseq-1.28.0  # Conda environment with Phyloseq
conda_dada2=qiime2-2019.7-dada2

# All primer sequences used to generate this data (Note: some include Illumina linkers) TODO: Make sure cutadapt handles these properly
forward_primers="AACMGGATTAGATACCCKG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAACMGGATTAGATACCCKG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGACCMGGATTAGATACCCKG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGAGGCAGCAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGYCAGCMGCCGCGGTAA"
reverse_primers="ACGTCATCCCCACCTTCC GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGACGTCATCCCCACCTTCC GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGGGTTGCGCTCGTTG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACHVGGGTWTCTAAT GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACNVGGGTWTCTAAT"

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

##############
# ANALYSIS PIPELINES
##############



# #####################
# # Default: Process raw sequence data down to BIOM feature table
# #####################
# 
# # # Clustering parameter to group OTUs; 0.99 = 99% identity, 0.97 = 97% identity, etc. Used in multiple steps
# workdir=1_AssignOtus               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi
# 
# otu_identity=0.99
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=$silva_seqs  # Reference sequences to use for OTU clustering
# defined_community=0_atcc_species.99.txt    # Defined community we're aiming for
# # conda activate $conda_qiime2   # Activate the qiime2 environment
# # bash 1_ProcessSequenceToBiom.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores
# # Look at how well coverage worked
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/clustered-seqs.$otu_identity.biom.txt --taxonomy $silva_taxonomy --targets $defined_community -o $workdir/1f_defined_community_metrics  -k $keyfile


#####################
# Try with the Zymobiomics standards; skipped first several steps and copied over 1d_* from the above to start there
#####################

# UPDATE: Eventually realized we hadn't used the Zymo standard
# workdir=1a_AssignOtus_ZymoTheoretical               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi
# otu_identity=0.90
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references="/home/jgwall/Projects/Microbiomes/MicrobiomeMethodsDevelopment/CompareSampleExtractionAndAmplification_Mohsen_Cecelia/2019 07 Cecelia Final Data/0_Troubleshooting/DefinedCommunityOtuCalling/Zymobiomics rRNA seqs/zymo_16s.fa"
# conda activate $conda_qiime2   # Activate the qiime2 environment
# bash 1_ProcessSequenceToBiom.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length "$references" $otu_identity $ncores


#####################
# Alternative: Use 97% cutoff instead
#####################

# workdir=1_AssignOtus_97percent               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi
# otu_identity=0.97
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=$SILVA/silva_132_97_16S.fna # Path to QIIME-compatible SILVA 16s RNA library at 99% identity
# silva_taxonomy=$SILVA/majority_taxonomy_7_levels.97.txt # Path to the corresponding taxonomy table for the above
# defined_community=0_atcc_species.97.txt    # Defined community we're aiming for
# conda activate $conda_qiime2   # Activate the qiime2 environment
# bash 1_ProcessSequenceToBiom.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores

# # # Look at how well coverage worked
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/clustered-seqs.$otu_identity.biom.txt --taxonomy $silva_taxonomy --targets $defined_community -o $workdir/1f_defined_community_metrics -k $keyfile



#####################
# Alternative: Do de novo clustering instead. NOTE: Just copied over 1d_derep-seqs.qza from the default folder since is the same up to that point and ran from there
#####################

# workdir=1_AssignOtus_denovo               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi
# otu_identity=0.99
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=NONE # Path to QIIME-compatible SILVA 16s RNA library at 99% identity, or NONE to do de novo clustering
# silva_taxonomy=NONE # Path to the corresponding taxonomy table for the above
# defined_community=0_atcc_species.99.txt    # Defined community we're aiming for
# conda activate $conda_qiime2   # Activate the qiime2 environment
# bash 1_ProcessSequenceToBiom.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores
# 
# # # Look at how well coverage worked
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/clustered-seqs.$otu_identity.biom.txt --taxonomy $silva_taxonomy --targets $defined_community -o $workdir/1f_defined_community_metrics -k $keyfile


# #####################
# # Alternative: Do de novo clustering instead, at 97%. NOTE: Just copied over 1d_derep-seqs.qza from the default folder since is the same up to that point and ran from there
# #####################
# 
# workdir=1_AssignOtus_denovo97               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi
# otu_identity=0.97
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=NONE # Path to QIIME-compatible SILVA 16s RNA library at 99% identity, or NONE to do de novo clustering
# silva_taxonomy=NONE # Path to the corresponding taxonomy table for the above
# defined_community=0_atcc_species.97.txt    # Defined community we're aiming for
# conda activate $conda_qiime2   # Activate the qiime2 environment
# bash 1_ProcessSequenceToBiom.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores
# 
# # # Look at how well coverage worked
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/clustered-seqs.$otu_identity.biom.txt --taxonomy $silva_taxonomy --targets $defined_community -o $workdir/1f_defined_community_metrics -k $keyfile


# ######################
# # Checking chimeras
# #####################

# This involved more manual work because I didn't want to hack my pipeline for a single test. I used the following command to get chimera sequences in the default folder:

# qiime vsearch uchime-ref --i-sequences 1e_clustered-seqs.0.99.qza --i-table 1e_clustered-table.0.99.qza --i-reference-sequences 1e_reference-seqs.qza  --p-threads 7 --output-dir 99_chimeras

# Then pulled out the names of the non-chimeric sequences and manually filtered for them in R with my 1_AssessAccuracy.r script. Turns out it only changed the matching percent by 0.1% at the Genus level 
#    and ~1.5% at the species level, so probably not making a huge difference. (Also only removed ~2% of the reads, so most things were not pegged as chimeric.)



# #####################
# # Alternative: Use only 5 OTUs (2 Acinetobacter and 3 Bacteroides) to see how different these really are; started by copying derep_seqs again since didn't change those parameters
# #####################
# 
# # Clustering parameter to group OTUs; 0.99 = 99% identity, 0.97 = 97% identity, etc. Used in multiple steps
# workdir=1_AssignOtus_only5               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi

# otu_identity=0.99
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=0_silva_99_short.5.fna  # Reference sequences to use for OTU clustering
# defined_community=0_atcc_species.99.txt    # Defined community we're aiming for
# conda activate $conda_qiime2   # Activate the qiime2 environment
# bash 1_ProcessSequenceToBiom.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores
# 
# # Look at how well coverage worked
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/clustered-seqs.$otu_identity.biom.txt --taxonomy $silva_taxonomy --targets $defined_community -o $workdir/1f_defined_community_metrics  -k $keyfile

# Have to manually run vsearch to get list of aligned sequences
#conda activate $conda_qiime2   # Activate the qiime2 environment
# vsearch --usearch_global $workdir/derep-seqs.fasta --id $otu_identity --db $references --uc $workdir/1g_uclust_otus.txt --strand plus --qmask none \
#    --notmatched $workdir/1g_uclust_nonmatched.fna --threads $ncores # --msaout $workdir/1g_uclust_aligned.msa -> this option didn't actually output
#cut -f9-10 $workdir/1g_uclust_otus.txt | grep -v "*$" | cut -f1 > $workdir/1g_uclust_seqnames.txt # Get names of seqs that matched; only get ones without "*" in last column
# conda activate  # go back to default, where seqkit is installed
# cat $workdir/derep-seqs.fasta | seqkit grep --pattern-file $workdir/1g_uclust_seqnames.txt | head -n 6000 > $workdir/1h_short_targets.fasta # Each seq is 6 lines, so this is first 1000
# clustalo -i $workdir/1h_short_targets.fasta --distmat-out=$workdir/1h_short_targets.dists.txt --full --force --threads $ncores > /dev/null  # Make distance matrix of reads
# Rscript -e "library(ape); dists=read.table('$workdir/1h_short_targets.dists.txt', skip=1, header=F, row.names=1); mytree = nj(as.dist(dists))" \
#   -e "clusters=kmeans(dists, centers=2); groups=split(names(clusters\$cluster), f=clusters\$cluster); subtrees = lapply(groups, keep.tip, phy=mytree)" \
#   -e "png('$workdir/1h_trees.png', width=12, height=4, units='in', res=300); par(mfrow=c(1,3))" \
#   -e "plot(mytree, type='unrooted', show.tip.label=FALSE, main='all')" \
#   -e "plot(subtrees[[1]], type='unrooted', show.tip.label=FALSE, main='group 1')" \
#   -e "plot(subtrees[[2]], type='unrooted', show.tip.label=FALSE, main='group 2')" \
#   -e "dev.off()"

# Final: Tree clearly separates the 2 true OTUs, but does not clearly separate the 2-3 groups within that QIIME is matching. Need to look at these in more detail.


# #####################
# # Alternative: Use only 3 OTUs (all Bacteroides) to see why are clustereing differently; started by copying derep_seqs again since didn't change those parameters
# #####################
# 
# # Clustering parameter to group OTUs; 0.99 = 99% identity, 0.97 = 97% identity, etc. Used in multiple steps
# workdir=1_AssignOtus_onlyBacteroides               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi
# 
# otu_identity=0.99
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=0_silva_99_short.bacteroides.fna  # Reference sequences to use for OTU clustering
# defined_community=0_atcc_species.99.txt    # Defined community we're aiming for
# conda activate $conda_qiime2   # Activate the qiime2 environment
# bash 1_ProcessSequenceToBiom.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores
# 
# # Look at how well coverage worked
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/clustered-seqs.$otu_identity.biom.txt --taxonomy $silva_taxonomy --targets $defined_community -o $workdir/1f_defined_community_metrics  -k $keyfile


# # Have to manually run vsearch to get list of aligned sequences
# conda activate $conda_qiime2   # Activate the qiime2 environment
# vsearch --usearch_global $workdir/derep-seqs.fasta --id $otu_identity --db $references --uc $workdir/1g_uclust_otus.txt --strand plus --qmask none \
#    --notmatched $workdir/1g_uclust_nonmatched.fna --threads $ncores # --msaout $workdir/1g_uclust_aligned.msa -> this option didn't actually output
# cut -f9-10 $workdir/1g_uclust_otus.txt | grep -v "*$" | cut -f1-2 > $workdir/1g_uclust_seqnames.txt # Get names of seqs that matched and which OTU they match; only get ones without "*" in last column

# Looking at these in Jalview, it looks like the separation into 3 OTUs is occuring based on a 2 nt locations (58 and 180); everything else is just noise



#####################
# Alternative: Time to try Deblur and see if it works. Copied over 1c_quality-filtered.qza fromd efault because is same up to that point rather than rerunning from scratch
#####################
 
# # Clustering parameter to group OTUs; 0.99 = 99% identity, 0.97 = 97% identity, etc. Used in multiple steps
# workdir=1_AssignOtus_deblur               # Folder for OTU-calling steps with QIIME
# if [ ! -e $workdir ]; then mkdir $workdir; fi
# 
# otu_identity=0.99
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=$silva_seqs  # Reference sequences to use for OTU clustering
# defined_community=0_atcc_species.99_deblur.txt    # Defined community we're aiming for
# #conda activate $conda_qiime2   # Activate the qiime2 environment
# # bash 1_ProcessSequenceToBiom_Deblur.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores
# 
# # TODO: Run quality assessment script on above to see how well it worked. (Eyeballing seemed good. Check results firmly and by the numbers.) NEED TO FIX TAXONOMY
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/deblur-seqs.biom.txt --taxonomy $workdir/taxonomy_formatted.tsv --targets $defined_community -o $workdir/1f_defined_community_metrics  -k $keyfile

#####################
# Alternative: Testing Deblur off one side only to see if it fixes throwing out as many sequences
#####################
 
# Clustering parameter to group OTUs; 0.99 = 99% identity, 0.97 = 97% identity, etc. Used in multiple steps
workdir=1_AssignOtus_deblur_nojoin               # Folder for OTU-calling steps with QIIME
if [ ! -e $workdir ]; then mkdir $workdir; fi

otu_identity=0.99
rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
references=$silva_seqs  # Reference sequences to use for OTU clustering
defined_community=0_atcc_species.99_deblur.txt    # Defined community we're aiming for
conda activate $conda_qiime2   # Activate the qiime2 environment
bash 1_ProcessSequenceToBiom_Deblur_nojoin.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores

# # TODO: Run quality assessment script on above to see how well it worked. (Eyeballing seemed good. Check results firmly and by the numbers.) NEED TO FIX TAXONOMY
# conda activate $conda_phyloseq
# Rscript 1_AssessAccuracy.r --biom $workdir/deblur-seqs.biom.txt --taxonomy $workdir/taxonomy_formatted.tsv --targets $defined_community -o $workdir/1f_defined_community_metrics  -k $keyfile


# # # #####################
# # # # Alternative: Testing Dada2 to see how it performs because Deblur is throwing out most of my reads in the actual data (esp. blocking oligos)
# # # #####################
# # #  
# # # 
# # # # # Create QIIME Conda environment (for reproducibility; only has to be done once, so leave commented out most of the time. This was taken from the QIIME2 installation page at https://docs.qiime2.org/2019.7/install/native/#install-qiime-2-within-a-conda-environment)
# # # wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
# # # conda env create -n $conda_dada2 --file qiime2-2019.7-py36-linux-conda.yml
# # # rm qiime2-2019.7-py36-linux-conda.yml # OPTIONAL CLEANUP
# # # #conda activate $conda_dada2
# # # #conda install r-essentials r-base #r-rcpp=1.0.1
# # #  
# # # # # # Clustering parameter to group OTUs; 0.99 = 99% identity, 0.97 = 97% identity, etc. Used in multiple steps
# # # # workdir=1_AssignOtus_dada2               # Folder for OTU-calling steps with QIIME
# # # # if [ ! -e $workdir ]; then mkdir $workdir; fi
# # # # # 
# # # # otu_identity=0.99
# # # # rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# # # # min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# # # # references=$silva_seqs  # Reference sequences to use for OTU clustering
# # # # defined_community=0_atcc_species.99_deblur.txt    # Defined community we're aiming for
# # # # conda activate $conda_dada2   # Activate the qiime2 environment
# # # # bash 1_ProcessSequenceToBiom_Dada2.sh $datadir $workdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $references $otu_identity $ncores
# # # # 
# # # # # TODO: Run quality assessment script on above to see how well it worked. (Eyeballing seemed good. Check results firmly and by the numbers.) NEED TO FIX TAXONOMY
# # # # conda activate $conda_phyloseq
# # # # Rscript 1_AssessAccuracy.r --biom $workdir/deblur-seqs.biom.txt --taxonomy $workdir/taxonomy_formatted.tsv --targets $defined_community -o $workdir/1f_defined_community_metrics  -k $keyfile


 



# TODO: I really need to get the actual 16s sequences to compare against. Check "History" tab of ATCC accessions for NCBI data. Will have to manually extract some
