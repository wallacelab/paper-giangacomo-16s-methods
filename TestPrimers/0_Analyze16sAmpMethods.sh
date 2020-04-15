#! /usr/bin/env bash 

# This is the master script to analyze the data for Giangacomo et al 2019, on comparing methods for extracting and amplifying microbiomes from plant-associated samples
# This script was run in Linux Mint 18.3 with KDE 64-bit, although the majority of actual computation is meant to take place within a Conda environment

############
# SETUP
############

# List of subdirectories
datadir=0_data                      # Folder for raw sequencing files and keyfile to go in
qiimedir=1_AssignOtus               # Folder for OTU-calling steps with QIIME
analdir=2_Analysis                  # Folder for initial phyloseq-derived analyses

# Make needed subdirectories
if [ ! -e $qiimedir ]; then mkdir $qiimedir; fi
if [ ! -e $analdir ]; then mkdir $analdir; fi


# Core files and options
keyfile=16s_primertest_keyfile.tsv    # QIIME-formatted keyfile of samples and metadata
analysis_keys=16s_primer_sets.tsv      # Tab-separated list of infomration for each Deblur analysis: dataset name (also the subdirectory within 0_data), forward primer, and reverse primer
SILVA=/home/jgwall/Projects/0_RawData/Silva_132_release # Base directory to SILVA data
silva_seqs=$SILVA/silva_132_99_16S.fna # Path to QIIME-compatible SILVA 16s RNA library at 99% identity
silva_taxonomy=$SILVA/majority_taxonomy_7_levels.99.txt # Path to the corresponding taxonomy table for the above
silva_tree=$SILVA/silva_99_otus.tre # SILVA phylogenetic tree
primers=0_primers.fa
ncores=7    # Number of processor cores to use

###################
# CREATE CONDA ENVIRONMENTS
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_qiime2=qiime2-2019-7      # Conda environment with QIIME 2019.7
conda_phyloseq=phyloseq-1.28.0  # Conda environment with Phyloseq 1.28.0

# # # Create QIIME Conda environment (for reproducibility; only has to be done once, so leave commented out most of the time. This was taken from the QIIME2 installation page at https://docs.qiime2.org/2019.7/install/native/#install-qiime-2-within-a-conda-environment)
# # wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
# # conda env create -n $conda_qiime2 --file qiime2-2019.7-py36-linux-conda.yml
# # rm qiime2-2019.7-py36-linux-conda.yml # OPTIONAL CLEANUP

# # # Create Phyloseq Conda environment (for reproducibility; only has to be done once, so leave commented out most of the time. This was taken/modified from https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/)
# note: this requires the bioconda channel to be included; see https://bioconda.github.io/user/install.html#set-up-channels for instructions
# conda create -n $conda_phyloseq r-essentials r-base bioconductor-phyloseq


# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

##############
# ANALYSIS PIPELINE
##############

conda activate $conda_qiime2   # Activate the qiime2 environment

# # Step 0: Import SILVA sequences for use in taxonomic classification later
# qiime tools import \
#   --type 'FeatureData[Sequence]' \
#   --input-path $silva_seqs \
#   --output-path $qiimedir/reference_seqs.qza
# 
# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-format HeaderlessTSVTaxonomyFormat \
#   --input-path $silva_taxonomy \
#   --output-path $qiimedir/reference_taxonomy.qza
# 
# makeblastdb -in $silva_seqs -out $qiimedir/reference_seqs -dbtype nucl


# # Step 1: Process raw sequence data down to BIOM feature table
# rev_trim=150            # How many bases to trim off the 3' end of Reverse reads (where most errors occur); some empirical checks early on found 150 to be the optimum (though some kept improving slgihtly to at least 200, others crashed and couldn't be joined)
# min_length=100          # Minimum length of sequence to keep after joining (msotly a sanity check, since true reads should be much longer)
# references=$silva_seqs  # Reference sequences to use for OTU clustering
# while read primer_set fwd rev deblur_length; do
#     workdir=$qiimedir/$primer_set
#     if [ ! -e $workdir ]; then mkdir $workdir; fi
#     
#     echo "### Processing $primer_set down to biom with Deblur ###"
#     bash 1_ProcessSequenceToBiom.sh $datadir/$primer_set $workdir $keyfile $fwd $rev $rev_trim $min_length $deblur_length \
#          $qiimedir/reference_seqs.qza $qiimedir/reference_taxonomy.qza $primers $qiimedir/reference_seqs $ncores
# #     break
#     
# done < $analysis_keys


# Step 2: Use the Phyloseq package in R to perform analysis and visualization on the raw data
sample_depth=500
otu_depth=3
otu_prevalence=2
rarefaction=500
bad_samples="BO57L2-5 BO57NOL2-5 DISC2-5 PNA3-1 PNA3-2 PNA3-3 PNA3-4 PNA3-5 PNA3-6" # Samples found to be bad in preliminary analyses
grouping="Genus"
conda activate $conda_phyloseq
bash ./2_AnalyzePrimerSets.sh $datadir $analdir $qiimedir $silva_taxonomy $silva_tree $keyfile $sample_depth $otu_depth $otu_prevalence $rarefaction "$bad_samples" $grouping



