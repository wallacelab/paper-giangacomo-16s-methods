#! /usr/bin/env bash 

##############
# This script was used to determine the optimal trim length to reduce the reverse read by to result in maximum read joining
# In short, some samples peaked at trimming 150 reads. Others continued to climb all the way out to 200, but this made the first group crash to nearly 0% joins, so 150 was chosen
##############



# This is the master script to analyze the data for Giangacomo et al 2019, on comparing methods for extracting and amplifying microbiomes from plant-associated samples
# This script was run in Linux Mint 18.3 with KDE 64-bit, although the majority of actual computation is meant to take place within a Conda environment

############
# SETUP
############

# List of subdirectories
datadir=0_data                      # Folder for raw sequencing files and keyfile to go in
scriptdir=0_scripts                 # Folder containing all necessary analysis scripts
qiimedir=1_AssignOtus               # Folder for OTU-calling steps with QIIME

# Make needed subdirectories
if [ ! -e $qiimedir ]; then mkdir $qiimedir; fi


# Core files and options
keyfile=$datadir/16s_primertest_keyfile.tsv    # QIIME-formatted keyfile of samples and metadata
silva_db=silva-132-99-nb-classifier.qza # Path to QIIME-compatible SILVA 16s RNA library
ncores=7    # Number of processor cores to use

# All primer sequences used to generate this data (Note: some include Illumina linkers) TODO: Make sure cutadapt handles these properly
forward_primers="AACMGGATTAGATACCCKG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAACMGGATTAGATACCCKG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGACCMGGATTAGATACCCKG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGAGGCAGCAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGYCAGCMGCCGCGGTAA"
reverse_primers="ACGTCATCCCCACCTTCC GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGACGTCATCCCCACCTTCC GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGGGTTGCGCTCGTTG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACHVGGGTWTCTAAT GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACNVGGGTWTCTAAT"


###################
# CREATE CONDA ENVIRONMENTS
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_qiime2=qiime2-2019-7      # Conda environment with QIIME 2019.7

# # # Create Conda environment (for reproducibility; only has to be done once, so leave commented out most of the time. This was taken from the QIIME2 installation page at https://docs.qiime2.org/2019.7/install/native/#install-qiime-2-within-a-conda-environment)
# # wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
# # conda env create -n $conda_qiime2 --file qiime2-2019.7-py36-linux-conda.yml
# # rm qiime2-2019.7-py36-linux-conda.yml # OPTIONAL CLEANUP

# Load functions required to be able to activate conda environments within scripts. 
# KEEP THESE UNCOMMENTED
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda
# # # export -f conda __conda_activate __conda_hashr __add_sys_prefix_to_path # Have to export these functions to be able to switch conda environments in subshells

##############
# ANALYSIS PIPELINE
##############

conda activate $conda_qiime2   # Activate the qiime2 environment



# Quick check to see which trim is best 
for rev_trim in 130 140 160 180; do
    tmpdir=${qiimedir}_trim$rev_trim
    if [ ! -e $tmpdir ]; then mkdir $tmpdir; fi
    bash 1_ProcessSequenceToBiom.sh $scriptdir $datadir $tmpdir $keyfile "$forward_primers" "$reverse_primers" $rev_trim $min_length $ncores 2>&1 | tee $tmpdir/log.txt
    echo $rev_trim > $tmpdir/joined_counts.txt
    grep "  Merged " $tmpdir/log.txt >> $tmpdir/joined_counts.txt
#     break
done
paste $qiimedir*/joined_counts.txt > 99_trim_counting.txt

conda deactivate
