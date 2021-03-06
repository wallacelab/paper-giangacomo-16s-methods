#! /usr/bin/env bash 

# Create the Conda environments used in running the specific environments. This script only needs to be run once, before starting any of the individual analyses


# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_qiime2=qiime2-2019-7      # Conda environment with QIIME 2019.7
conda_phyloseq=phyloseq-1.28.0  # Conda environment with Phyloseq 1.28.0

##############
# Conda environment creation as done originally - DO NOT USE
##############

# # Create QIIME Conda environment (for reproducibility; only has to be done once, so leave commented out most of the time. This was taken from the QIIME2 installation page at https://docs.qiime2.org/2019.7/install/native/#install-qiime-2-within-a-conda-environment)
# wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
# conda env create -n $conda_qiime2 --file qiime2-2019.7-py36-linux-conda.yml
# rm qiime2-2019.7-py36-linux-conda.yml # OPTIONAL CLEANUP

# Create Phyloseq Conda environment; This was taken/modified from https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/
# note: this requires the bioconda channel to be included; see https://bioconda.github.io/user/install.html#set-up-channels for instructions
# conda create -n $conda_phyloseq r-essentials=3.6.0 r-base=3.6.3 bioconductor-phyloseq=1.28.0 r-ape=5.3 r-dplyr=0.8.3 r-ggplot2=3.2.1 r-tidyr=1.0.0 r-vegan=2.5.5 bioconductor-deseq2

# # Export the environment files (recommended to use these instead when creating a new conda environment)
# conda activate $conda_qiime2
# conda env export > $conda_qiime2.yml
# conda activate $conda_phyloseq
# conda env export > $conda_phyloseq.yml

##############
# Conda environment creation with the supplied files - USE THESE ONES
##############

conda env create -n $conda_qiime2 --file qiime2-2019-7.yml
conda env create -n $conda_phyloseq --file phyloseq-1.28.0.yml


##############
# To be able to load these conda environments within scripts, you need to include the following line of code:
##############

. $(conda info --root)/etc/profile.d/conda.sh 

