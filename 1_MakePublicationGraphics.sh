#! /usr/bin/env bash 

# Create the final publication graphics

# TODO: Find taxa most discriminated against by the different extraction/primer methods

##############
# SETUP
##############

# Directories for each set of experiments
extractdir=TestExtraction
primerdir=TestPrimers
figdir=Figures

# Create figure directory if needed
if [ ! -e $figdir ]; then mkdir $figdir; fi


##############
# CONDA ENVIRONMENTS
##############

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_qiime2=qiime2-2019-7      # Conda environment with QIIME 2019.7
conda_phyloseq=phyloseq-1.28.0  # Conda environment with Phyloseq 1.28.0

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

conda activate $conda_phyloseq

##############
# FIGURE CREATION
# Note: After making these files, small formatting and layout tweaks were done to the SVG files to prepare them for final publication
##############

# # Figure - PCoA of extraction methods
# dataset=$extractdir/2_Analysis/2f_otu_table.no_organelles.RDS
# rarefaction=1000
# Rscript Extractions_PCoA.r -i $dataset --rarefaction $rarefaction -o $figdir/ExtractionPCoA

# # Figure - Extraction alpha diversity
# dataset=$extractdir/2_Analysis/2f_otu_table.no_organelles.RDS
# Rscript Extractions_AlphaDiversity.r -i $dataset -o $figdir/ExtractionAlphaDiversity

# # Figure - Fraction total/unique OTUs captured by each extraction method
# dataset=$extractdir/2_Analysis/2f_otu_table.no_organelles.RDS
# Rscript Extractions_FractionOtusCaptured.r -i $dataset -o $figdir/ExtractionUniqueSharedOtus --group-by Genus
# 
# # Figure - Fraction organelle DNA by primer set
# dataset=$primerdir/2_Analysis/2b_filtered_data.phyloseq.RDS
# Rscript Primers_OrganelleContamination.r -i $dataset -o $figdir/PrimerOrganelleContamination
# 
# # Figure - Primer set PCoA
# dataset=$primerdir/2_Analysis/2f_otu_table.no_organelles.RDS
# rarefaction=500
# Rscript Primers_PCoA.r -i $dataset -o $figdir/PrimerPCoA --rarefaction $rarefaction
# 
# Figure - Primer set fraction OTUs captured 
dataset=$primerdir/2_Analysis/2f_otu_table.no_organelles.RDS
Rscript Primers_FractionOtusCaptured.r -i $dataset -o $figdir/PrimerUniqueSharedOtus --group-by Genus
