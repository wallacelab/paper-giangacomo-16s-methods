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
troubleshootdir=0_Troubleshooting   # Specifically for checking species identity by BLAST and Kraken

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

# TODO: with standard renaming, see if can combine scripts and only use one for both extractions and amplifications instead of current duplication

# # # Extractions # # # 

extractions=$extractdir/2_Analysis/2f_otu_table.no_organelles.RDS

# # Figure - PCoA of extraction methods
# rarefaction=1000
# Rscript Extractions_PCoA.r -i $extractions --rarefaction $rarefaction -o $figdir/ExtractionPCoA

# # Figure - Extraction alpha diversity
# Rscript Extractions_AlphaDiversity.r -i $extractions -o $figdir/ExtractionAlphaDiversity

# # Figure - Fraction total/unique OTUs captured by each extraction method
# Rscript Extractions_FractionOtusCaptured.r -i $extractions -o $figdir/ExtractionUniqueSharedOtus --group-by Genus

# Figure - Taxa discriminated against by each method
levels="Phylum Class Order Family Genus"
alpha=0.01
Rscript Both_TaxaDiscrimination.r -i $extractions -o $figdir/ExtractionTaxaDiscrimination --reference PowerSoil --levels $levels --fix-zeros --alpha $alpha --type extraction

# # Supplemental Figure - Confirming species identity of extraction samples (specifically, that maize-powersoil is actually maize and not Arabidopsis)
# blast_results=$troubleshootdir/CheckSpeciesByBlast/1_*.blast.txt
# kraken_results=$troubleshootdir/KrakenCheckExtractions/0a_kraken_report.*.txt
# min_cutoff=0.05 # Taxon has to be at least this fraction of total in at least 1 sample to be displayed (=weed out the rare stuff)
# taxonomy=~/Projects/0_RawData/Silva_132_release/majority_taxonomy_7_levels.99.txt
# keyfile=$extractdir/16s_extractions_keyfile.tsv
# rds_file=$extractdir/2_Analysis/2f_otu_table.no_organelles.RDS
# Rscript Extractions_SpeciesCheckKrakenBlast.r --blast-results $blast_results --kraken-results $kraken_results --min-cutoff $min_cutoff --taxonomy $taxonomy \
#     --keyfile $keyfile --rds-file $rds_file -o $figdir/ExtractionSpeciesCheck



# # # Primer Sets # # #

with_organelles=$primerdir/2_Analysis/2b_filtered_data.phyloseq.RDS
no_organelles=$primerdir/2_Analysis/2f_otu_table.no_organelles.RDS

# # Figure - Fraction organelle DNA by primer set
# Rscript Primers_OrganelleContamination.r -i $with_organelles -o $figdir/PrimerOrganelleContamination

# # Figure - Primer set PCoA TODO: Moiddle graphic legend doesn't match others for some reason.
# rarefaction=500
# Rscript Primers_PCoA.r -i $no_organelles -o $figdir/PrimerPCoA --rarefaction $rarefaction

# # Figure - Primer set fraction OTUs captured 
# Rscript Primers_FractionOtusCaptured.r -i $no_organelles -o $figdir/PrimerUniqueSharedOtus --group-by Genus

# # Figure - Taxa discriminated against by each method
# levels="Phylum Class Order Family Genus"
# alpha=0.01
# Rscript Both_TaxaDiscrimination.r -i $extractions -o $figdir/ExtractionTaxaDiscrimination --reference PNAs --levels $levels --fix-zeros --alpha $alpha


