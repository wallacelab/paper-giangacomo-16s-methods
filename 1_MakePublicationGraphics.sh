#! /usr/bin/env bash 

# Create the final publication graphics
# Note: After making these files, small formatting and layout tweaks were done to the SVG files to prepare them for final publication


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
# FIGURE CREATION - Extractions
##############

extractions=$extractdir/2_Analysis/2f_otu_table.no_organelles.RDS
mytype='extraction'
rarefaction=1000

# # Figure - PCoA of extraction methods
# Rscript Extractions_PCoA.r -i $extractions --rarefaction $rarefaction -o $figdir/ExtractionPCoA --type $mytype

# # Figure - Extraction alpha diversity
# Rscript Extractions_AlphaDiversity.r -i $extractions -o $figdir/ExtractionAlphaDiversity --type $mytype

# # Figure - Fraction total/unique OTUs captured by each extraction method
# Rscript Extractions_FractionOtusCaptured.r -i $extractions -o $figdir/ExtractionUniqueSharedOtus --group-by Genus --type $mytype
# 
# # Figure - Taxa discriminated against by each method
levels="Domain Phylum Class Order Family Genus"
fit_by_mean="Domain"
alpha=0.01
# Rscript Both_TaxaDiscrimination.r -i $extractions -o $figdir/ExtractionTaxaDiscrimination.rarefy --reference PowerSoil --levels $levels --fix-zeros --alpha $alpha --type $mytype --mean-fits $fit_by_mean --rarefy --seed 1

# # Supplemental Figure - Confirming species identity of extraction samples (specifically, that maize-powersoil is actually maize and not Arabidopsis)
# blast_results=$troubleshootdir/CheckSpeciesByBlast/1_*.blast.txt
# kraken_results=$troubleshootdir/KrakenCheckExtractions/0a_kraken_report.*.txt
# min_cutoff=0.05 # Taxon has to be at least this fraction of total in at least 1 sample to be displayed (=weed out the rare stuff)
# taxonomy=~/Projects/0_RawData/Silva_132_release/majority_taxonomy_7_levels.99.txt
# keyfile=$extractdir/16s_extractions_keyfile.tsv
# Rscript Extractions_SpeciesCheckKrakenBlast.r --blast-results $blast_results --kraken-results $kraken_results --min-cutoff $min_cutoff --taxonomy $taxonomy \
#     --keyfile $keyfile --rds-file $extractions -o $figdir/ExtractionSpeciesCheck --type $mytype


##############
# FIGURE CREATION - Primer Amplification Sets
##############

with_organelles=$primerdir/2_Analysis/2b_filtered_data.phyloseq.RDS
no_organelles=$primerdir/2_Analysis/2f_otu_table.no_organelles.RDS
mytype="amplification"
rarefaction=500
# 
# # Figure - Fraction organelle DNA by primer set
# Rscript Primers_OrganelleContamination.r -i $with_organelles -o $figdir/PrimerOrganelleContamination --type $mytype
# 
# # Figure - Primer set PCoA
# Rscript Primers_PCoA.r -i $no_organelles -o $figdir/PrimerPCoA --rarefaction $rarefaction --type $mytype
# 
# # Figure - Primer set fraction OTUs captured 
# Rscript Primers_FractionOtusCaptured.r -i $no_organelles -o $figdir/PrimerUniqueSharedOtus --group-by Genus --type $mytype
# 
# # Figure - Taxa discriminated against by each method
levels="Domain Phylum Class Order Family Genus"
fit_by_mean="Domain"
alpha=0.01
Rscript Both_TaxaDiscrimination.r -i $no_organelles -o $figdir/PrimerTaxaDiscrimination.rarefy --reference Universal --levels $levels --fix-zeros --alpha $alpha --type $mytype --sample-type "Soil 1" "Soil 2" --mean-fits $fit_by_mean --rarefy --seed 1



