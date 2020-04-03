#! /bin/bash

#############
# Set up conda environment
#############

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_qiime2=qiime2-2019-7      # Conda environment with QIIME 2019.7

# # # Create QIIME Conda environment (for reproducibility; only has to be done once, so leave commented out most of the time. This was taken from the QIIME2 installation page at https://docs.qiime2.org/2019.7/install/native/#install-qiime-2-within-a-conda-environment)
# # wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
# # conda env create -n $conda_qiime2 --file qiime2-2019.7-py36-linux-conda.yml
# # rm qiime2-2019.7-py36-linux-conda.yml # OPTIONAL CLEANUP

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

conda activate $conda_qiime2


 

################
# Run UniFrac calculations via QIIME & export for comparison
################

# Download data
wget https://docs.qiime2.org/2018.11/data/tutorials/moving-pictures/table-deblur.qza # OTU Table
wget https://docs.qiime2.org/2018.11/data/tutorials/moving-pictures/rooted-tree.qza # Tree

# Variables
orig_table=table-deblur.qza
filtered_table=table-filtered.qza
tree=rooted-tree.qza

# Remove OTUs not present in phylogeny
qiime phylogeny filter-table \
    --i-table $orig_table \
    --i-tree $tree \
    --o-filtered-table $filtered_table

# Calculate both weighted and unweighted unifrac
for metric in weighted_unifrac unweighted_unifrac; do

    qiime diversity beta-phylogenetic \
        --i-table $filtered_table \
        --i-phylogeny $tree \
        --p-metric $metric \
        --o-distance-matrix $metric.qza \
        --verbose

    # Export to text
    qiime tools export \
        --input-path $metric.qza \
        --output-path $metric

    # Move distance matrix file to more sensible location
    mv $metric/distance-matrix.tsv $metric.tsv
    rmdir $metric
done


# Export table and tree for comparison with R functions
qiime tools export --input-path $filtered_table --output-path ./
qiime tools export --input-path $tree --output-path ./
 
# Convert biom table to text (less likely to have import issues)
biom convert -i feature-table.biom -o feature-table.biom.txt --to-tsv
