#! /usr/bin/env bash 

biom=1_AssignOtus/feature-table.biom
keyfile=16s_primertest_keyfile.tsv    # QIIME-formatted keyfile of samples and metadata
taxonomy=/home/jgwall/Projects/0_RawData/Silva_132_release/majority_taxonomy_7_levels.txt # Path to the corresponding taxonomy table for the above


# Set up directory
workdir=99_phinch
if [ ! -e $workdir ]; then mkdir $workdir; fi

# # Convert taxonomy
# echo -e "#OTUID\ttaxonomy" > $workdir/taxonomy.txt
# cat $taxonomy >> $workdir/taxonomy.txt


# Convert file
# biom convert -i $biom -o $workdir/feature-table.biom.tsv --to-tsv
biom convert -i $workdir/feature-table.biom.tsv -o $workdir/feature-table.biom --to-hdf5 --sample-metadata-fp $keyfile --observation-metadata-fp $workdir/taxonomy.txt #--process-obs-metadata sc_separated


# biom convert -i $workdir/feature-table.biom.tsv -o $workdir/feature-table.biom --to-hdf5  
