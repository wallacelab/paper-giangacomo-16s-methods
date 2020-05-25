#! /bin/bash

# Perform quality control and filtering of 16s reads 

# Arguments
datadir=$1      # Directory of raw sequence reads
workdir=$2      # Working directory to output results
qiimedir=$3     # Directory for QIIME exported results (assumed to have default QIIME names like "dna-sequences.fasta", "feature-table.biom", etc.
taxonomy=$4     # Taxonomy key for the sequences to be analyzed (assumed to be from SILVA)
silva_tree=$5   # SILVA (or other) phylogenetic tree of the OTUs
keyfile=$6      # QIIME-formatted keyfile for the samples
sample_depth=$7 # Minimum total reads to keep a sample
otu_depth=$8    # Minimum total reads to keep an OTU
otu_prevalence=$9 # Minimum number of samples an OTU has to appear in to be kept
rarefaction=${10}    # Level to rarefy samples to before doing UniFrac analysis
bad_samples="${11}"  # List of bad samples to be filtered out
group_by=${12}

phylo_filtered=$workdir/2b_filtered_data.phyloseq.RDS

# # Make a combined phylogenetic tree across the datasets
# Rscript 2_MakePhylogeneticTreeFromBlastResults.r --blast $qiimedir/*/1g_rep-seqs.blast_results.txt --tree $silva_tree --taxonomy $taxonomy -o $workdir/2_combine_phylogenies
# 
# # Convert the QIIME-exported files into a Phyloseq object
# Rscript 2a_ConvertQiimeToPhyloseqObject.r --biom $qiimedir/*/deblur-seqs.biom.txt --tree $workdir/2_combine_phylogenies.tre \
#     --taxonomy $qiimedir/*/taxonomy_formatted.tsv --key $keyfile -o $workdir/2a_combined_data.phyloseq --samples-to-remove $bad_samples 
# 
# # Filter the phyloseq data based on sample & otu depth and otu prevalence
# Rscript 2b_FilterPhyloseqData.r -i $workdir/2a_combined_data.phyloseq.RDS --sample-depth $sample_depth --otu-depth $otu_depth --otu-prevalence $otu_prevalence \
#     -o $phylo_filtered --outgraphics $workdir/2b_filtered_data.reports.png --group-by $group_by
#     
# # # Determine how much plastid sequence is in each sample
# Rscript 2c_AssessPlastidContamination.r -i $phylo_filtered -o $workdir/2c_plastid_contamination
#   
# # # Determine how much each extraction method distorts the community taxonomy
# Rscript 2d_AssessCommunityDistortion_taxonomy.r -i $phylo_filtered -o $workdir/2d_taxonomy_distortion
# 
# # Determine how much each primer set distorts the community using Principal Coordinates analysis
# # NOTE: Phyloseq's UniFrac calculation appears to be wrong, so using rbiom instead
# Rscript 2e_AssessCommunityDistortion_PCoA.r -i $phylo_filtered --rarefaction $rarefaction -o $workdir/2e_pcoa_distortion
# 
# # Remove organelle sequences and recheck MDS plots
# Rscript 2f_RemoveOrganelles.r -i $phylo_filtered --outprefix $workdir/2f_otu_table.no_organelles
# Rscript 2e_AssessCommunityDistortion_PCoA.r -i $workdir/2f_otu_table.no_organelles.RDS --rarefaction $rarefaction -o $workdir/2g_pcoa_distortion.no_organelles
# 
# # Plot alpha diversity and shared/unique OTUs
# Rscript 2h_AnalyzeAlphaDiversity.r -i $phylo_filtered -o $workdir/2h_filtered
# Rscript 2h_AnalyzeAlphaDiversity.r -i $workdir/2f_otu_table.no_organelles.RDS -o $workdir/2h_filtered.no_organelles

# Sanity check to make sure each sample comes from what we think it does (esp. since have outliers); use mitochondria and chloroplast seqs to identify origins
#    Overall looks good. M99 (maize) looks a bit like an outlier with lots of Camelina salvia mitochondria (like Arabidopsis), but the levels are still 20x lower than actual Arabidopsis samples, so probably okay.
# Rscript 2i_CheckOrganelleSpecies.r -i $phylo_filtered -o $workdir/2i_organelle_check.png

# Repeat last step with vsearch-clustered OTUs, which in preliminary analyses gave clearer identification of each sample
# Rscript 2a_ConvertQiimeToPhyloseqObject.r --biom $qiimedir/*/vsearch-seqs.biom.txt --tree $silva_tree \
#      --taxonomy $taxonomy --key $keyfile -o $workdir/2j_combined_vsearch.phyloseq --samples-to-remove $bad_samples 
# Rscript 2i_CheckOrganelleSpecies.r -i $workdir/2j_combined_vsearch.phyloseq.RDS -o $workdir/2j_organelle_check.vsearch.png
# Ultimately the above species checks are ambiguous, depending on whether you look at identity, sequence count, etc. The BLAST and Kraken results in the troubleshooting directory seem more reliable

# Look for taxa being discriminated against
Rscript 2k_FindDiscriminatedTaxa.r -i $workdir/2f_otu_table.no_organelles.RDS -o $workdir/2k_discrimination --reference MoBioPowerSoil --levels Phylum Class Order Family Genus Species
# TODO: Adding DESeq seems to have broken my conda environment. *sigh* Fix that next time
