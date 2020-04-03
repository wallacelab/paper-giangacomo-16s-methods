#! /bin/bash

# Perform quality control and filtering of 16s reads

# Arguments
datadir=$1      # Directory of raw sequence reads
workdir=$2      # Working directory to output results
keyfile=$3      # QIIME-formatted keyfile
forward_primers=$4  # Space-separated list of forward primers to remove from reads
reverse_primers=$5  # Space-separated list of reverse primers to remove from reads
rev_trim=$6         # How many bases to trim off the 3' end of reverse reads
min_length=$7    # Minimum length of sequence to keep after joining
references="$8"    # Classifer artifact to use for adding taxonomy to reads
otu_identity=$9  # Fraction identity for OTU cluster (so, 0.99 = 99%)
num_cores=${10}  # number of CPU cores to use in cutadapt
#silva_taxonomy=${11} # Taxonomy file for use with reference OTUs
#defined_community=${12}

# Set up subdirectories
trimdir=$workdir/1_trimmed_seqs
filtdir=$workdir/1a_filtered_seqs

if [ ! -e $trimdir ]; then mkdir $trimdir; fi
if [ ! -e $filtdir ]; then mkdir $filtdir; fi


# #######################
# # Step 1 - Basic sequence cleaning
# #######################
# 
# # # # # # # # Trim the last X bases from the reverse reads, since this seems to be the source of most read-joining errors
# # # # # # # cp $datadir/*_R1_*.fastq.gz $trimdir        # Just copy over the forward reads
# # # # # # # for rev in $datadir/*_R2_*.fastq.gz; do     # Trim reverse reads 
# # # # # # #     filename=`basename $rev`
# # # # # # #     cutadapt --cut -$rev_trim --cores $num_cores -o "$trimdir/$filename" "$rev"
# # # # # # # done
# # # # # # # 
# # # # # # # # Filter out reads with no data
# # # # # # # for fwd in $trimdir/*_R1_*.fastq.gz; do     # Trim reverse reads 
# # # # # # #     rev=${fwd/_R1_/_R2_}
# # # # # # #     fwd_out=`basename $fwd`
# # # # # # #     rev_out=`basename $rev`
# # # # # # #     cutadapt --cores $num_cores -m 1 -o "$filtdir/$fwd_out" -p "$filtdir/$rev_out" "$fwd" "$rev"    # "-m 1" means "at least 1 bp long", so filter out failed sequences
# # # # # # # done
# 
# 
# # # Import sequences into QIIME (assumes Casava 1.8 paired-end demultiplexed data, which is what the original data was in; see https://docs.qiime2.org/2019.7/tutorials/importing/#importing-seqs )
# echo "Importing FASTQ sequences into QIIME artifact"
# qiime tools import \
#   --type 'SampleData[PairedEndSequencesWithQuality]' \
#   --input-path $datadir \
#   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#   --output-path $workdir/1_demux-paired-end.qza
#   
# 
# # # Use cutadapt to remove the primer sequences from forward and reverse reads
# echo "Removing primers with cutadapt using $num_cores processing cores"
# qiime cutadapt trim-paired \
#     --i-demultiplexed-sequences $workdir/1_demux-paired-end.qza \
#     --o-trimmed-sequences $workdir/1a_cutadapt-paired-end.qza \
#     --p-front-f $forward_primers \
#     --p-front-r $reverse_primers \
#     --p-cores $num_cores

#     
# # # # # # # # Join paired ends with vsearch
# # # # # # # echo "Joining paired reads"
# # # # # # # qiime vsearch join-pairs \
# # # # # # #   --i-demultiplexed-seqs $workdir/1a_cutadapt-paired-end.qza \
# # # # # # #   --o-joined-sequences $workdir/1b_joined-paired-end.qza \
# # # # # # #   --p-minmergelen $min_length \
# # # # # # #   --verbose 
# # # # # # #   
# # # # # # # # Summary/visualization of joined pairs
# # # # # # # echo "Summarizing joined pairs"
# # # # # # # qiime demux summarize \
# # # # # # #   --i-data $workdir/1b_joined-paired-end.qza \
# # # # # # #   --o-visualization $workdir/1b_joined-paired-end.summary.qzv
# # # # # # #   
# # # # # # # # Quality filter joined pairs
# # # # # # # echo "Quality filtering joined pairs"
# # # # # # # qiime quality-filter q-score-joined \
# # # # # # #   --i-demux $workdir/1b_joined-paired-end.qza \
# # # # # # #   --o-filtered-sequences $workdir/1c_quality-filtered.qza \
# # # # # # #   --o-filter-stats $workdir/1c_quality-filtered.stats.qza


#######################
# Step 2 - Running Deblur to 
#######################

# # Run deblur
# qiime deblur denoise-16S \
#   --i-demultiplexed-seqs $workdir/1a_cutadapt-paired-end.qza  \
#   --p-trim-length 250 \
#   --p-sample-stats \
#   --p-jobs-to-start $num_cores \
#   --o-representative-sequences $workdir/1d_rep-seqs.qza \
#   --o-table $workdir/1d_deblur_table.qza \
#   --o-stats $workdir/1d_deblur_stats.qza

# # Classify sequences taxonomically
# classifier=/home/jgwall/Projects/0_RawData/Silva_132_release/silva-132-99-515-806-nb-classifier.qza
# qiime feature-classifier classify-sklearn \
#   --i-classifier $classifier \
#   --i-reads $workdir/1d_rep-seqs.qza \
#   --o-classification $workdir/1d_rep-seqs.taxonomy.qza
# 
# qiime metadata tabulate \
#   --m-input-file $workdir/1d_rep-seqs.taxonomy.qza \
#   --o-visualization $workdir/1d_rep-seqs.taxonomy.qzv
# 
# 
# # Summarize the resulting feature table  
# echo "Summarizing feature table"
# qiime feature-table summarize \
#   --i-table $workdir/1d_deblur_table.qza \
#   --o-visualization $workdir/1d_deblur_table.visualization.qzv \
#   --m-sample-metadata-file $keyfile
   
  
###################
# Step 3 - Extract data for other programs to access
###################


# # Extract needed data to put everything into Phyloseq (=feature table and phylogenetic tree)
# qiime tools export --input-path $workdir/1d_deblur_table.qza --output-path $workdir
# qiime tools export --input-path $workdir/1d_rep-seqs.qza --output-path $workdir
# qiime tools export --input-path $workdir/1d_rep-seqs.taxonomy.qza --output-path $workdir 
# 
# # Rename more sensibly
# mv $workdir/feature-table.biom $workdir/deblur-seqs.biom
# mv $workdir/dna-sequences.fasta $workdir/rep-seqs.fasta
# 
# # Convert biom table to text (less likely to have import issues into phyloseq)
# biom convert -i $workdir/deblur-seqs.biom -o $workdir/deblur-seqs.biom.txt --to-tsv
# 
# # Reformat taxonomy to match SILVA style
# cut -f1-2 $workdir/taxonomy.tsv | tail -n +2 > $workdir/taxonomy_formatted.tsv


# # Make summary graphic of how many reads made it into the final table (for use in adjusting parameters)
qiime tools export --input-path $workdir/1d_deblur_stats.qza --output-path $workdir
mv $workdir/stats.csv $workdir/deblur-stats.csv
python3 1g_GraphDeblurReadRetention.py -i $workdir/deblur-stats.csv -o $workdir/deblur-stats.mapped_reads.png
