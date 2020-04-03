#! /bin/bash

# Perform quality control and filtering of 16s reads

# Arguments
datadir=$1      # Directory of raw sequence reads
workdir=$2      # Working directory to output results
keyfile=$3      # QIIME-formatted keyfile
forward_primer=$4  # Forward primer to remove from reads
reverse_primer=$5  # Feverse primer to remove from reads
rev_trim=$6         # How many bases to trim off the 3' end of reverse reads
min_length=$7    # Minimum length of sequence to keep after joining
deblur_length=$8 # Length to trim sequences to for Deblur
ref_seqs=$9    # Reference sequences used to train the classifier
ref_taxonomy=${10} # Reference taxonomy used to train the classifier
primers=${11}   # Fasta file of primers to check for correct amplicons
blast_db=${12}  # BLAST database of OTUs to use for phylogenetic trees
num_cores=${13}  # number of CPU cores to use in cutadapt

# Set up subdirectories
trimdir=$workdir/1_trimmed_seqs
filtdir=$workdir/1a_filtered_seqs
primerdir=$workdir/1d_primer_checks

if [ ! -e $trimdir ]; then mkdir $trimdir; fi
if [ ! -e $filtdir ]; then mkdir $filtdir; fi
if [ ! -e $primerdir ]; then mkdir $primerdir; fi


# #######################
# # Step 1 - Basic sequence cleaning
# #######################
# 
# # Trim the last X bases from the reverse reads, since this seems to be the source of most read-joining errors
# cp $datadir/*_R1_*.fastq.gz $trimdir        # Just copy over the forward reads
# for rev in $datadir/*_R2_*.fastq.gz; do     # Trim reverse reads 
#     filename=`basename $rev`
#     cutadapt --cut -$rev_trim --cores $num_cores -o "$trimdir/$filename" "$rev"
# done
# 
# # Filter out reads with no data
# for fwd in $trimdir/*_R1_*.fastq.gz; do     # Trim reverse reads 
#     rev=${fwd/_R1_/_R2_}
#     fwd_out=`basename $fwd`
#     rev_out=`basename $rev`
#     cutadapt --cores $num_cores -m 1 -o "$filtdir/$fwd_out" -p "$filtdir/$rev_out" "$fwd" "$rev"    # "-m 1" means "at least 1 bp long", so filter out failed sequences
# done
# 
# 
# # # Import sequences into QIIME (assumes Casava 1.8 paired-end demultiplexed data, which is what the original data was in; see https://docs.qiime2.org/2019.7/tutorials/importing/#importing-seqs )
# echo "Importing FASTQ sequences into QIIME artifact"
# qiime tools import \
#   --type 'SampleData[PairedEndSequencesWithQuality]' \
#   --input-path $filtdir \
#   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#   --output-path $workdir/1_demux-paired-end.qza
#   
# 
# # # Use cutadapt to remove the primer sequences from forward and reverse reads
# echo "Removing primers with cutadapt using $num_cores processing cores"
# qiime cutadapt trim-paired \
#     --i-demultiplexed-sequences $workdir/1_demux-paired-end.qza \
#     --o-trimmed-sequences $workdir/1a_cutadapt-paired-end.qza \
#     --p-front-f "^$forward_primer" \
#     --p-front-r "^$reverse_primer" \
#     --p-cores $num_cores
# 
#     
# # Join paired ends with vsearch
# echo "Joining paired reads"
# qiime vsearch join-pairs \
#   --i-demultiplexed-seqs $workdir/1a_cutadapt-paired-end.qza \
#   --o-joined-sequences $workdir/1b_joined-paired-end.qza \
#   --p-minmergelen $min_length \
#   --verbose 
#   
# # Summary/visualization of joined pairs
# echo "Summarizing joined pairs"
# qiime demux summarize \
#   --i-data $workdir/1b_joined-paired-end.qza \
#   --o-visualization $workdir/1b_joined-paired-end.summary.qzv
#   
# # Quality filter joined pairs
# echo "Quality filtering joined pairs"
# qiime quality-filter q-score-joined \
#   --i-demux $workdir/1b_joined-paired-end.qza \
#   --o-filtered-sequences $workdir/1c_quality-filtered.qza \
#   --o-filter-stats $workdir/1c_quality-filtered.stats.qza
# 
# # Extract joined sequences to get length distributions
# echo "Extracting joined and filtered pairs for length assessment"
# qiime tools export \
#     --input-path $workdir/1c_quality-filtered.qza \
#     --output-path $workdir/1c_quality-filtered-seqs
# 
# # Get legnth distributions (sanity check that Deblur trim lengths are correct)
# echo "Determining length distributions"
# python3 1c_GetFastqLengthDistributions.py -i $workdir/1c_quality-filtered-seqs/*.fastq.gz -o $workdir/1c_quality-filtered-seqs.lengths
# rm -rf $workdir/1c_quality-filtered-seqs    # Delete extracted directory; no longer needed
# 
# # Sanity check to make sure primers in original sequences are correct (and so can filter out problem samples if they aren't)
# echo "Sanity checking primers used (note that BO3/4 and Universal/PNAs are nearly identical so often misidentified)"
# for fastq in $datadir/*.fastq.gz; do
#     # Get sample name
#     sample=`basename $fastq`
#     sample=${sample/_001.fastq.gz/}
#     sample=${sample/_S*_L001/}
#     echo $sample
#   
#     # Run cutadapt and get primer report (cut to first 50 bp so focus on primer region)
#     zcat $fastq | cut -c1-50 | cutadapt -a file:$primers -o /dev/null --report minimal --info-file $primerdir/1_$sample.txt -
# 
#     cut -f8 $primerdir/1_$sample.txt | sort | uniq -c > $primerdir/1a_$sample.collated.txt
#     rm $primerdir/1_$sample.txt   # Big file, no longer needed
#     
# #     break
#     
# done
# 
# # Collate in a python file and plot heatmap of results
# python3 1c_PlotPrimersInSamples.py -i $primerdir/1a_*.collated.txt -o $workdir/1c_primer_sets 



#######################
# Step 2 - Run Deblur to get down to ASVs
#######################

# Step 2a - ASV Assignment

# # # Run deblur
# qiime deblur denoise-16S \
#   --i-demultiplexed-seqs $workdir/1c_quality-filtered.qza  \
#   --p-trim-length $deblur_length \
#   --p-sample-stats \
#   --p-min-size 1 \
#   --p-min-reads 2 \
#   --p-jobs-to-start $num_cores \
#   --o-representative-sequences $workdir/1d_rep-seqs.qza \
#   --o-table $workdir/1d_deblur_table.qza \
#   --o-stats $workdir/1d_deblur_stats.qza

# # Summarize the resulting feature table  
# echo "Summarizing feature table"
# qiime feature-table summarize \
#   --i-table $workdir/1d_deblur_table.qza \
#   --o-visualization $workdir/1d_deblur_table.visualization.qzv \
#   --m-sample-metadata-file $keyfile


# Step 2b - Taxonomic assignment

# # Extract part of silva reads that correspond to this amplicon region
# echo "Extracting and training taxonomic classifier (this will take a while)"
# qiime feature-classifier extract-reads \
#     --i-sequences $ref_seqs \
#     --p-f-primer $forward_primer \
#     --p-r-primer $reverse_primer \
#     --p-min-length 100 \
#     --p-max-length 1000 \
#     --o-reads $workdir/1e_ref-seqs.trimmed.qza \
#     --verbose
#   
# # Plot size distribution of reads of the above
# qiime tools export \
#     --input-path $workdir/1e_ref-seqs.trimmed.qza \
#     --output-path $workdir 
# mv $workdir/dna-sequences.fasta $workdir/1e_ref-seqs.trimmed.fa
# grep -v "^>" $workdir/1e_ref-seqs.trimmed.fa | awk '{print length}' > $workdir/1e_ref-seqs.trimmed_lengths.txt
# Rscript -e "lengths=scan('$workdir/1e_ref-seqs.trimmed_lengths.txt'); png('$workdir/1e_ref-seqs.trimmed_lengths.png')" \
#     -e "hist(lengths, col='darkblue', breaks=50, main='Trimmed Length Distribution'); dev.off()"
# rm $workdir/1e_ref-seqs.trimmed.fa
# 
# # Train Classifier
# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads $workdir/1e_ref-seqs.trimmed.qza \
#   --i-reference-taxonomy $ref_taxonomy \
#   --o-classifier $workdir/1e_classifier.qza
# 
# # Classify sequences taxonomically
# qiime feature-classifier classify-sklearn \
#   --i-classifier $workdir/1e_classifier.qza \
#   --i-reads $workdir/1d_rep-seqs.qza \
#   --p-n-jobs $num_cores \
#   --o-classification $workdir/1f_rep-seqs.taxonomy.qza
# 
# qiime metadata tabulate \
#   --m-input-file $workdir/1f_rep-seqs.taxonomy.qza \
#   --o-visualization $workdir/1f_rep-seqs.taxonomy.qzv

  
###################
# Step 3 - Extract data for other programs to access
###################


# # Extract needed data to put everything into Phyloseq (=feature table and phylogenetic tree)
# qiime tools export --input-path $workdir/1d_deblur_table.qza --output-path $workdir
# qiime tools export --input-path $workdir/1d_rep-seqs.qza --output-path $workdir
# qiime tools export --input-path $workdir/1f_rep-seqs.taxonomy.qza --output-path $workdir 
# qiime tools export --input-path $workdir/1d_deblur_stats.qza --output-path $workdir
# 
# # Rename more sensibly
# mv $workdir/feature-table.biom $workdir/deblur-seqs.biom
# mv $workdir/dna-sequences.fasta $workdir/rep-seqs.fasta
# mv $workdir/stats.csv $workdir/deblur-stats.csv
# 
# # Convert biom table to text (less likely to have import issues into phyloseq)
# biom convert -i $workdir/deblur-seqs.biom -o $workdir/deblur-seqs.biom.txt --to-tsv
# 
# # Reformat taxonomy to match SILVA style
# cut -f1-2 $workdir/taxonomy.tsv | tail -n +2 > $workdir/taxonomy_formatted.tsv
# 
# # Make summary graphic of how many reads made it into the final table (for use in adjusting parameters)
# python3 1g_GraphDeblurReadRetention.py -i $workdir/deblur-stats.csv -o $workdir/deblur-stats.mapped_reads.png



###################
# Step 4 - BLAST representative seqs so can get a unified taxonomic tree across all datasets
###################

# blastn -query $workdir/rep-seqs.fasta -db $blast_db -out $workdir/1g_rep-seqs.blast_results.txt -evalue 1 -outfmt 6 -num_alignments 100 -num_threads $num_cores
