#! /bin/bash

# Perform quality control and filtering of 16s reads

# Arguments
scriptdir=$1    # Directory of support scripts
datadir=$2      # Directory of raw sequence reads
workdir=$3      # Working directory to output results
keyfile=$4      # QIIME-formatted keyfile
forward_primers=$5  # Space-separated list of forward primers to remove from reads
reverse_primers=$6  # Space-separated list of reverse primers to remove from reads
rev_trim=$7         # How many bases to trim off the 3' end of reverse reads
min_length=$8   # Minimum length of sequence to keep after joining
truncqual=$9    # Truncate sequences at the first base with this score or lower in the join-pairs segment
num_cores=${10} # number of CPU cores to use in cutadapt

# Set up subdirectories
trimdir=$workdir/1_trimmed_seqs

if [ ! -e $trimdir ]; then mkdir $trimdir; fi


#######################
# Step 1 - Basic sequence cleaning
#######################

# Trim the last X bases from the reverse reads, since this seems to be the source of most read-joining errors
# Just copy over the forward reads
cp $datadir/*_R1_*.fastq.gz $trimdir
# Trim reverse reads TODO: Make the length a parameter set by the master script
for rev in $datadir/*_R2_*.fastq.gz; do
    filename=`basename $rev`
    cutadapt --cut -$rev_trim --cores $num_cores -o "$trimdir/$filename" "$rev"
done


# # Import sequences into QIIME (assumes Casava 1.8 paired-end demultiplexed data, which is what the original data was in; see https://docs.qiime2.org/2019.7/tutorials/importing/#importing-seqs )
echo "Importing FASTQ sequences into QIIME artifact"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $trimdir \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $workdir/1_demux-paired-end.qza
  

# # Use cutadapt to remove the primer sequences from forward and reverse reads
echo "Removing primers with cutadapt using $num_cores processing cores"
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences $workdir/1_demux-paired-end.qza \
    --o-trimmed-sequences $workdir/1a_cutadapt-paired-end.qza \
    --p-front-f $forward_primers \
    --p-front-r $reverse_primers \
    --p-cores $num_cores

    
# Join paired ends with vsearch
echo "Joining paired reads"
qiime vsearch join-pairs \
  --i-demultiplexed-seqs $workdir/1a_cutadapt-paired-end.qza \
  --o-joined-sequences $workdir/1b_joined-paired-end.qza \
  --p-minmergelen $min_length \
  --p-truncqual $truncqual \
  --verbose 
  
