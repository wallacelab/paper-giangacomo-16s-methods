#! /bin/bash 

# Several of the samples have the wrong length amplicons. I'm trying to check them for primers to see if I can figure out what happened (ie, if a simple mislabel or something)

datadir=../../0_data
all=0_primers.fa
fwd=0_primers.fwd.fa
rev=0_primers.rev.fa
maxprocs=7

workdir=1_PrimerCheck
if [ ! -e $workdir ]; then mkdir $workdir; fi

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED
conda activate qiime2-2019-7


# Use Cutadapt to get adapter sequences
for fastq in $datadir/*/*.fastq.gz; do
    
    # Get sample name
    sample=`basename $fastq`
    sample=${sample/_001.fastq.gz/}
    sample=${sample/_S*_L001/}
    
    # Run cutadapt and get primer report (cut to first 50 bp so focus on primer region)
    zcat $fastq | cut -c1-50 | cutadapt -a file:$all -o /dev/null --report minimal --info-file $workdir/1_$sample.txt -

    cut -f8 $workdir/1_$sample.txt | sort | uniq -c > $workdir/1a_$sample.collated.txt
    rm $workdir/1_$sample.txt   # Big file, no longer needed
    
#     break
    
done

# Collate in a python file and plot heatmap of results
python3 2_PlotPrimersInSamples.py -i $workdir/1a_*.collated.txt -o 2_primer_sets 
