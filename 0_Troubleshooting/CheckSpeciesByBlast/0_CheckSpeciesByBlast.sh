#! bin/bash

# Sanity check what samples each extraction is supposed to be from

reads=../../TestExtraction/0_data/*/*_R1_*.fastq.gz   # Just the first read of each
silva=/home/jgwall/Projects/0_RawData/Silva_132_release/silva_132_99_16S.fna

# # Make SILVA Blast DB
# makeblastdb -in $silva -dbtype nucl -out 1_silva 

# num_reads=1000    # Number of reads to BLAST
# for read in $reads; do
#     sample=`basename $read`
#     sample=${sample/_*/}
#     echo "Processing $sample"
#     
#     zcat $read | sed -n '1~4s/^@/>/p;2~4p' | head -n `expr $num_reads \* 2` | blastn -db 1_silva -outfmt 6 -num_threads 7 -num_alignments 20 > 1_${sample}.blast.txt    # sed converts from FASTQ to FASTA
#     
#     #break
# done
 
Rscript 0a_TallyBlastData.r -i 1_*.blast.txt -o 1a_hit_tallies.txt --outgraphic 1a_hit_heatmap.png --min-total-score 250
