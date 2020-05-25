#! bin/bash

# Sanity check what samples each extraction is supposed to be from

KRAKEN=/home/jgwall/Software/Metagenomics/kraken2/bin/kraken2
kraken_db=/home/jgwall/Software/Metagenomics/kraken2/kraken_plants

reads=../../TestExtraction/0_data/*/*_R1_*.fastq.gz   # Just the first read of each

# for read in $reads; do
#     sample=`basename $read`
#     sample=${sample/_*/}
#     echo "Processing $sample"
#     
#     $KRAKEN -db $kraken_db --threads 7 --report 0a_kraken_report.$sample.txt --gzip-compressed $read
# #     break
# done


python3 0_ConsolidateKraken.py -i 0a_kraken_report.*.txt -o 0b_kraken_reports --keyfile ../../TestExtraction/16s_extractions_keyfile.tsv #--debug
