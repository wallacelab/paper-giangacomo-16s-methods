#! /bin/bash 


# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

conda activate qiime2-2019-7


SILVA= # Base directory to SILVA data
silva_seqs=$SILVA

vsearch --usearch_global 1_AssignOtus/derep-seqs.fasta --id 0.99 --db /home/jgwall/Projects/0_RawData/Silva_132_release/silva_132_99_16S.fna --uc 99_vsearch --strand plus --qmask none --notmatched 99_nomatch --threads 7
