#! /bin/bash

# Make blocking oligos according to Agler et al 2016, starting from scratch (basically) so I can make absolutely sure we're right

PRIMER3=/home/jgwall/Software/Sequencing/primer3-2.3.7/src/primer3_core

# Directory setup
seqdir=0_Sequences
scriptdir=0_Scripts
workdir=2_BlockingOligos
blastdir=$workdir/2a_BlastDBs
if [ ! -e $workdir ]; then mkdir $workdir; fi
if [ ! -e $blastdir ]; then mkdir $blastdir; fi

# Master variables
greengenes=/usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta
maize_mito=$seqdir/Zea_mays.AGPv4.dna.chromosome.Mt.fa
maize_chloro=$seqdir/Zea_mays.AGPv4.dna.chromosome.Pt.fa
arabidopsis_mito=$seqdir/arabidopsis_TAIR10_mitochondria.fa
arabidopsis_chloro=$seqdir/arabidopsis_TAIR10_chloroplast.fa
ecoli_16s=$seqdir/ecoli_16s.fa
block_mito=$seqdir/agler_mito_blocking_olgos.fa
block_chloro=$seqdir/agler_plastid_blocking_olgos.fa
maxprocs=7

###
# First, pull out 16s gene locations from maize and arabidopsis organelles
###

# # Make Blast DBs of maize and arabidopsis organelles
# makeblastdb -in $maize_mito -title MaizeMitochondria -out $blastdir/maize_mito -dbtype nucl   
# makeblastdb -in $maize_chloro -title MaizeChloroplast -out $blastdir/maize_chloro -dbtype nucl   
# makeblastdb -in $arabidopsis_mito -title ArabidopsisMitochondria -out $blastdir/arabidopsis_mito -dbtype nucl   
# makeblastdb -in $arabidopsis_chloro -title ArabidopsisChloroplast -out $blastdir/arabidopsis_chloro -dbtype nucl   

# # BLAST e_coli 16s and blocking oligos against each to find best fit; default output is "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# for organism in maize arabidopsis; do
#   for organelle in mito chloro; do
#     name=${organism}_${organelle}
#     task=blastn
#     blastn -query $ecoli_16s -db $blastdir/$name -out $workdir/2a_${name}_16s_hits.txt -max_target_seqs 20 -num_threads $maxprocs -outfmt 6 -task $task
#     blastn -query $block_mito -db $blastdir/$name -out $workdir/2a_${name}_mito_block_hits.txt -max_target_seqs 20 -num_threads $maxprocs -outfmt 6 -task $task
#     blastn -query $block_chloro -db $blastdir/$name -out $workdir/2a_${name}_chloro_block_hits.txt -max_target_seqs 20 -num_threads $maxprocs -outfmt 6 -task $task
# #     break
#   done
# #   break
# done
# echo "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tr ' ' '\t' > $workdir/2b_combined_blast.txt
# cat $workdir/2a_*hits.txt | sort >> $workdir/2b_combined_blast.txt 

# Quickly plot relative locations
# Rscript $scriptdir/2_QuickPlotOligoLocations.r -i $workdir/2b_combined_blast.txt -o $workdir/2b_combined_blast.locations.png
Rscript $scriptdir/2_QuickPlotOligoLocations.r -i $workdir/2b_combined_blast.manual.txt -o $workdir/2b_combined_blast.locations.manual.png  # Manual modification to split different regions within plastids, etc.


# # Step 1 - get 30mers from target (maize) chloroplast region
# python3 $scriptdir/1_DivideTargetIntoKmers.py -i $plastid -o $workdir/1_target_kmers.fa -k 30

# # Step 2 - BLAST against Greengenes
# # Make BLAST database
# makeblastdb -in $greengenes -title Greengenes -out $workdir/1a_greengenes -dbtype nucl   
# # Perform BLAST
# blastn -query $workdir/1_target_kmers.fa -db $workdir/1a_greengenes -out $workdir/1b_kmer_hits.txt -perc_identity 50 -word_size 12 -dust no -evalue 1000 -max_target_seqs 10000 -num_threads $maxprocs -outfmt 6
# pigz $workdir/1b_kmer_hits.txt

# # Step 3 - Graph hits
# python3 $scriptdir/1c_TallyKmerHits.py -i $workdir/1b_kmer_hits.txt.gz -o $workdir/1c_kmer_hits.tally

# Step 4 - Determine best oligos to use - Need Tm, hairpin likelihood, and location
# Figure out where Agler's fall in maize plastid and how close a match they are
# cat $plastid 0_agler_plastid_blocking_olgos.fa > $workdir/1d_seqs_to_align.fa
# clustalo -i $workdir/1d_seqs_to_align.fa -o $workdir/1d_agler_seqs_aligned_to_maize.clu --threads=$maxprocs --verbose --outfmt clu --force --output-order tree-order --seqtype DNA

# Based on this, the forward Agler oligo lines up against kmer 257 (GAGGTGGAAGGCCTACGGGTCGTCAACTTC) and the reverse against kmer 553 (CGAAAGCACTCTGCTGGGCCGACACTGACA)
  # these have 335 and 758 hits in Greengenes, respectively. Their neighbors have slightly lower, but is going for a few greengenes hits worth changing the original protocol?
  
  
# # Line up the mitochondria oligos
# cat $mitochondria 0_agler_mito_blocking_olgos.fa > $workdir/1e_seqs_to_align.mito.fa
# clustalo -i $workdir/1e_seqs_to_align.mito.fa -o $workdir/1e_agler_seqs_aligned_to_maize.mito.clu --threads=$maxprocs --verbose --outfmt clu --force --output-order tree-order --seqtype DNA

# TODO: Mitochondria are turning out to be much trickier. Far more divergence between arabidopsis and maize. Is it targeting ITS regions? Also amplicon seems large, ~600 bp. Where do the universal primers sit on this?