#! /bin/bash

# Make blocking oligos according to Agler et al 2016

# TODO: Need to clean this pipeline up. A lot of cruft accumulating

PRIMER3=/home/jgwall/Software/Sequencing/primer3-2.3.7/src/primer3_core

# Maize chloroplast v3/v4 region from preliminary data (random chloroplast OTU)
# >1131894 M115_12346
# GTGTCAGCCGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTGGAGTACGGTAGGGGCAGAGGGAATTTCCAGTGGAGCGGTGAAATGCGTAGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGGATTAGAAACCCTCGTAGTCC

# Corresponding region cut out of the maize AGPv4 plastid genome
# >maize_plastid
# GTGCCAGCAGCCGCGGTAAGACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTCAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTGGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCATTGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCAAATGGGATTAGAGACCCCAGTAGTCC

# Same, but longer. The Agler et al equivalen from figure S2 seems to start at ~nt 200 in this example
# >maize_plastid
# CTTGGGAGGGGAACAACAACTGGAAACGGTTGCTAATACCCCGTAGGCTGAGGAGCAAAAGGAGAAATCCGCCCAAGGAGGGGCTCGCGTCTGATTAGCTAGTTGGTGAGGCAATAGCTTACCAAGGCGATGATCAGTAGCTGGTCCGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCGCAATGGGCGAAAGCCTGACGGAGCAATGCCGCGTGGAGGTGGAAGGCCTACGGGTCGTCAACTTCTTTTCTCGGAGAAGAAACAATGACGGTATCTGAGGAATAAGCATCGGCTAACTCTGTGCCAGCAGCCGCGGTAAGACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTCAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTGGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCATTGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCAAATGGGATTAGAGACCCCAGTAGTCC


# TODO: Note: Agler et al actually went back ~200 nt from the EMP, so pull a longer sequence

# Directory setup
scriptdir=0_Scripts
seqdir=0_Sequences
blastdir=1_BLAST
if [ ! -e $blastdir ]; then mkdir $blastdir; fi

# file locations
plastid=$seqdir/maize_plastid_16s.fa
mitochondria=$seqdir/maize_mito_16s.fa
primers=$seqdir/16s_primers.fa
greengenes=/usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta
maxprocs=7

# # Step 1 - get 30mers from target (maize) chloroplast region
# python3 $scriptdir/1_DivideTargetIntoKmers.py -i $plastid -o $blastdir/1_target_kmers.fa -k 30

# # Step 2 - BLAST against Greengenes
# # Make BLAST database
# makeblastdb -in $greengenes -title Greengenes -out $blastdir/1a_greengenes -dbtype nucl   
# # Perform BLAST
# blastn -query $blastdir/1_target_kmers.fa -db $blastdir/1a_greengenes -out $blastdir/1b_kmer_hits.txt -perc_identity 50 -word_size 12 -dust no -evalue 1000 -max_target_seqs 10000 -num_threads $maxprocs -outfmt 6
# pigz $blastdir/1b_kmer_hits.txt

# # Step 3 - Graph hits
# python3 $scriptdir/1c_TallyKmerHits.py -i $blastdir/1b_kmer_hits.txt.gz -o $blastdir/1c_kmer_hits.tally

# Step 4 - Determine best oligos to use - Need Tm, hairpin likelihood, and location
# Figure out where Agler's fall in maize plastid and how close a match they are
# cat $plastid 0_agler_plastid_blocking_olgos.fa > $blastdir/1d_seqs_to_align.fa
# clustalo -i $blastdir/1d_seqs_to_align.fa -o $blastdir/1d_agler_seqs_aligned_to_maize.clu --threads=$maxprocs --verbose --outfmt clu --force --output-order tree-order --seqtype DNA

# Based on this, the forward Agler oligo lines up against kmer 257 (GAGGTGGAAGGCCTACGGGTCGTCAACTTC) and the reverse against kmer 553 (CGAAAGCACTCTGCTGGGCCGACACTGACA)
  # these have 335 and 758 hits in Greengenes, respectively. Their neighbors have slightly lower, but is going for a few greengenes hits worth changing the original protocol?
  
  
# # Line up the mitochondria oligos
# cat $mitochondria 0_agler_mito_blocking_olgos.fa > $blastdir/1e_seqs_to_align.mito.fa
# clustalo -i $blastdir/1e_seqs_to_align.mito.fa -o $blastdir/1e_agler_seqs_aligned_to_maize.mito.clu --threads=$maxprocs --verbose --outfmt clu --force --output-order tree-order --seqtype DNA

# TODO: Mitochondria are turning out to be much trickier. Far more divergence between arabidopsis and maize. Is it targeting ITS regions? Also amplicon seems large, ~600 bp. Where do the universal primers sit on this?

####
# Mitochondria
####

# # Step 1 - get 30mers from target (maize) mitochondria region
# python3 $scriptdir/1_DivideTargetIntoKmers.py -i $mitochondria -o $blastdir/1f_mito_kmers.fa -k 30

# Step 2 - BLAST against Greengenes
# blastn -query $blastdir/1f_mito_kmers.fa -db $blastdir/1a_greengenes -out $blastdir/1g_mito_hits.txt -perc_identity 50 -word_size 12 -dust no -evalue 1000 -max_target_seqs 10000 -num_threads $maxprocs -outfmt 6
# pigz $blastdir/1g_mito_hits.txt

# # # Step 3 - Graph hits
# python3 $scriptdir/1c_TallyKmerHits.py -i $blastdir/1g_mito_hits.txt.gz -o $blastdir/1g_mito_hits.tally

# # # # # # # Get locations of primers and blocking oligos
# # # # # # # cat $mitochondria $primers $seqdir/agler_mito_blocking_olgos.fa > $blastdir/1h_mito_seqs_to_align.fa
# # # # # # 
# # # # # # # $mitochondria 
# # # # # # # clustalo -i $blastdir/1h_mito_seqs_to_align.fa -o $blastdir/1h_mito_seqs.clu --threads=$maxprocs --verbose --outfmt clu --force --output-order tree-order --seqtype DNA


####
# Determine the best alignment for everything
####

# # Find best alignment for each thing
# cat $primers $seqdir/agler_mito_blocking_olgos.fa $seqdir/agler_plastid_blocking_olgos.fa $seqdir/pna_seqs.fa > $blastdir/1h_seqs_to_align.fa

# # Test Arabidopsis stuff
# python3 $scriptdir/1h_AlignSeqs.py -t $seqdir/arabidopsis_chloro_16s.fa -q $blastdir/1h_seqs_to_align.fa -o $blastdir/1j_arabidopsis_chloro.txt --sort-alignment #--debug
# python3 $scriptdir/1h_AlignSeqs.py -t $seqdir/arabidopsis_mito_16s.fa -q $blastdir/1h_seqs_to_align.fa -o $blastdir/1j_arabidopsis_mito.txt --sort-alignment #--debug


# # Now maize again
# python3 $scriptdir/1h_AlignSeqs.py -t $seqdir/maize_chloro_16s.fa -q $blastdir/1h_seqs_to_align.fa -o $blastdir/1j_maize_chloro.txt --sort-alignment #--debug
# python3 $scriptdir/1h_AlignSeqs.py -t $seqdir/maize_mito_16s.fa -q $blastdir/1h_seqs_to_align.fa -o $blastdir/1j_maize_mito.txt --sort-alignment #--debug


####
# Mitochondria redux
####

# Okay, now that I have a clearer idae of how everything fits together, time to design some actual blocking oligos for the Maize 799-1192 mitochondria region
mito_target=$seqdir/maize_mito_target_region.fa
# # Step 1 - get 30mers from target (maize) mitochondria region
# python3 $scriptdir/1_DivideTargetIntoKmers.py -i $mito_target -o $blastdir/1k_mito_target_kmers.fa -k 30

# # Step 2 - BLAST against Greengenes
# blastn -query $blastdir/1k_mito_target_kmers.fa -db $blastdir/1a_greengenes -out $blastdir/1l_mito_target_hits.txt -perc_identity 50 -word_size 12 -dust no -evalue 1000 -max_target_seqs 10000 -num_threads $maxprocs -outfmt 6
# pigz $blastdir/1l_mito_target_hits.txt

# # # # Step 3 - Graph hits
# python3 $scriptdir/1c_TallyKmerHits.py -i $blastdir/1l_mito_target_hits.txt.gz -k $blastdir/1k_mito_target_kmers.fa -o $blastdir/1m_mito_target_hits.tally


# Went in to file manually to pull sequences that looked good. Time to check if they work
cat $primers $seqdir/agler_mito_blocking_olgos.fa $seqdir/agler_plastid_blocking_olgos.fa $seqdir/pna_seqs.fa $seqdir/wallace_mito_blocking_oligos.maize.fa > $blastdir/1n_seqs_to_align.fa
python3 $scriptdir/1h_AlignSeqs.py -t $seqdir/maize_mito_16s.fa -q $blastdir/1n_seqs_to_align.fa -o $blastdir/1n_maize_mito.txt --sort-alignment #--debug
