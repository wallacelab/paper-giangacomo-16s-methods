__author__ = 'jgwall'

import argparse
import gzip
import matplotlib.pyplot as plt
from Bio import SeqIO
from primer3 import calcTm, calcHairpin, calcHomodimer

debug = False


def main():
    args = parse_args()
    print("Tabulating BLAST hits in",args.infile)

    # Count hits
    hits = dict()
    kmer_order = list()  # Save the order these were all loaded in
    for line in get_filehandle(args.infile, "rt"):
        query=line.strip().split('\t')[0]
        if query not in hits:
            hits[query]=0
            kmer_order.append(query)
        hits[query]+=1
    print("\tFound",len(hits),"queries to tally")

    # Load additional characteristics of each kmer if provided sequence data
    characteristics = {kmer:"" for kmer in hits}    # Initialize so that if it's not loaded, it doesn't affect the output string.
    if args.kmer_seqs:
        characteristics = calc_characteristics(args.kmer_seqs)

    # Make text output
    outtext = args.outprefix + ".txt"
    print("Writing kmer counts to",outtext)
    OUT = open(outtext, "wt")
    if args.kmer_seqs:
        OUT.write("kmer\thits\tseq\ttm\ttm_hairpin\ttm_homodimer\tseq_rc\ttm_hairpin_rc\ttm_homodimer_rc\n")
    else:
        OUT.write("kmer\thits\n")
    for kmer in kmer_order:
        OUT.write(kmer + "\t" + str(hits[kmer]) + characteristics[kmer] + "\n")
    OUT.close()

    # Make graphical output
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111, ylabel="Number of BLAST hits to each kmer")
    hit_count = [hits[k] for k in kmer_order]
    xvals = list(range(len(hit_count)))
    ax.bar(left=xvals, height=hit_count, linewidth=0)

    # Prettify axes
    pad = len(xvals)/10
    ax.set_xlim(min(xvals)-pad, max(xvals)+pad)

    # Output graphic
    outpng = args.outprefix + ".png"
    print("Making graphical output in", outpng)
    fig.savefig(outpng)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Tabular-formatted BLAST output file with the query sequence in the first column")
    parser.add_argument("-o", "--outprefix", help="Output prefix for output files")
    parser.add_argument("-k", "--kmer-seqs", help="FASTA file of all the kmer sequences")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def calc_characteristics(infile):
    print("Using sequence file",infile,"to calculate Tm and structure potential")
    seqs = list(SeqIO.parse(infile, "fasta"))
    chars = dict()
    for myseq in seqs:
        # Forward
        forward = str(myseq.seq)
        hairpinF = calcHairpin(forward)
        homoF = calcHomodimer(forward)

        #Reverse complement
        reverse = str(myseq.seq.reverse_complement())
        hairpinR = calcHairpin(reverse)
        homoR = calcHomodimer(reverse)
        # print(hairpinF,'\n', homoF,'\n', hairpinR,'\n', homoR, '\n', calcTm(forward))

        # Values to save
        tm = calcTm(forward)
        hairpinF = hairpinF.tm if hairpinF.structure_found else "NA"
        hairpinR = hairpinR.tm if hairpinR.structure_found else "NA"
        homoF = homoF.tm if homoF.structure_found else "NA"
        homoR = homoR.tm if homoR.structure_found else "NA"

        output = [str(x) for x in [forward, tm, hairpinF, homoF, reverse, hairpinR, homoR]]
        chars[myseq.id] = "\t" + "\t".join(output)
        # print(chars[myseq.id])
    return chars



def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    return open(file, mode)

if __name__ == '__main__': main()