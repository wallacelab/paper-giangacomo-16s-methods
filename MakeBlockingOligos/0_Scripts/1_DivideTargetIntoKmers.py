__author__ = 'jgwall'

import argparse
import gzip
from Bio import SeqIO

debug = False


def main():
    args = parse_args()
    print("Cutting sequence in",args.infile,"into kmers of size",args.kmer_size)

    # Load sequence
    IN = get_filehandle(args.infile, "rt")
    seq = SeqIO.parse(IN, format="fasta")
    seq = list(seq)
    IN.close()

    # Make sure only one record
    if len(seq) >1:
        print("Error: Input file should have only 1 sequence but",len(seq),"were found")
    seq=str(seq[0].seq)

    # Slice into kmers
    kmers = list()
    for start in range(len(seq)-args.kmer_size):
        end = start + args.kmer_size
        kmers.append(seq[start:end])
    print(len(kmers),"kmers created from input sequence")

    # Write out
    print("Writing results to", args.outfile)
    OUT = get_filehandle(args.outfile, "wt")
    for k, myseq in zip(range(len(kmers)), kmers):
        OUT.write(">" + args.output_name + "_" + str(k) + "\n") # Name
        OUT.write(myseq + "\n")    # Kmer sequence
    OUT.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="FASTA-formatted input file of the target sequence (generally a chloroplast 16s region)")
    parser.add_argument("-o", "--outfile", help="Output file of FASTA-formatted kmers made from the input sequence")
    parser.add_argument("-k", "--kmer-size", type=int, default=30, help="Size of kmers to divide the target sequence into")
    parser.add_argument("-n", "--output-name", default="kmer", help="Name prefix for all output kmer sequences")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    return open(file, mode)

if __name__ == '__main__': main()