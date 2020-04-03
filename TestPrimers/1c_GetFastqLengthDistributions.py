__author__ = 'jgwall'

import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import re
import statistics

debug = False


def main():
    args = parse_args()
    print("Calculating sequence lengths in", len(args.infiles),"input files")

    lengths=dict()
    for infile in args.infiles:
        print("\tProcessing",infile)
        IN = get_filehandle(infile, "rt")
        sample = re.sub(string=infile, pattern=".+/", repl="")
        sample = re.sub(string=sample, pattern=".fastq.gz$", repl="")

        lengths[sample] = list()
        for seqname in IN:  # Sequence name (not important)
            sequence = IN.readline().strip()    # Sequence (important)
            IN.readline()   # "+" separator
            IN.readline()   # Quality scores
            lengths[sample].append(len(sequence))

            # Debugging; break once have 500 samples
            if debug and len(lengths[sample]) > 500:
                break


    # Process data
    samples = sorted(lengths.keys())

    # Write out raw data
    OUT = open(args.outprefix + ".raw_values.txt", "w")
    OUT.write("sample\tlength\n")
    for sample in samples:
        for mylength in lengths[sample]:
            OUT.write(sample + "\t" + str(mylength) + "\n")
    OUT.close()

    # Write out medians and quantiles per sample
    medians = [statistics.median(lengths[s]) for s in samples]
    q10 = [np.quantile(lengths[s], 0.1) for s in samples]   # 90% of samples larger than this
    q05 = [np.quantile(lengths[s], 0.05) for s in samples]  # 95% larger than this
    q01 = [np.quantile(lengths[s], 0.01) for s in samples]  # 99% larger than this
    
    
    OUT = open(args.outprefix + ".quantiles.txt", "w")
    OUT.write("sample\tmedian_length\tkeep_99\tkeep_95\tkeep_90\n")
    for mysample, mymedian, myq01, myq05, myq10 in zip(samples, medians, q01, q05, q10):
        values = [str(x) for x in [mysample, mymedian, myq01, myq05, myq10]]
        OUT.write("\t".join(values) + "\n")
    OUT.close()

    # Write out a graphic showing the distribution of values
    plotdata = [lengths[s] for s in samples]
    xticks = range(len(plotdata))
    plt.violinplot(plotdata, showmedians=True, positions=xticks)
    plt.xticks(ticks = xticks, labels=samples, rotation=90)
    plt.subplots_adjust(bottom=0.5)
    plt.savefig(args.outprefix + ".length_distributions.png", dpi=300)




def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="List of FASTQ files to get lengths from")
    parser.add_argument("-o", "--outprefix", help="Prefix to use for output files")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    else:
        return open(file, mode)

if __name__ == '__main__': main()
