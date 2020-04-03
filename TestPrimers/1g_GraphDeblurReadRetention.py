__author__ = 'jgwall'

import argparse
import matplotlib.pyplot as plt
import pandas as pd

debug = False


def main():
    args = parse_args()

    # Load data
    print("Plotting Deblur mapping percentages\n")
    stats = pd.read_csv(args.infile)

    # Make graphic
    ax = plt.gca()
    xvals=range(len(stats))
    ax.bar(x=xvals, height=stats['reads-raw'], color="salmon", tick_label=stats['sample-id'], label="Unmapped reads")
    ax.bar(x=xvals, height=stats['reads-hit-reference'], color="blue", tick_label=stats['sample-id'], label="Mapped reads")
    ax.legend()

    plt.savefig(args.outfile, dpi=300)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Stats CSV file exported from deblur step of QIIME2")
    parser.add_argument("-o", "--outfile", help="Output graphic file")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()