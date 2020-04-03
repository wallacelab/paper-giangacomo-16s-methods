__author__ = 'jgwall'

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import re

debug = False


def main():
    args = parse_args()
    print("Parsing primer counts")

    # If debug, trim to just first few files
    if debug:
        args.infiles = args.infiles[:10]

    # Load files
    data=dict()
    for i in args.infiles:
        mydata = pd.read_csv(i, header=None, sep=' ', skipinitialspace=True, names=['count','primer'], index_col='primer')
        mysample = re.sub(pattern=".+1a_(.+)_R[12].collated.txt", repl="\\1", string=i)
        mydata = mydata.rename(index=str).rename(columns={"count":mysample}, index={"nan":"_None_"})  # Have to convert NaN index to string before replacing

        # Add samples together if one half already loaded
        if mysample in data:
            data[mysample] = data[mysample] + mydata
        else:
            data[mysample] = mydata

    # Concatenate into a single DataFrame
    data = pd.concat(data.values(), axis=1, sort=True)
    data = data.fillna(value=0)

    # Write out concatenated frame
    print("Saving data frame")
    data.to_csv(args.outprefix + ".combined.tsv", sep='\t')

    # Scale each column to 0-1
    data_norm = data.copy()
    for mycol in data_norm.columns:
        data_norm[mycol] = data_norm[mycol] / sum(data_norm[mycol])

    # Plot heatmap
    print("Plotting heatmap")
    sns.heatmap(data_norm, linewidth=0.1)
    plt.ylim([-0.5,len(data_norm)])
    plt.xticks(ticks = np.arange(len(data_norm.columns)) + 0.5, labels=data_norm.columns, size='x-small')
    plt.yticks(size='x-small')
    plt.gca().invert_yaxis()

    # Save figure
    fig = plt.gcf()
    fig.set_size_inches(w=0.25 * len(data_norm.columns) + 2, h=5)
    plt.savefig(args.outprefix + ".heatmap.png", dpi=150)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="List of processed cutadapt output files")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()
