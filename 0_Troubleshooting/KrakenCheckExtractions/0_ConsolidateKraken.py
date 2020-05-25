__author__ = 'jgwall'

import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn

debug = False


def main():
    args = parse_args()

    # Load data
    print("Consolidating Kraken results for",len(args.infiles),"input files")
    results = [load_summaries(i) for i in args.infiles]
    results=pd.concat(results)

    # Write raw results
    out_all = args.outprefix + ".raw.txt"
    print("Writing combined results to", out_all)
    results.to_csv(out_all, sep='\t', index=None)

    # Filter to just focus level
    filtered = results.loc[results['tax_level']==args.filter_level, :]
    out_filt = args.outprefix + ".filt_long.txt"
    print("Writing filtered results to", out_filt)
    filtered.to_csv(out_filt, sep='\t', index=None)

    # Make wide format and combine with sample data
    key = pd.read_csv(args.keyfile, sep='\t', comment='#', index_col='sample-id')
    filtered_wide = filtered.pivot_table(index='sample', columns='tax_name', values='cum_percent', fill_value=0)
    # Add sample metadata
    filtered_wide['sample_type'] = key.loc[filtered_wide.index, 'sample-type']
    filtered_wide['method'] = key.loc[filtered_wide.index, 'treatment']
    filtered_wide = filtered_wide.set_index(keys=['sample_type', 'method'], drop=True, append=True)
    # Output
    out_wide = args.outprefix + ".filt_wide.txt"
    print("Writing filtered results in wide format to", out_wide)
    filtered_wide.to_csv(out_wide, sep='\t')

    # Plot heatmap
    plotdata = filtered_wide.sort_index(level=['sample_type', 'method'])
    plotdata = plotdata.loc[:, plotdata.max(axis=0) > 5]
    fig, ax = plt.subplots(figsize=[10,15])
    fig.subplots_adjust(left=0.3)
    seaborn.heatmap(plotdata, ax=ax)
    fig.savefig(args.outprefix + ".heatmap.png", dpi=150)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="List of Kraken summary files")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-f", "--filter-level", default='F', help="Taxonomic level to filter focused results to")
    parser.add_argument("-k", "--keyfile", help="QIIME keyfile")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_summaries(infile):
    # Get sample name
    name = re.search("\.(.+)\.txt$", string=infile).group(1)

    # Load data
    if debug: print("\tLoading",infile)
    data = pd.read_csv(infile, sep='\t', header=None, names=['cum_percent', 'cum_hits', 'exact_hits', 'tax_level', 'tax_id', 'tax_name'])
    data.insert(loc=0, column='sample', value=name) # Add sample name

    # Strip whitespace from taxon names
    data['tax_name'] = [t.strip() for t in data['tax_name']]

    return data

if __name__ == '__main__': main()