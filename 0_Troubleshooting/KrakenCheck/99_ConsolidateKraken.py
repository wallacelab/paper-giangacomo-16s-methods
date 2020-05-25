__author__ = 'jgwall'

import argparse
import re
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Consolidating Kraken results for",len(args.infiles),"input files")
    results = [load_summaries(i) for i in args.infiles]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="List of Kraken summary files")
    parser.add_argument("-o", "--outfile", help="Output file of tallied results")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_summaries(infile):
    # Get sample name
    name = re.search("\.(.+)\.txt$", string=infile).group(1)
    print(name)


if __name__ == '__main__': main()