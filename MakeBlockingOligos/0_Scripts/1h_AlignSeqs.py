__author__ = 'jgwall'

import argparse
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt

debug = False

def main():
    args = parse_args()

    # Load data
    target = list(SeqIO.parse(args.targetfile, "fasta"))
    queries = [seq for seq in SeqIO.parse(args.queryfile, "fasta")]

    # Check target length and trim down
    if len(target) != 1:
        print("Warning! Target file has",len(target),"sequences in it but only the first will be used")
    target = target[0]
    print("Target has length",len(target),"and will align",len(queries),"queries against it")

    # Align
    alignments = [align_seqs(target, q) for q in queries]

    # Sort if requested
    if args.sort_alignments:
        def sortalign(x): return x['start']
        alignments = sorted(alignments, key=sortalign)

    # Write out
    OUT = open(args.outfile, "w")
    for a in alignments:
        # Print out to file
        name = a['name']
        if a['revcomp']: name += "(reverse_complement)"
        OUT.write(name + "\n")
        OUT.write(format_alignment(*a['raw']) + "\n")
    OUT.close()

    # Print graphic
    fig = plt.figure()
    ax = fig.add_axes((0.1, 0.1, 0.65, 0.8), title="Oligo lineup")

    # Plot target
    # ax.hlines(xmin=0, xmax=len(target), y=-1, linestyles='solid', linewidth=3, label=target.id)
    arrow_length = len(target)/100
    ax.arrow(x=0, y=-1, dx=len(target), dy=0, linestyle='solid', linewidth=3, label=target.id, head_width=0.2, head_length=arrow_length)

    # Plot queries
    for i in range(len(alignments)):
        a = alignments[i]
        # ax.hlines(xmin=a['start'], xmax=a['stop'], y=i, linestyles='solid', linewidth=2, label=a['name'])
        xdist = (a['stop']-a['start'])
        if a['revcomp']: xdist *= -1
        color = 'red' if a['revcomp'] else 'blue'
        ax.arrow(x=a['start'], y=i, dx=xdist, dy=0, linestyle='solid', linewidth=3, label=target.id, head_width=0.2,
                 head_length=arrow_length, color=color)
        ax.axhline(y=i, zorder=-20, linestyle='dashed', linewidth=0.5, color='gray')

    #Add legend and format
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='x-small')

    ax.set_xlim(left=arrow_length, right=len(target)+arrow_length)
    ax.set_ylim(bottom=-2, top=len(alignments)+1)
    ax.invert_yaxis()

    #
    yticks = [-1] +  list(range(len(alignments)))
    ylabels = [target.id] + [a['name'] for a in alignments]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize="x-small")
    ax.yaxis.set_ticks_position('right')  # left, right, both, or none
    ax.tick_params(labeltop='off', labelright='on', labelleft='off')


    fig.savefig(args.outfile + ".png", dpi=100)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--targetfile", help="Sequence to align other sequences to")
    parser.add_argument("-q", "--queryfile", help="Sequence(s) to align against target sequence")
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("-s", "--sort-alignments", default=False, action="store_true", help="Whether to sort alignments")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def align_seqs(target, query):
    targetID, queryID, scoreID, startID, stopID = 0, 1, 2, 3, 4  # Parts of the alginment tuple
    target_gapopen, query_gapopen = -100, -1
    target_gapextend, query_gapextend = -100, -0.2
    print("\tAligning",query.id)
    # forward = pairwise2.align.localms(target.seq, query.seq, 2, -1, gapopen, gapextend)
    # revcomp = pairwise2.align.localms(target.seq, query.seq.reverse_complement(), 2, -1, gapopen, gapextend)

    forward = pairwise2.align.localmd(target.seq, query.seq, 2, -1, target_gapopen, query_gapopen, target_gapextend, query_gapextend)
    revcomp = pairwise2.align.localmd(target.seq, query.seq.reverse_complement(), 2, -1, target_gapopen, query_gapopen, target_gapextend, query_gapextend)

    # Find best alignment
    best = forward[0]
    for a in forward[1:]:
        print("\t\tAlignment has score",a[scoreID],"versus",best[scoreID])
        if a[scoreID] > best[scoreID]:
            print("\t\t\tNew best located")
            best = a

    # Check reverse complement
    reversed=False
    for r in revcomp:
        print("\t\tRevcomp alignment has score",r[scoreID],"versus",best[scoreID])
        if r[scoreID] > best[scoreID]:
            print("\t\t\tNew best located")
            best = r
            reversed=True


    # Make and return output as an easier to use dictonary which includes the original alignment
    output= { 'target':best[targetID], 'query':best[queryID], 'score':best[scoreID], 'start':best[startID], \
            'stop':best[stopID], 'revcomp':reversed, 'name':query.id, 'raw':best }
    return output

if __name__ == '__main__': main()