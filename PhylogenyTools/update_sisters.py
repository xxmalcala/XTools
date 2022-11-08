#!/usr/bin/env python3
import sys
import numpy as np
from collections import defaultdict

"""
Script 'simplifies' the output of the phylogeny walking "sisters" script
(as packaged with the phylogenomic pipeline, PhyloToL).

Specifically, this script ignores all the columns for each taxon-code (columns 15+),
as well as merges near identical taxa (e.g. "bacterial bins";
Am_ab_Mabd + Am_ar_Mabd --> Am_a_Mabd) into a single row.

NOTE: This script ASSUMES the output table is using 10-character taxon codes, which
are part of PhyloToL's methodology (current version: v4.3)

Last-Updated: 2022-11-08
Original-Author: Xyrus [xmaurer-alcala@amnh.org]
Most-Recent-Updat(er): Xyrus [see above]
"""

def fix_sisters(sister_csv):
    sisters_out_tsv = sister_csv.replace('csv','Cleaned.tsv')
    # Make a dictionary where keys have an "empty" array.
    sisters = defaultdict(lambda: np.zeros(11))

    # Parse the relevant columns of the sisters-summary file.
    sister_lines = [i.split(',')[:13] for i in open(sister_csv).readlines()]
    header = '\t'.join(sister_lines[0][:10]+sister_lines[0][11:])

    for line in sister_lines[1:]:
        # Simplify the taxon name.
        taxon = f'{line[0][:4]}{line[0][5:]}'
        # Adding the tip counts to the simplified taxon name.
        sisters[taxon] += np.array(line[1:10]+line[11:]).astype(int)

    with open(sisters_out_tsv, 'w+') as w:
        w.write(f'{header}\tProportion-Same\n')
        for k, v in sisters.items():
            # Calculate the "proportion-same" for each taxon.
            ps = f'{v[0]/v[-1]:.3f}'
            str_counts = '\t'.join([f'{i}' for i in v])
            w.write(f'{k}\t{str_counts}\t{ps}\n')

if __name__ == '__main__':
    try:
        sister_csv = sys.argv[1]
    except IndexError:
        print('Usage:\n\n    python3 update_sisters.py [SISTER-SUMMARY-CSV]\n')
        sys.exit(1)

    fix_sisters(sister_csv)
