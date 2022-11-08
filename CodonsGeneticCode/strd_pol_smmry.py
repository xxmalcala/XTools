#!/usr/bin/env python3

import shutil, sys
from collections import defaultdict
from itertools import groupby
from operator import itemgetter

"""Use to aid in identification of absolute strand polarity and possible
polycistronic mRNAs from a given genome's GTF file."""

# Parses the GTF file and returns a dictionary with the scaffold/chromosome (key)
# and a list of tuples with relevant information (value; scaffold name, gene-id, strand).
# Note: this only looks at protein-coding genes, excluding psuedogenized genes.
def parse_gtf(gtf_file):
    gtf_data = defaultdict(list)
    for line in open(gtf_file).readlines():
        if (line[0] != '#' and
            line.split('\t')[2] == 'gene' and
            'pseudo' not in line):

            scf = line.split('\t')[0]
            strd = line.split('\t')[6]
            gene_id = line.split('gene_id "')[1].split('"')[0]

            gtf_data[scf].append((scf, gene_id, strd))

    return gtf_data

# Summarizes the output of the gtf-parsing. Also identifies "groups" of consecutive
# genes on the same strand.

def strand_summarize(gtf_data, taxon):
    with open(f'{taxon}.PCG_Strand_Summary.tsv','w+') as w:
        w.write('Genomic-Scaffold\tGene-ID\tStrand\tGroup-Number\n')
        for k, v in gtf_data.items():
            if v:
                grp_n = 1
                strd_smmry = [list(group) for key, group in groupby(v, itemgetter(2))]
                for i in strd_smmry:
                    for j in i:
                        w.write(f'{k}\t{j[1]}\t{j[2]}\t{grp_n}\n')
                    grp_n += 1

if __name__ == '__main__':
    if len(sys.argv[1:]) == 2:
        gtf_file = sys.argv[1]
        taxon = sys.argv[2]
    else:
        print('Usage:\n\n    python3 strd_pol_smmry.py [GTF-FILE] [TAXON]\n')
        sys.exit(1)

    gtfd = parse_gtf(gtf_file)
    strand_summarize(gtfd, taxon)
