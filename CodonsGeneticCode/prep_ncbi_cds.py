#!/usr/bin/env python3

import shutil, sys
from pathlib import Path
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC


"""Script is intended to tighten up the names from the "cds_from_genomic.fna" files
found on NCBI's RefSeq/GenBank datasets. This is a companion script to cbias.py and
strd_pol_smmry.py, but has its own uses as well."""

def rename_cds(fasta_file, taxons):
    orig_new = {}
    with open(f'{taxon}.WGS.CDS.fas','w+') as w:
        for i in SeqIO.parse(fasta_file,'fasta'):
            # Likely just need locus_tag, but being extra cautious.
            if ('protein_id' in i.description and
                '[psuedo=true]' not in i.description and
                'locus_tag' in i.description):
                g_id = i.description.split('locus_tag=')[1].split(']')[0]
                orig_new[i.id] = g_id
                w.write(f'>{g_id}\n{i.seq}\n')

    with open(f'{taxon}.WGS.CDS.Naming.tsv','w+') as w:
        w.write('CDS-ID\tBrief-Name\n')
        for k, v in orig_new.items():
            w.write(f'>{k}\t{v}\n')

if __name__ == '__main__':
    if len(sys.argv[1:]) == 2:
        fasta_file = sys.argv[1]
        taxon = sys.argv[2]
    else:
        print('Usage:\n\n    python3 prep_ncbi_cds.py [FASTA-FILE] [TAXON]\n')
        sys.exit(1)

    rename_cds(fasta_file, taxon)
