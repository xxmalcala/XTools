#!/usr/bin/env python3

import shutil, sys
from pathlib import Path
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

def prep_folders(taxon):
    Path(f'{taxon}_CodonUsageBias/SpreadSheets').mkdir(parents=True, exist_ok=True)
    Path(f'{taxon}_CodonUsageBias/Original').mkdir(parents=True, exist_ok=True)
    return f'{taxon}_CodonUsageBias'


def getCDNtable(gCode):
    # Returns the appropriate codon table to be used for the ENc calculations.
    # Universal codon table, with 6-fold degenerate codons split
    # into four-fold and two-fold groups.
    universal_no6fold = {
    	'GCT': ['A', 'four', 0], 'GCC': ['A', 'four', 0], 'GCA': ['A', 'four', 0],
    	'GCG': ['A', 'four', 0], 'CGT': ['R', 'four', 0], 'CGC': ['R', 'four', 0],
    	'CGG': ['R', 'four', 0], 'CGA': ['R', 'four', 0], 'AGA': ['R_', 'two', 0],
    	'AGG': ['R_', 'two', 0], 'AAT': ['N', 'two', 0], 'AAC': ['N', 'two', 0],
    	'GAT': ['D', 'two', 0], 'GAC': ['D', 'two', 0], 'TGT': ['C', 'two', 0],
    	'TGC': ['C', 'two', 0], 'CAA': ['Q', 'two', 0], 'CAG': ['Q', 'two', 0],
    	'GAA': ['E', 'two', 0], 'GAG': ['E', 'two', 0], 'GGT': ['G', 'four', 0],
    	'GGC': ['G', 'four', 0], 'GGA': ['G', 'four', 0], 'GGG': ['G', 'four', 0],
    	'CAT': ['H', 'two', 0], 'CAC': ['H', 'two', 0], 'ATT': ['I', 'three', 0],
    	'ATC': ['I', 'three', 0], 'ATA': ['I', 'three', 0], 'ATG': ['M', 'one', 0],
    	'TTA': ['L_', 'two', 0], 'TTG': ['L_', 'two', 0], 'CTT': ['L', 'four', 0],
    	'CTC': ['L', 'four', 0], 'CTA': ['L', 'four', 0], 'CTG': ['L', 'four', 0],
    	'AAA': ['K', 'two', 0], 'AAG': ['K', 'two', 0], 'TTT': ['F', 'two', 0],
    	'TTC': ['F', 'two', 0], 'CCT': ['P', 'four', 0], 'CCC': ['P', 'four', 0],
    	'CCA': ['P', 'four', 0], 'CCG': ['P', 'four', 0], 'TCT': ['S', 'four', 0],
    	'TCC': ['S', 'four', 0], 'TCA': ['S', 'four', 0], 'TCG': ['S', 'four', 0],
    	'AGT': ['S_', 'two', 0], 'AGC': ['S_', 'two', 0], 'ACT': ['T', 'four', 0],
    	'ACC': ['T', 'four', 0], 'ACA': ['T', 'four', 0], 'ACG': ['T', 'four', 0],
    	'TGG': ['W', 'one', 0], 'TAT': ['Y', 'two', 0], 'TAC': ['Y', 'two', 0],
    	'GTT': ['V', 'four', 0], 'GTC': ['V', 'four', 0], 'GTA': ['V', 'four', 0],
    	'GTG': ['V', 'four', 0], 'TAA': ['*', 'none', 0], 'TGA': ['*', 'none', 0],
    	'TAG': ['*', 'none', 0], 'XXX': ['_missing', 'none', 0]}
    # Universal codon table, with 6-fold degenerate codons kept
    # whole, no splitting! Traditional Universal codon table.
    universal_6fold = {
    	'GCT': ['A', 'four', 0], 'GCC': ['A', 'four', 0], 'GCA': ['A', 'four', 0],
    	'GCG': ['A', 'four', 0], 'CGT': ['R', 'six', 0], 'CGC': ['R', 'six', 0],
    	'CGG': ['R', 'six', 0], 'CGA': ['R', 'six', 0], 'AGA': ['R', 'six', 0],
    	'AGG': ['R', 'six', 0], 'AAT': ['N', 'two', 0], 'AAC': ['N', 'two', 0],
    	'GAT': ['D', 'two', 0], 'GAC': ['D', 'two', 0], 'TGT': ['C', 'two', 0],
    	'TGC': ['C', 'two', 0], 'CAA': ['Q', 'two', 0], 'CAG': ['Q', 'two', 0],
    	'GAA': ['E', 'two', 0], 'GAG': ['E', 'two', 0], 'GGT': ['G', 'four', 0],
    	'GGC': ['G', 'four', 0], 'GGA': ['G', 'four', 0], 'GGG': ['G', 'four', 0],
    	'CAT': ['H', 'two', 0], 'CAC': ['H', 'two', 0], 'ATT': ['I', 'three', 0],
    	'ATC': ['I', 'three', 0], 'ATA': ['I', 'three', 0], 'ATG': ['M', 'one', 0],
    	'TTA': ['L', 'six', 0], 'TTG': ['L', 'six', 0], 'CTT': ['L', 'six', 0],
    	'CTC': ['L', 'six', 0], 'CTA': ['L', 'six', 0], 'CTG': ['L', 'six', 0],
    	'AAA': ['K', 'two', 0], 'AAG': ['K', 'two', 0], 'TTT': ['F', 'two', 0],
    	'TTC': ['F', 'two', 0], 'CCT': ['P', 'four', 0], 'CCC': ['P', 'four', 0],
    	'CCA': ['P', 'four', 0], 'CCG': ['P', 'four', 0], 'TCT': ['S', 'six', 0],
    	'TCC': ['S', 'six', 0], 'TCA': ['S', 'six', 0], 'TCG': ['S', 'six', 0],
    	'AGT': ['S', 'six', 0], 'AGC': ['S', 'six', 0], 'ACT': ['T', 'four', 0],
    	'ACC': ['T', 'four', 0], 'ACA': ['T', 'four', 0], 'ACG': ['T', 'four', 0],
    	'TGG': ['W', 'one', 0], 'TAT': ['Y', 'two', 0], 'TAC': ['Y', 'two', 0],
    	'GTT': ['V', 'four', 0], 'GTC': ['V', 'four', 0], 'GTA': ['V', 'four', 0],
    	'GTG': ['V', 'four', 0], 'TAA': ['*', 'none', 0], 'TGA': ['*', 'none', 0],
    	'TAG': ['*', 'none', 0], 'XXX': ['_missing', 'none', 0]}
    # Blepharisma (table 4) genetic code codon table, with 6-fold degenerate
    # codons kept whole, no splitting!
    blepharisma_6fold = {**universal_6fold,
        'TGA': ['W', 'two', 0], 'TGG': ['W', 'two', 0],
        'TAA': ['*', 'two', 0], 'TAG': ['*', 'two', 0]}
    # Blepharisma (table 4) genetic code codon table, with 6-fold degenerate
    # codons split into four-fold and two-fold groups.
    blepharisma_no6fold = {**universal_no6fold,
        'TGA': ['W', 'two', 0], 'TGG': ['W', 'two', 0],
        'TAA': ['*', 'two', 0], 'TAG': ['*', 'two', 0]}
    # Chilodonella genetic code codon table, with 6-fold degenerate
    # codons kept whole, no splitting!
    chilo_6fold = {**universal_6fold,
        'CAA': ['Q', 'four', 0], 'CAG': ['Q', 'four', 0],
        'TAA': ['*', 'one', 0], 'TAG': ['Q', 'four', 0],
        'TGA': ['Q', 'four', 0]}
    # Chilodonella genetic code codon table, with 6-fold degenerate
    # codons split into four-fold and two-fold groups.
    # Note that this also splits four-fold degenerate codons that OUGHT to
    # be in "different" functional categories (e.g. CAG =/= TAG)
    chilo_no6fold = {**universal_no6fold,
        'TAA': ['*', 'one', 0], 'TAG': ['Q_', 'one', 0],
        'TGA': ['Q_', 'one', 0]}
    # Ciliate (table 6) genetic code codon table, with 6-fold degenerate
    # codons kept whole, no splitting! Traditional ciliate codon table.
    ciliate_6fold = {**universal_6fold,
        'CAA': ['Q', 'four', 0], 'CAG': ['Q', 'four', 0],
        'TAA': ['Q', 'four', 0], 'TAG': ['Q', 'four', 0],
        'TGA': ['*', 'one', 0]}
    # Ciliate (table 6) genetic code codon table, with 6-fold degenerate
    # codons split into four-fold and two-fold groups.
    # Note that this also splits four-fold degenerate codons that OUGHT to
    # be in "different" functional categories (e.g. CAA =/= TAA)
    ciliate_no6fold = {**universal_no6fold,
        'TAA': ['Q_', 'two', 0], 'TAG': ['Q_', 'two', 0],
        'TGA': ['*', 'one', 0]}
    # Euplotes codon table, with 6-fold degenerate codons kept
    # whole, no splitting! Traditional Universal codon table.
    euplotes_6fold = {**universal_6fold,
        'TGA': ['C', 'three', 0], 'TGT': ['C', 'three', 0],
        'TGC': ['C', 'three', 0], 'TAA': ['*', 'two', 0],
        'TAG': ['*', 'two',0]}
    # Euplotes genetic code codon table, with 6-fold degenerate codons
    # split into four-fold and two-fold groups.
    euplotes_no6fold = {**universal_no6fold,
        'TGA': ['C', 'three', 0], 'TGT': ['C', 'three', 0],
        'TGC': ['C', 'three', 0], 'TAA': ['*', 'two', 0],
        'TAG': ['*', 'two',0]}
    # Mesodinium/Myrionecta (table 29) genetic code codon table, with 6-fold
    # degenerate codons kept whole, no splitting! Traditional ciliate codon table.
    mesodinium_6fold = {**universal_6fold,
        'TAA': ['Y', 'four', 0], 'TAT': ['Y', 'four', 0],
        'TAG': ['Y', 'four', 0], 'TAC': ['Y', 'four', 0],
        'TGA': ['*', 'one', 0]}
    # Mesodinium/Myrionecta (table 29) genetic code codon table, with 6-fold
    # degenerate codons split into four-fold and two-fold groups.
    mesodinium_no6fold = {**universal_no6fold,
        'TAA': ['Y', 'four', 0], 'TAT': ['Y', 'four', 0],
        'TAG': ['Y', 'four', 0], 'TAC': ['Y', 'four', 0],
        'TGA': ['*', 'one', 0]}
    # Peritrich (table 30) genetic code codon table, with 6-fold degenerate
    # codons kept whole, no splitting! Traditional ciliate codon table.
    peritrich_6fold = {**universal_6fold,
        'GAA': ['E', 'four', 0], 'GAG': ['E', 'four', 0],
        'TAA': ['E', 'four', 0], 'TAG': ['E', 'four', 0],
        'TGA': ['*', 'one', 0]}
    # Peritrich (table 30) genetic code codon table, with 6-fold degenerate
    # codons split into four-fold and two-fold groups.
    # Note that this also splits four-fold degenerate codons that OUGHT to
    # be in "different" functional categories (e.g. CAA =/= TAA)

    peritrich_no6fold = {**universal_no6fold,
        'TAA': ['E_', 'two', 0], 'TAG': ['E_', 'two', 0],
        'TGA': ['*', 'one', 0]}

    cdnTableDict = {1:[universal_no6fold,universal_6fold],
        4:[blepharisma_no6fold, blepharisma_6fold],
        6:[ciliate_no6fold,ciliate_6fold],
        10:[euplotes_no6fold,euplotes_6fold],
        29:[mesodinium_no6fold,mesodinium_6fold],
        30:[peritrich_no6fold,peritrich_6fold],
        'chilodonella':[chilo_no6fold,chilo_6fold],
        'chilo':[chilo_no6fold,chilo_6fold]}

    return cdnTableDict[gCode]

def mapCdns(seq, cdnTable):
    # Updates the codon counts for a given sequence to the respective codon
    # count table (e.g. with or without 6-fold degeneracy).
    codons = [seq[n:n+3] for n in range(0, len(seq)-len(seq)%3, 3)]
    amb_cdn = 0
    for c in codons:
        try:
            cdnTable[c][-1] += 1
        except:
            amb_cdn += 1

    if cdnTable['TCC'][1] == 'six':
        return cdnTable, amb_cdn

    else:
        return cdnTable

class GCeval():
    """
    Returns %GC values from DNA sequences of various types.
    """
    def gcTotal(seq):
        # This function returns global GC content
        return round(GC(seq), 4)

    def gc1(seq):
        # This function return the GC content of the first position of a codon
        return round(GC(''.join([seq[n] for n in range(0, len(seq), 3)])), 4)

    def gc2(seq):
        # This function return the GC content of the second position of a codon
        return round(GC(''.join([seq[n] for n in
            range(1, len(seq), 3)])), 4)

    def gc3(seq):
        # This function return the GC content of the third position of a codon
        return round(GC(''.join([seq[n] for n in
            range(2, len(seq[2:]), 3)])), 4)

    def gc3_4F(cdnTbl):
    #     # This function return the GC content of the third position of four-fold
    #     # degenerate codons
        FrFold = round(GC(''.join([k[-1]*v[-1] for k, v in cdnTbl.items() if
            'four' in v[1]])), 4)
        return FrFold

class CalcCUB:
    """
    Returns the Effective Number of Codons used (observed and expected)
    following the equations originally from Wright 1990.
    """
    def calc_exp_Wright_ENc(gc3):
        # Calculates the expected ENc from a sequence's GC3 under Wright 1990
        if gc3 > 1:
            # If GC3 looks as though it is > 1 (e.g. 100%), converts to a float ≤ 1.
            # Calculations expect a value between 0 and 1
            gc3 = gc3/100
        exp_enc = 2 + gc3 + (29 / ((gc3**2) + (1 - gc3)**2))
        return round(exp_enc, 4)
    def null_ENc_GC3():
        # Calculates the expected ENc from the null distribution of GC3
        # values (0, 100% GC)
        null = [CalcCUB.calc_exp_Wright_ENc(n) for n in np.arange(0,.51,0.01)]
        null += null[:-1][::-1]
        return [str(i)+'\t'+str(j) for i, j in zip([n for n in range(0, 101)],null)]

    def calc_Wright_ENc(cdnTable):

        # Follows Wright's (1990) calculations for determining ENc scores.
        def faCalcWright(aa_counts):
            # Returns the codon homozygosity (fa) for a given "type" of AA (e.g.
            # 2-fold degeneracy).
            counts = [i[2] for i in aa_counts]
            # n_aa --> number of this particular AA
            n_aa = sum(counts)
            # fa --> codon homozygosity
            try:
                fa = (((n_aa*sum([(i/float(n_aa))**2 for i in counts]))-1)/(n_aa-1))

            except:
                fa = 0

            return fa

        def ENcWright_by_Degen(fa_data):
            # Same as used in Wright 1990, averages the homozygosity across all codons
            # of a given class (e.g. 2-fold degeneracy)
            # Codons without any degeneracy (e.g. ATG == M) have 100% homozygosity
            # and provide a "base" for the ENc score
            enc = 2
            for k, v in fa_data.items():
                non_zero_vals, non_zero_sum = len([i for i in v if i != 0]), sum([i for i in v if i != 0])
                try:
                    f_aa = non_zero_sum/non_zero_vals
                except:
                    f_aa = 1
                enc += k/f_aa
            return enc

        # Determines the number of degenerate groups to use (i.e. whether 6-Fold
        # degeneracy is present).
        degen_cdns = {}
        for k, v in cdnTable.items():
            if v[1] not in degen_cdns.keys():
                degen_cdns[v[1]] = [v[0]]
            else:
                if v[0] not in degen_cdns[v[1]]:
                    degen_cdns[v[1]] += [v[0]]

        # Calculates codon homozygosity (fa) for each amino acid. Groups the
        # resulting values based on the amino acids degeneracy (e.g. 'two-fold').
        fa_cdns = {len(v):[] for k, v in degen_cdns.items() if 'one' not in k}
        for k, v in degen_cdns.items():
            # Skip codons lacking degeneracy
            if 'one' in k:
                continue
            for aa in v:
                aa_counts = [cdnTable[k] for k in cdnTable.keys() if cdnTable[k][0] == aa]
                fa_cdns[len(v)] += [faCalcWright(aa_counts)]

        enc_val = min(61, round(ENcWright_by_Degen(fa_cdns),4))
        return enc_val

def check_gcode(gcode):
    # supported_codes = {'universal':1, '1':1, 'blepharisma':4, '4':4, 'ciliate':6,
    #     '6':6, 'euplotes':10, '10':10, 'mesodinium':29, '29':29,
    #     'peritrich':30, '30':30, 'bacterial':11, '11':11, 'alt-yeast':12, '12':12}

    supported_codes = {'universal':1, '1':1, 'blepharisma':4, '4':4, 'ciliate':6,
        '6':6, 'euplotes':10, '10':10, 'mesodinium':29, '29':29,
        'peritrich':30, '30':30}

    if gcode not in supported_codes.keys():
        print('[ERROR]: Unsupported genetic code provided. Defaulting to Universal')
        return 1
    else:
        return supported_codes[gcode.lower()]


def CalcRefFasta(fasta_file, gcode):
    seqDB = {i.description:i for i in SeqIO.parse(fasta_file,'fasta')}
    GenCDNtable = {}
    for k, v in seqDB.items():
        codon_tables = getCDNtable(gcode)
        v.cdnCounts_6F, v.amb_cdn = mapCdns(f'{v.seq}', codon_tables[1])
        v.cdnCounts_No6F = mapCdns(v.seq, codon_tables[0])
        v.gc_total = GCeval.gcTotal(v.seq)
        v.gc1 = GCeval.gc1(v.seq)
        v.gc2 = GCeval.gc2(v.seq)
        v.gc3 = GCeval.gc3(v.seq)
        v.gc12 = (v.gc1+v.gc2)/2
        v.gc3_4F = GCeval.gc3_4F(v.cdnCounts_No6F)
        v.exp_ENc = CalcCUB.calc_exp_Wright_ENc(v.gc3_4F)
        v.obsENc_6F = CalcCUB.calc_Wright_ENc(v.cdnCounts_6F)
        v.obsENc_No6F = CalcCUB.calc_Wright_ENc(v.cdnCounts_No6F)
    return seqDB



if __name__ == '__main__':
    args = sys.argv[1:]
    if len(sys.argv[1:]) == 2:
        cds_fasta = sys.argv[1]
        taxon = sys.argv[2]
        gcode = 1
    elif len(sys.argv[1:]) == 3:
        cds_fasta = sys.argv[1]
        taxon = sys.argv[2]
        gcode = sys.argv[3]
    else:
        print('Usage:\n\n    python3 cbias.py [FASTA-FILE-CDS] [TAXON-NAME] ' \
            '[TRANSLATION-TABLE (default = 1)]\n')
        sys.exit(1)

    gcode_good = check_gcode(gcode)

    out_folder = prep_folders(taxon)
    shutil.copy2(cds_fasta, f'{out_folder}/Original/')
    seq_info = CalcRefFasta(cds_fasta, gcode_good)
    null_gc3 = CalcCUB.null_ENc_GC3()

    with open(f'{out_folder}/SpreadSheets/{taxon}.ENc_GC3.tsv','w+') as w:
        w.write(f'Gene\tLength\tTotal GC\tGC1\tGC2\tGC3\tGC12\tGC3-Four-Fold-Degen' \
            f'\tExpected Wright ENc\tObserved Wright 6-Fold ENc\t' \
            f'Observed Wright No 6-Fold ENc\n')

        for k, v in seq_info.items():
            w.write(f'{k}\t{len(v.seq)}\t{v.gc_total}\t{v.gc1}\t{v.gc2}\t{v.gc3}' \
                f'\t{v.gc12}\t{v.gc3_4F}\t{v.exp_ENc}\t{v.obsENc_6F}\t{v.obsENc_No6F}\n')

    with open(f'{out_folder}/SpreadSheets/{taxon}.Null_ENc_GC3.tsv','w+') as w:
        w.write('GC3\tENc\n')
        w.write('\n'.join(null_gc3))
