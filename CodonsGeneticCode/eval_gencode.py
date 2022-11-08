#!/usr/bin/env python3

import os,subprocess, sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML

from collections import defaultdict
from pathlib import Path


def prepare_folders(taxon):
    output_folder = f'{taxon}_Eval_Genetic_Code/'
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    return output_folder


def dmnd_blastx_compare(fasta_file, taxon, out_folder, db_loc, evalue=1e-15, threads=2, num_hits=25):
    dmnd_cmd = f'diamond blastx ' \
        f'-p {threads} ' \
        f'-k {num_hits} ' \
        f'-e {evalue} ' \
        f'-f 5 ' \
        f'-q {fasta_file} ' \
        f'-d {db_loc} ' \
        f'-o {out_folder}{taxon}.GenCodeEval.Hits.xml'
    subprocess.call(dmnd_cmd, shell=True)
    return f'{out_folder}{taxon}.GenCodeEval.Hits.xml'


def get_codons(seq, start, end, frame):
    seq = seq[start:end]
    if frame < 0:
        seq = f'{Seq(seq).reverse_complement()}'
    codons = [seq[n:n+3] for n in range(0, len(seq),3)]
    return codons


def parse_fasta(fasta_file):
    nuc_seqs = {f'{i.id}':f'{i.seq.upper()}' for i in SeqIO.parse(fasta_file,'fasta')}
    return nuc_seqs


def parse_xml(xml_file, nuc_seqs):
    total_codons_eval = 0
    codon_assignments = defaultdict(list)

    for record in NCBIXML.parse(open(xml_file)):
        gene_name = record.query
        for aln in record.alignments:
            for hsp in aln.hsps:
                q_st = hsp.query_start-1
                q_end = hsp.query_end -1
                q_aln = hsp.query
                s_aln = hsp.sbjct
                true_pos = -1
                codons_by_pos = get_codons(nuc_seqs[gene_name], q_st, q_end, hsp.frame[0])
                for n in range(len(q_aln)):
                    if q_aln[n] != '-':
                        true_pos += 1
                        total_codons_eval += 1
                        if s_aln[n] != '-':
                            codon_assignments[codons_by_pos[true_pos]].append(s_aln[n])

    return codon_assignments, total_codons_eval


def keep_strong(codon_assignments, total_codons_eval, proportion = 0.1):
    strong_assignments = defaultdict(list)
    best_cdn_guess = {}
    for k, v in codon_assignments.items():
        for aa in list(set(v)):
            if v.count(aa)/len(v) >= 0.1:
                strong_assignments[k] += [aa]*v.count(aa)
        if 100*len(v)/total_codons_eval >= .05:
            strong_AAs = strong_assignments[k]
            if strong_AAs:
                best_cdn_guess[k] = list(max(set(strong_AAs), key=strong_AAs.count))
                best_cdn_guess[k] += [f'{100*len(v)/total_codons_eval:.3f}']
            else:
                best_cdn_guess[k] = ['Stop',f'{100*len(v)/total_codons_eval:.3f}']
        else:
            best_cdn_guess[k] = ['Stop',f'{100*len(v)/total_codons_eval:.3f}']
    return strong_assignments, best_cdn_guess



def codon_usage_table(taxon, out_folder, codon_assignments, total_codons_eval):
    codons={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S',
    'TCG':'S','TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C',
    'TGA':'Stop','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P',
    'CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I',
    'ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N',
    'AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V',
    'GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G',
    'GGG':'G'}

    strong_cdn_assign, best_cdn_guesses = keep_strong(codon_assignments, total_codons_eval)
    with open(f'{out_folder}{taxon}.CodonAssignment.New.Summary.tsv','w+') as w:
        w.write('Codon\tAssignment\tCodon_Frequency\tNotes\n')
        for k, v in codons.items():
            cdn_info = best_cdn_guesses[k]
            if cdn_info[0] == v:
                w.write(f'{k}\t'+'\t'.join(best_cdn_guesses[k])+'\n')
            else:
                w.write(f'{k}\t'+'\t'.join(best_cdn_guesses[k])+'\tReassigned\n')

if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
        taxon = sys.argv[2]
        db_loc = sys.argv[3]
    except:
        print('Usage:\n    python eval_gencode.py [FASTA-FILE-CDS] [TAXON] [AA-DATABASE]\n')
        sys.exit(1)

    out_folder = prepare_folders(taxon)
    xml_file = dmnd_blastx_compare(fasta_file, taxon, out_folder, db_loc, 1e-15, 24)
    nuc_seqs = parse_fasta(fasta_file)

    codon_assignments, total_codons = parse_xml(xml_file, nuc_seqs)
    codon_usage_table(taxon, out_folder, codon_assignments, total_codons)
