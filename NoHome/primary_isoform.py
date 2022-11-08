#!/usr/bin/env python3

'''This script captures and returns the largest isoform for each gene from a
fasta file of CDSs, or folder of fasta files, sourced from Ensembl.'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com
# Last Modified: 2020-04-02
# usage: python primary_isoform.py [folder/fasta]

# Dependencies:
# Python3, BioPython


import os
import sys
from collections import defaultdict
from Bio import SeqIO


def prep_out_folder():
	out_folder = 'primary_transcripts'
	if not os.path.isdir(out_folder):
		os.mkdir(out_folder)


# Captures the unique gene name from the Ensembl or GenBank CDSs.
def get_gene_name(seq_name):
	if 'protein_id=' in seq_name:
		acc = seq_name.split('protein_id=')[-1].split(']')[0]
		gname = acc
	elif 'gene:' in seq_name:
		acc = seq_name.split()[0]
		gname = seq_name.split('gene:')[1].split()[0]
	elif 'gene=' in seq_name:
		acc = seq_name.split()[0]
		gname = seq_name.split('gene=')[1].split()[0]
	# for GenBank sourced CDSs!
	elif 'locus_tag=' in seq_name:
		acc = seq_name.split()[0].split('_')[-2]
		gname = seq_name.split('locus_tag=')[1].split(']')[0]
	else:
		print(seq_name)
		sys.exit()
	return gname, acc

# Check to make sure that the CDSs being compared are complete, with both
# valid start and stop codons, and without potential frameshifts.
def check_valid_seq(seq):
	ends = ['TGA','TAG','TAA']
	starts_bac = ['ATG','ATT','ATC','ATA','TTG','CTG','GTG']
	starts = ['ATG']
	if seq[:3] in starts and seq[-3:] in ends and len(seq) % 3 == 0:
		return 'Valid'


# Parses the FASTA files and returns just the largest isoform for each gene.
# This is not always the best option, but simplifies phylogenomic analyses.
def capture_isoforms(fasta_file):
	reduced_fasta = f'{".".join(fasta_file.split("/")[-1].split(".")[:-1])}'\
					f'.LargeIso.fasta'
	seq_len_db = defaultdict(int)
	seq_acc_db = defaultdict(str)

	for i in SeqIO.parse(fasta_file,'fasta'):
		if 'pseudo=true' in i.description.lower():
			continue
		gene_name, acc_name = get_gene_name(i.description)
		if len(i.seq) > seq_len_db[gene_name]:
			if check_valid_seq(f'{i.seq.upper()}'):
				seq_len_db[gene_name] = len(i.seq)
				seq_acc_db[gene_name] = f'>{acc_name}\n{i.seq}\n'

	with open(f'primary_transcripts/{reduced_fasta}','w+') as w:
		w.write(''.join(seq_acc_db.values()))

if __name__ == '__main__':
	try:
		fasta_to_prep = sys.argv[1]

	except IndexError:
		print(f'\nMissing inputs. Usage is:\n'
		f'\n\tpython primary_isoform.py [folder/fasta-file]\n')
		sys.exit()

	prep_out_folder()

	if os.path.isdir(fasta_to_prep):
		for f in os.listdir(fasta_to_prep):
			capture_isoforms(f'{fasta_to_prep}/{f}')
	else:
		capture_isoforms(fasta_to_prep)
