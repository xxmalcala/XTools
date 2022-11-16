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
	if 'gene:' in seq_name:
		acc = seq_name.split()[0]
		gname = seq_name.split('gene:')[1].split()[0]
	elif 'locus_tag=' in seq_name:
		acc = '_'.join(seq_name.split()[0].split('_')[-3:-1])
		gname = seq_name.split('locus_tag=')[1].split(']')[0]
	elif 'gene=' in seq_name:
		acc = '_'.join(seq_name.split()[0].split('_')[-3:-1])
		gname = seq_name.split('gene=')[1].split(']')[0]
	# for GenBank sourced CDSs!
	elif 'jgi' in seq_name:
		acc = seq_name.split('|')[-1].replace('-','_')
		gname = acc.split('.t')[0].replace('-','_')
	else:
		# print(seq_name)
		return None, None
	return gname, acc

# Check to make sure that the CDSs being compared are complete, with both
# valid start and stop codons, and without potential frameshifts.
def check_valid_seq(seq):
	ends = ['TGA','TAG','TAA']
	starts = ['ATG','ATT','ATC','ATA','TTG','CTG','GTG']
	if seq[:3] in starts and seq[-3:] in ends and len(seq) % 3 == 0:
		return 'Valid'

# Parses the FASTA files and returns just the largest isoform for each gene.
# This is not always the best option, but simplifies phylogenomic analyses.
def capture_isoforms(fasta_file):
	reduced_fasta = f'{".".join(fasta_file.split("/")[-1].split(".")[:-1])}' \
					f'.LargeIso.fasta'
	seq_len_db = defaultdict(int)
	seq_acc_db = defaultdict(str)
	seqs_seen = 0
	for i in SeqIO.parse(fasta_file,'fasta'):
		seqs_seen += 1
		gene_name, acc_name = get_gene_name(i.description)
		if gene_name and len(i.seq) > seq_len_db[gene_name]:
			if check_valid_seq(f'{i.seq.upper()}'):
				print('y')
				seq_len_db[gene_name] = len(i.seq)
				seq_acc_db[gene_name] = f'>{acc_name}\n{i.seq}\n'
	final_seq_count = len(seq_acc_db)
	if final_seq_count > 0:

		with open(f'primary_transcripts/{reduced_fasta}','w+') as w:
			w.write(''.join(seq_acc_db.values()))

		print(f'For {fasta_file.split("/")[-1]}: There were {seqs_seen} ' \
			f'originally, with {final_seq_count} unique loci prepared.')
		return 'success', reduced_fasta

	else:
		print(f'{fasta_file.split("/")[-1]} is likely not from Ensembl nor NCBI!')
		return 'fail', reduced_fasta

def force_parse(fasta_file, final_fasta):
	with open(f'primary_transcripts/{final_fasta}','w+') as w:
		for i in SeqIO.parse(f'{fasta_file}','fasta'):
			check_valid = check_valid_seq(f'{i.seq.upper()}')
			if check_valid == 'Valid':
				seq_name = i.id.split()[0].split('|')[0].replace("-","_")
				w.write(f'>{seq_name}\n{i.seq.upper()}\n')

if __name__ == '__main__':
	if len(sys.argv[1:]) == 1:
		fasta_to_prep = sys.argv[1]
		force_name = False
	elif len(sys.argv[1:]) == 2:
		fasta_to_prep = sys.argv[1]
		force_name = True
	else:
		print(f'\nMissing inputs. Usage is:\n'
		f'\n\tpython primary_isoform.py [folder/fasta-file]\n')
		sys.exit()

	prep_out_folder()

	bad_files = []
	if os.path.isdir(fasta_to_prep):

		for f in os.listdir(fasta_to_prep):
			if f.split('.')[-1].lower() in ['fa','fas','fasta','fna']:

				eval = capture_isoforms(f'{fasta_to_prep}/{f}')
				if eval[0] == 'fail':
					if force_name:
						print('Forced processing. Note: NO ISOFORM selection will be '\
						'performed.')
						force_parse(f'{fasta_to_prep}/{f}', eval[1])
					bad_files.append(f)
		if bad_files:
			with open('Primary_Isoform.Error_Files.txt','w+') as w:
				w.write('\n'.join(bad_files))
	else:
		eval = capture_isoforms(fasta_to_prep)
		if eval[0] == 'fail':
			if force_name:
				print('Forced processing. Note: NO ISOFORM selection will be '\
				'performed.')
				force_parse(fasta_to_prep, eval[1])
			bad_files.append(fasta_to_prep)
