#!/usr/bin/env python3


import glob, os, sys
import numpy as np
from collections import defaultdict
from pathlib import Path
from ete3 import Tree


def check_args():
    pass

def parse_tree(tree_file):
    t = Tree(tree_file)
    all_taxa = list(set([node.name[:10] for node in t.get_leaves()]))
    bad_taxon_codes = [c for c in all_taxa if [len(j) for j in c.split('_')] != [2,2,4]]

    if bad_taxon_codes:
        print('ERROR: Taxon-names need to conform to the "PhyloToL" style.')
        print('\nExample: Homo sapiens (Opisthokonta, metazoa, homo sapiens) --> Op_me_Hsap')
        sys.exit(1)

    mnr_clades = list(set([taxon[:5] for taxon in all_taxa]))
    mjr_clades = list(set([taxon[:2] for taxon in all_taxa if taxon[:2] != 'EE']))

    if len(mnr_clades) < 4:
        print(f'Ignoring {tree_file} as it contains fewer than 4 "minor clades"')
        return None

    else:
        return t


def adjust_tree_root(tree):
    major_clades = ['Ba','Za', 'Op','Pl','Am','Ex','Sr']
    clade_sizes = {i:[] for i in ['BaZa', 'Op','Pl','Am','Ex','Sr']}
    blen_dist = []
    for node in tree.iter_descendants("postorder"):
        blen_dist.append(node.dist)
        mjr_c = [i.name[:2] for i in node.get_leaves() if i.name[:2] in major_clades]
        if len(mjr_c) > 1:
            if mjr_c.count('Ba') + mjr_c.count('Za') >= len(mjr_c)-1:
                clade_sizes['BaZa'].append((node, len(mjr_c)))
                # break
            else:
                for clade in major_clades[2:]:
                    if mjr_c.count(clade) >= len(mjr_c)-1:
                        clade_sizes[clade].append((node, len(mjr_c)))
    for k, v in clade_sizes.items():
        if v:
            root_node = max(v, key = lambda x: x[-1])[0]
            break
    tree.set_outgroup(root_node)

    return tree, blen_dist


def check_same_taxon(node, reps = 0):
    taxon_node_names = list(set([leaf.name[:10] for leaf in node.up.get_leaves()]))
    if len(taxon_node_names) == 1:
        parent = node.up
        return parent, True, reps + 1
    else:
        return node.up, False, reps


def check_sisters(tree, gene_fam, blen_dist, blen_mode):

    if blen_mode == 'average':
        thresh_blen = np.mean(blen_dist)
    elif blen_mode == 'median':
        thresh_blen = np.median(blen_dist)

    tree_phylo_sister_summary = defaultdict(list)

    for node in tree.get_leaves():
        taxon_eval = [node, True, 0]
        query_taxon = node.name[:10]
        mjr = 'same-major'
        mnr = 'same-minor'
        qblen_type = 'long'
        while taxon_eval[1]:
            taxon_eval = check_same_taxon(taxon_eval[0], taxon_eval[2])

        sister_seqs = [taxon.name for taxon in taxon_eval[0].get_leaves()]
        sister_taxa = list(set([i[:10] for i in sister_seqs if i[:10] != query_taxon]))

        mjr_c = list(set([i[:2] for i in sister_taxa]))
        mnr_c = list(set([i[:5] for i in sister_taxa]))

        if len(set(mjr_c)) > 1:
            mjr = 'non-monophyletic'

        elif len(set(mjr_c)) == 1 and mjr_c[0] != query_taxon[:2]:
            mjr = mjr_c[0]

        if len(set(mnr_c)) > 1:
            mnr = 'non-monophyletic'

        elif len(set(mnr_c)) == 1 and mnr_c[0] != query_taxon[:5]:
            mnr = mnr_c[0]

        if node.dist < thresh_blen:
            qblen_type = 'short'

        if len(sister_taxa) < 10:
            tree_phylo_sister_summary[query_taxon].append(
                [gene_fam,
                node.name,
                mjr,
                mnr,
                ','.join(sister_taxa),
                ','.join(sister_seqs),
                f'{node.dist:.4f}',
                f'{thresh_blen:.4f}',
                qblen_type])

    return tree_phylo_sister_summary


def score_tree(tree_file, gene_fam, reroot_folder, blen_mode):
    tree = parse_tree(tree_file)
    if not tree:
        return None

    tree, blen_dist = adjust_tree_root(tree)
    save_tree(tree, gene_fam, reroot_folder)

    return check_sisters(tree, gene_fam, blen_dist, blen_mode)


def check_many_trees(tree_folder, blen_mode):
    reroot_folder = f'{tree_folder.rstrip("/")}_ReRooted/Rooted_Trees/'
    Path(reroot_folder).mkdir(parents = True, exist_ok = True)

    comp_summary = defaultdict(list)

    all_tree_files = glob.glob(f'{tree_folder}/*.tre*')

    for tree_file in all_tree_files:
        gene_tree = tree_file.split('/')[-1]
        if gene_tree.startswith('RAxML'):
            gene_fam = gene_tree.split('.')[1]
        else:
            gene_fam = gene_tree.split('.')[0]

        detailed_summary = score_tree(
                            tree_file,
                            gene_fam,
                            reroot_folder,
                            blen_mode)

        for k, v in detailed_summary.items():
            comp_summary[k] += v

    return comp_summary, reroot_folder.split('/Rooted_Trees/')[0]


def save_tree(tree, gene_fam, reroot_folder):
    out_tree = f'{reroot_folder}{gene_fam}.ReRooted.nwk'
    tree.write(format=1, outfile= out_tree)


def summarize_results(comp_summary, out_prefix, blen_filt = True, verbose = True):
    priority = ['BaZa', 'Op','Pl','Am','Ex','Sr']

    brief_summary = {}

    for k, v in comp_summary.items():
        unique_gfs = []
        category = []
        short_blen = []

        total_seqs = len(v)

        for i in v:
            unique_gfs.append(i[0])
            if (blen_filt and i[3][:2] and i[-1] == 'short') or (not blen_filt):
                category.append(i[3][:2])
                short_blen.append(i[-1])

        brief_summary[k] = [category.count('sa'), category.count('Op'),
            category.count('Pl'), category.count('Am'), category.count('Ex'),
            category.count('Sr'), category.count('EE'), category.count('Ba'),
            category.count('Za'), category.count('no'), total_seqs,
            f'{category.count("sa")/total_seqs:.3f}']

    with open(f'{out_prefix}.PhylogenSisters.Summary.tsv','w+') as w:
        w.write('Taxon\tSame-Minor-Clade\tOp\tPl\tAm\tEx\tSr\tEE\tBa\tZa\t' \
                'Non-Monophyletic\tTotal-Tips\tProportion-Same-Minor\n')
        for k, v in brief_summary.items():
            temp = '\t'.join([f'{i}' for i in v])
            w.write(f'{k}\t{temp}\n')

    with open(f'{out_prefix}.PhylogenSisters.tsv', 'w+') as w:
        w.write('Taxon\tGene-Family\tSeq-Name\tMajor-Clade\tMinor-Clade\tSister-Taxa\t' \
            'Sister-Seqs\tBranch-Length\tBranch-Length-Threshold\tBranch-Category\n')
        for k, v in comp_summary.items():
            for i in v:
                temp = '\t'.join([f'{j}' for j in i])
                w.write(f'{k}\t{temp}\n')


if __name__ == '__main__':
    if len(sys.argv[1:]) == 1:
        tree_folder = sys.argv[1]
        blen_eval = 'average'

    elif len(sys.argv[1:]) == 2:
        tree_folder = sys.argv[1]
        blen_eval = sys.argv[2].lower()

        if blen_eval not in ['average', 'median']:
            print('Warning: Branch length evaluation is neither: AVERAGE nor MEDIAN.')
            print('Defaulting to AVERAGE')
            blen_eval = 'average'

    else:
        print('Usage:\n\n    python count_sisters.py [FOLDER-WITH-TREES] ' \
            '[AVERAGE/MEDIAN (for classifying branch-lengths)]\n')
        sys.exit(1)


    trees_summary, out_folder = check_many_trees(tree_folder, blen_eval)

    out_prefix = out_folder
    summarize_results(trees_summary, out_prefix, blen_filt = False, verbose = True)
