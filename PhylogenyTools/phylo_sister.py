#!/usr/bin/env python3

import glob, os, sys
import numpy as np
from collections import defaultdict
from pathlib import Path
from ete3 import Tree

"""This script requires a FOLDER of newick formatted phylogenetic trees (note that
ETE3 can handle bootstrapped or Sh-aLRT values at nodes, not both!), then "walks"
through each tree reporting a summary of the sequence/taxon as well as its
phylogenetic sisters.

Option to report just data from single-sister relationships (e.g. sequence sister
to a single other sequence versus clade of sequences) to be done later.

Outputs include a large PER sequence summary table as well as a summary by taxon.

IMPORTANT! Taxon names are currently expected to follow PhyloToL naming conventions.
"MajorClade_MinorClade_Taxon", for example (Opisthokonta, Metazoa, Homo sapiens --> Op_me_Hsap).

Major clade codes are currently limited to:
Opisthokonta: Op
Archaeplastida: Pl
Amoebozoa: Am
Excavata: Ex
SAR: Sr
Other eukaryotes (e.g. cryptophytes, "haptista"): EE
Bacteria: Ba
Archaea: Za
"""

# Args to add (to follow Katz lab options):
# single-sister (default being ALL)
# Anything else?! Perhaps multi-cell? Not sure the limit on number (e.g. 1-10)
# Taxon name length changes or "code" table ...
def check_args():
    pass

# Parse tree, check taxon naming convention and return tree-object if "parseable".
def parse_tree(tree_file):
    bad_taxon_codes = None
    t = Tree(tree_file)

    all_taxa = list(set([node.name[:10] for node in t.get_leaves()]))

    mnr_clades = set([taxon[:5] for taxon in all_taxa])
    mjr_clades = set([taxon[:2] for taxon in all_taxa if taxon[:2] in
        ['Ba','Za', 'Op','Pl','Am','Ex','Sr']])

    if [c for c in all_taxa if [len(j) for j in c.split('_')] != [2,2,4]]:
        bad_taxon_codes = True

    elif not (mjr_clades or mnr_clades):
        bad_taxon_codes = True

    if bad_taxon_codes:
        print('ERROR: Taxon-names need to conform to the "PhyloToL" style.')
        print('\nExample: Homo sapiens (Opisthokonta, metazoa, homo sapiens) --> Op_me_Hsap')
        sys.exit(1)

    return t

# Re-root the phylogenetic tree prior to traversal. Returns rooted phylogenetic tree.
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

# Checks taxon names of all descendants from node. Moves "up" if same taxon and returns
# node name, whether the leaves are the same taxon, and number of iterations.
def check_same_taxon(node, reps = 0):
    taxon_node_names = set([leaf.name[:10] for leaf in node.up.get_leaves()])
    if len(taxon_node_names) == 1:
        parent = node.up
        return parent, True, reps + 1

    else:
        return node, False, reps

# Performs the phylogenetic sister evaluations. Summarizes the data on a per
# sequence basis, also making inferences on type of branch length (i.e. short
# vs long) based on an overly simple metric (mean/median phyloenetic branch lengths).
def check_sisters(tree, gene_fam, blen_dist, blen_mode):
    tree_phylo_sister_summary = defaultdict(list)

    if blen_mode == 'average':
        thresh_blen = np.mean(blen_dist)
    elif blen_mode == 'median':
        thresh_blen = np.median(blen_dist)

    for node in tree.get_leaves():
        taxon_eval = [node, True, 0]
        query_taxon = node.name[:10]
        mjr = 'same-major'
        mnr = 'same-minor'
        qblen_type = 'long'

        while taxon_eval[1]:
            taxon_eval = check_same_taxon(taxon_eval[0], taxon_eval[2])

        sister_seqs = [taxon.name for taxon in taxon_eval[0].up.get_leaves()]
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

        if node.up.dist < thresh_blen:
            qblen_type = 'short'

        if len(sister_taxa) > 20:
            sister_taxa = sister_seqs = ['Too-Many']

        tree_phylo_sister_summary[query_taxon].append(
            [gene_fam,
            node.name,
            mjr,
            mnr,
            ','.join(sister_taxa),
            ','.join(sister_seqs),
            f'{node.up.dist:.5f}',
            f'{thresh_blen:.5f}',
            qblen_type])

    return tree_phylo_sister_summary

# Solely performs the individual tree preparation, root adjustment, and scoring.
def score_tree(tree_file, gene_fam, reroot_folder, blen_mode):
    tree = parse_tree(tree_file)
    if not tree:
        return None

    tree, blen_dist = adjust_tree_root(tree)
    save_tree(tree, gene_fam, reroot_folder)

    return check_sisters(tree, gene_fam, blen_dist, blen_mode)

# Collects the phylogenetic tree summary for all trees in a given folder.
def check_many_trees(tree_folder, blen_mode):
    reroot_folder = f'{tree_folder.rstrip("/")}_ReRooted_Trees/'
    Path(reroot_folder).mkdir(parents = True, exist_ok = True)

    comp_summary = defaultdict(list)

    all_tree_files = glob.glob(f'{tree_folder}/*.tre*')

    for tree_file in all_tree_files:
        gene_tree = tree_file.split('/')[-1].replace('_postguidance','')
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

# Saves the re-rooted phylogenies.
def save_tree(tree, gene_fam, reroot_folder):
    out_tree = f'{reroot_folder}{gene_fam}.ReRooted.tre'
    tree.write(format=1, outfile= out_tree)

# Generates summary tables using data from all the phylogenies in a given folder.
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

    summarize_results(trees_summary, out_folder.rstrip('/'), blen_filt = False, verbose = True)
