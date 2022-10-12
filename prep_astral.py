#!/usr/bin/env python3

# Author: Xyrus
# Contact: xmaurer-alcala@amnh.org
# Updated: 2022-10-12

"""Comments to come..."""

import glob, os, sys
from pathlib import Path
from ete3 import Tree
from collections import defaultdict
import shutil

def fix_tree_naming(tree_file):
    t = Tree(tree_file)
    seq_names = {leaf.name:leaf.name[:4]+leaf.name[5:10] for leaf in t.get_leaves()}
    for leaf in t.get_leaves():
        leaf.name = seq_names[leaf.name]
    t.write(format=1, outfile=f'{tree_file}.fixed.nwk')


def make_family_file(project_name, initial_tree_dir):
    tree_files = glob.glob(f'{initial_tree_dir}/OG*.tre*')

    for f in tree_files:
        fix_tree_naming(f)

    os.system(f'cat {initial_tree_dir}/*fixed.nwk* > {project_name}.ForAstral.tre')

if __name__ == '__main__':
    if len(sys.argv[1:]) == 2:
        initial_tree_dir = sys.argv[1]
        project_name = sys.argv[2]
    else:
        print('Usage:\n    python3 prep_astral.py [FOLDER-WITH-PHYLOGENIES] [PROJECT-NAME]\n')
        sys.exit(1)
    # gs_tree_dir = prep_dirs(project_name)
    make_family_file(project_name, initial_tree_dir)
