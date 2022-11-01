#!/usr/bin/env python3

"""Comments to come..."""

import glob, os, sys
from pathlib import Path
from ete3 import Tree
from collections import defaultdict
import shutil

def parse_tree(tree_file):
    mapping_info = defaultdict(list)
    seq_names = [leaf.name for leaf in Tree(tree_file).get_leaves()]
    for s in seq_names:
        f_name = f'{s[:4]}{s[5:10]}'
        mapping_info[f_name].append(s)
    return mapping_info

def prep_dirs(project_name):
    grax_fam_dir = f'{project_name}_SpRax/'
    Path(f'{grax_fam_dir}Trees/').mkdir(parents=True, exist_ok=True)
    Path(f'{grax_fam_dir}Mapping/').mkdir(parents=True, exist_ok=True)
    return f'{grax_fam_dir}Trees/'



def make_family_file(gs_tree_dir, initial_tree_dir):
    map_dir = gs_tree_dir.replace('Trees/','Mapping/')
    fam_file = ['[FAMILIES]']
    fam_num = 1
    tree_files = glob.glob(f'{initial_tree_dir}/OG*.tre*')
    for f in tree_files:
        shutil.copy2(f, gs_tree_dir)
        mappings = parse_tree(f)
        map_file = f'{map_dir}{f.split("/")[-1].split(".tre")[0]}.link'
        with open(map_file, 'w+') as w:
            for k, v in mappings.items():
                w.write(f'{k}:{";".join(v)}\n')
        fam_file += [f'- family {fam_num}',
                    f'starting_gene_tree = {f}',
                    f'mapping = {map_file}']
        fam_num += 1
    with open('families_speciesrax.txt', 'w+') as w:
        w.write('\n'.join(fam_file))

def sp_rax_cmd(project_name):
    sprx_cmd = 'generax --families families_speciesrax.txt ' \
        '--strategy SKIP ' \
        '--si-strategy HYBRID ' \
        '--species-tree MiniNJ ' \
        '--rec-model UndatedDL ' \
        '--per-family-rates ' \
        '--prune-species-tree ' \
        '--si-estimate-bl ' \
        '--si-quartet-support ' \
        f'--prefix {project_name}_SpRax'
    print('Example Species-Rax commandline below (change "--rec-model UndatedDL" ' \
        'to "--rec-model UndatedDTL" if you want to allow TRANSFERS!):\n')
    print(f'    {sprx_cmd}\n')


if __name__ == '__main__':
    if len(sys.argv[1:]) == 2:
        initial_tree_dir = sys.argv[1]
        project_name = sys.argv[2]
    else:
        print('Usage:\n    python3 prep_grax.py [FOLDER-WITH-PHYLOGENIES] [PROJECT-NAME]\n')
        sys.exit(1)
    gs_tree_dir = prep_dirs(project_name)
    make_family_file(gs_tree_dir, initial_tree_dir)
    sp_rax_cmd(project_name)
