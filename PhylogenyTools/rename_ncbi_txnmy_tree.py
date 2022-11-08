
from ete3 import Tree
import sys

def get_names(node_file):
    names = {i.split()[0]:i.split('\t')[2] for i in open(node_file).readlines()}
    return names

def rename_tree(tree_file, node_file):
    names = get_names(node_file)
    t = Tree(tree_file)
    for leaf in t.get_leaves():
        try:
            ncbi_name = names[leaf.name]
            leaf.name = ncbi_name
        except:
            pass
    # return t
    t.write(format=1, outfile=tree_file.replace('.tre','.Renamed.tre'))

if __name__ == '__main__':
    try:
        tree_file, node_file = sys.argv[1:]
    except:
        print('usage: python rename.py tree-file.tre names.dmp')
        sys.exit()
