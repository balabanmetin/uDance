import treeswift as ts
import collections
import json
import json
import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description='add back duplicates')
parser.add_argument('--rm-duplicate', type=str)
parser.add_argument('--tree-file', type=str)
parser.add_argument('--out-dir', type=str)
args = parser.parse_args()

tree = ts.read_tree_newick(args.tree_file)
rm_duplicate = pd.read_csv(args.rm_duplicate, sep=', ', index_col=0, header=None, engine='python')
rm_map = set(rm_duplicate.index)
for l in tree.traverse_leaves():
    if l.label in rm_map:
        node = ts.Node(label=None, edge_length=l.edge_length)
        parent = l.parent
        l.edge_length = 0
        parent.remove_child(l)
        parent.add_child(node)
        node.add_child(l)
        for item in rm_duplicate.loc[l.label]:
            if not item:
                break 
            child = ts.Node(label=item, edge_length=0)
            node.add_child(child)
tree.write_tree_newick(f'{args.out_dir}/didactic_full.nwk')

