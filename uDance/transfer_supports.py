import treeswift as ts
import sys

#$1 source
#$2 destination

t1 = ts.read_tree_newick(sys.argv[1])
t2 = ts.read_tree_newick(sys.argv[2])

for n1,n2 in zip(t1.traverse_postorder(),t2.traverse_postorder()):
    if n1.is_root() or n2.is_root():
        if not n1.is_root() and n2.is_root():
            raise Exception("Postorder traversals of the nodes of the tree are not identical!")
        continue
    if n1.is_leaf() and n2.is_leaf():
        if n1.label == n2.label:
            continue
        else:
            raise Exception(f"Postorder traversals of the nodes of the tree are not identical! The order changed at {n1.label} {n2.label}")
    if (n1.is_leaf() and not n2.is_leaf()) or (not n1.is_leaf() and n2.is_leaf()):
        raise Exception(f"Postorder traversals of the nodes of the tree are not identical! The order changed at {n1.label} {n2.label}")
    else:
        n2.label = n1.label
print(t2)

