import treeswift as ts

treestr = input().strip()
t = ts.read_tree(treestr, schema="newick")
print(t.diameter())
