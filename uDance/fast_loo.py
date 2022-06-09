#!/usr/bin/env python3

from collections import deque
import treeswift as ts
import sys
import functools
import heapq


class PrioritySet(object):
    def __init__(self):
        self.heap = []
        self.set = set()

    def add(self, d, pri):
        if not d in self.set:
            heapq.heappush(self.heap, (pri, d))
            self.set.add(d)

    def get(self):
        pri, d = heapq.heappop(self.heap)
        self.set.remove(d)
        return d

    def __len__(self):
        return len(self.heap)

def set_levels(tree):
    root = tree.root
    root.level = 0
    q = deque()
    q.append(root)
    while len(q) != 0:
        n = q.popleft()
        for c in n.children:
            c.level = n.level + 1
        q.extend(n.children)


def distance_between(n1, n2):
    ps = PrioritySet()
    for n in [n1, n2]:
        ps.add(n.parent, -n.level)  # adding minus to give priority to the largest level

    count = 0
    while len(ps) > 1:
        x = ps.get()
        if x.true_split:
            count += 1
        ps.add(x.parent, -x.parent.level)
    mrca = ps.get()
    if mrca.true_split:
        count += 1
    return count

t = ts.read_tree_newick(sys.argv[1])
# find a backbone
for n in t.traverse_postorder(internal=False):
    if not n.label.endswith("-query"):
        rrt = n
        break

t.reroot(n)
t.suppress_unifurcations()
set_levels(t)
labs = [l for l in t.labels(internal=False) if not l.endswith("-query")]
l2n = t.label_to_node(selection="leaves")
for e in t.traverse_postorder():
    if e.is_leaf():
        if e.label.endswith("-query"):
            e.all_query = 1
        else:
            e.all_query = 0
    else:
        e.all_query = functools.reduce(lambda x, y: x and y, [c.all_query for c in e.children])

for e in t.traverse_postorder(leaves=False):
    if sum([1-f.all_query for f in e.children]) >= 2:
        e.true_split = True
    else:
        e.true_split = False

def dist_calc():
    for l in labs:
        if l in l2n and l + "-query" in l2n:
            yield l, distance_between(l2n[l], l2n[l + "-query"])-1


dists = dict(dist_calc())
for k, v in dists.items():
    print(k + "\t" + str(int(v)))
