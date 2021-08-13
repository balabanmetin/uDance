from collections import deque
ZERO_LEN = 10**(-5)

# initialize properties of input tree and return set containing taxa of leaves
def prep(tree, support):
    tree.resolve_polytomies(); tree.suppress_unifurcations()
    for node in tree.traverse_postorder():
        node.color = -1
        if not hasattr(node, 'placements'):
            node.placements = []
            node.resolved_randomly = True
        else:
            node.resolved_randomly = False
        if node.is_leaf():
            pass
        else:
            try:
                node.confidence = float(str(node))
            except:
                node.confidence = 100. # give edges without support values support 100
            if node.confidence < support: # don't allow low-support edges
                node.edge_length = float('inf')


def paint(node, color):
    s = deque()
    s.append(node)
    while len(s) != 0:
        n = s.pop()
        if n.color < 0:
            n.color = color
            s.extend(n.children)


def min_tree_coloring_sum(tree, thr):
    prep(tree, 0)
    color = 0
    for current in tree.traverse_postorder():
        if current == tree.root:
            current.weight = float("inf")
        else:
            current.weight = len(current.placements)
        if current.is_leaf():
            current.weight += 1  # add leaf to it
        else:
            left, right = current.children
            if left.weight + right.weight + current.weight <= thr:
                current.weight += left.weight + right.weight
            elif left.weight + right.weight <= thr:
                paint(left, color); paint(right, color); color += 1
            else:
                heavier, lighter = (left, right) if left.weight > right.weight else (right, left)
                paint(heavier, color); color += 1
                if lighter.weight + current.weight <= thr:
                    current.weight += lighter.weight
                else:
                    paint(lighter, color); color += 1


def min_tree_coloring_sum_max(tree, thr, max_thr):
    prep(tree, 0)
    color = 0
    for current in tree.traverse_postorder():
        if current == tree.root:
            current.weight = float("inf")
        else:
            current.weight = len(current.placements)
        if current.is_leaf():
            current.weight += 1  # add leaf to it
            current.farthest = 0
        else:
            left, right = current.children
            if left.weight + right.weight + current.weight <= thr or \
                    left.edge_length + left.farthest + right.edge_length + right.farthest < max_thr or \
                    (left.edge_length <= ZERO_LEN and len(left.placements) > 0) or \
                    (right.edge_length <= ZERO_LEN and len(right.placements) > 0) or \
                    (not current.is_root() and current.edge_length <= ZERO_LEN and len(current.placements) > 0):
                current.weight += left.weight + right.weight
                current.farthest = max(left.edge_length + left.farthest, right.edge_length + right.farthest)
            elif left.weight + right.weight <= thr:
                paint(left, color); paint(right, color); color += 1
                current.farthest = 0
            else:
                heavier, lighter = (left, right) if left.weight > right.weight else (right, left)
                paint(heavier, color); color += 1
                current.farthest = lighter.farthest + lighter.edge_length
                if lighter.weight + current.weight <= thr:
                    current.weight += lighter.weight
                else:
                    paint(lighter, color); color += 1
                    current.farthest = 0
