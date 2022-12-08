import dendropy as dy


def compute_bipartition_alignment(tree_str):
    dytree = dy.Tree.get(data=tree_str, schema="newick")
    dytree.encode_bipartitions()

    bmatrix = []
    for b in dytree.bipartition_encoding[:-1]:
        bmatrix.append(str(b))

    bmatrix_inverted = list(map(lambda y: "".join(list(y)), zip(*bmatrix)))[::-1]

    labels = [">" + t for t in dytree.taxon_namespace.labels()]

    assert len(labels) == len(bmatrix_inverted)

    return "\n".join(map(lambda x: "\n".join(list(x)), zip(labels, bmatrix_inverted))) + "\n"
