from itertools import groupby


def tc_parser(tc_outfp):
    with open(tc_outfp, "r") as tc_output:
        tc_output.readline()
        lines = map(lambda x: x.strip().split('\t'), tc_output.readlines())
    lines_sorted = sorted(lines, key=lambda x: x[1])
    clusters = [(key, [i[0] for i in list(group)]) for key, group in groupby(lines_sorted, lambda x: x[1])]
    return clusters
