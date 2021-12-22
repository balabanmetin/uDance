
def expand_dedupe_newick(treestr, dups):
    for i in dups:
        if len(i) >= 2:
            if treestr.find("("+i[0]+":") != -1:
                ex = "((" + ":0,".join(i) + ":0)1:"
                treestr = treestr.replace("("+i[0]+":", ex, 1)
            elif treestr.find(","+i[0]+":") != -1:
                ex = ",(" + ":0,".join(i) + ":0)1:"
                treestr = treestr.replace(","+i[0]+":", ex, 1)
    return treestr
