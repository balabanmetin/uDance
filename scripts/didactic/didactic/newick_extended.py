import treeswift


# store bracket open/close for convenience in label parsing
BRACKET = {
    '[': ']', # square bracket
    '{': '}', # curly bracket
    "'": "'", # single-quote
    '"': '"', # double-quote
}

INVALID_NEWICK = "Tree not valid Newick tree"


def read_tree_newick(newick):
    '''Read a tree from a Newick string or file

    Args:
        ``newick`` (``str``): Either a Newick string or the path to a Newick file (plain-text or gzipped)

    Returns:
        ``Tree``: The tree represented by ``newick``. If the Newick file has multiple trees (one per line), a ``list`` of ``Tree`` objects will be returned
    '''
    if not isinstance(newick, str):
        try:
            newick = str(newick)
        except:
            raise TypeError("newick must be a str")

    ts = newick.strip()

    try:
        t = treeswift.Tree(); t.is_rooted = ts.startswith('[&R]')
        if ts[0] == '[':
            ts = ']'.join(ts.split(']')[1:]).strip(); ts = ts.replace(', ',',')
        n = t.root; i = 0
        while i < len(ts):
            # end of Newick string
            if ts[i] == ';':
                #if i != len(ts)-1 or n != t.root:
                if i != len(ts)-1:
                    raise RuntimeError(INVALID_NEWICK)

            # go to new child
            elif ts[i] == '(':
                c = treeswift.Node(); n.add_child(c); n = c

            # go to parent
            elif ts[i] == ')':
                n = n.parent

            # go to new sibling
            elif ts[i] == ',':
                n = n.parent; c = treeswift.Node(); n.add_child(c); n = c

            # edge length
            elif ts[i] == ':':
                i += 1; ls = ''
                while ts[i] != ',' and ts[i] != ')' and ts[i] != ';':
                    if ts[i] == "{":
                        ei=''
                        i += 1
                        while ts[i] != "}":
                            ei += ts[i]
                            i += 1
                        n.edge_index = int(ei)
                        break
                    ls += ts[i]; i += 1
                if ls[0] == '[':
                    n.edge_params = ']'.join(ls.split(']')[:-1]); ls = ls.split(']')[-1]
                n.edge_length = float(ls)

            # node label
            else:
                label = ''; bracket = None
                while bracket is not None or ts[i] in BRACKET or (ts[i] != ':' and ts[i] != ',' and ts[i] != ';' and ts[i] != ')'):
                    if ts[i] in BRACKET and bracket is None:
                        bracket = ts[i]
                    elif bracket is not None and ts[i] == BRACKET[bracket]:
                        bracket = None
                    label += ts[i]; i += 1
                i -= 1; n.label = label
            i += 1
    except Exception as e:
        print(e)
        raise RuntimeError("Failed to parse string as Newick: %s"%ts)
    return t
