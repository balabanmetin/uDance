import dendropy,sys

T = dendropy.Tree.get(path=sys.argv[1], schema='newick', preserve_underscores=True)
T.resolve_polytomies()
print(T.as_string(schema='newick'))


