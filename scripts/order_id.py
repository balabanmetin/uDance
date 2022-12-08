import json
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='order id by time')
parser.add_argument('--backbone', type=str)
parser.add_argument('--all', type=str)
parser.add_argument('--out-dir', type=str)
args = parser.parse_args()
with open(args.backbone, 'r') as f:
    a = set(f.read().split('\n'))
    a.remove('')
with open(args.all, 'r') as f:
    b = set(f.read().split('\n'))
    b.remove('')

if len(a.intersection(b)) != len(a):
    if len(a.difference(b)) > 1 or a.difference(b).pop() != 'ROOT':
        print('Missing backbone sequences in alignment!')
        print(f'Missing: {a.difference(b)}')
        sys.exit()

query_id = list(b.difference(a))
with open(f'{args.all}', 'w') as f:
    s = list(a) + list(query_id)
    f.write("\n".join(s)+'\n')

with open(f'{args.out_dir}/query_id.txt', 'w') as f:
    f.write("\n".join(list(query_id))+'\n')

with open(f'{args.out_dir}/backbone_id.txt', 'w') as f:
        f.write("\n".join(list(a))+'\n')
