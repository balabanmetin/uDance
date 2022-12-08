import json
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='order id by time')
parser.add_argument('--out-dir', type=str)
args = parser.parse_args()

subdir = []
for f in os.listdir(args.out_dir):
    if os.path.isdir(os.path.join(args.out_dir, f)) and f.isnumeric():
        subdir.append(os.path.join(args.out_dir, f))

check = False
while not check:
    cnt = 0
    for d in subdir:
        if os.path.isfile(os.path.join(d, 'seqs', 'RUN.treefile')):
            cnt += 1
    if cnt == len(subdir):
        check = True
