#! /usr/bin/env python
import sys

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def readDupDict(f):
    res=dict()
    for line in f.readlines():
        x = list(line.strip().split(","))
        res[x[0]]=x[1:]
    return res

if __name__ == "__main__":
    f=open(sys.argv[1])
    g=open(sys.argv[2])
    dct = readDupDict(g)
    g.close()
    
    reqs = dict()
    for name, seq, qual in readfq(f):
        reqs[name]=seq
        if name in dct:
            for i in dct[name]:
                reqs[i]=seq
    
    f.close()
    res=[]
    for k,v in reqs.items():
        res.append(">"+k)
        res.append(v)
    with open(sys.argv[3],"w") as of:
        of.write("\n".join(res))
        of.write("\n")
