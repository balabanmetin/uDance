import numpy as np


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


def fasta2dic(ref_fp, prot_flag, mask_flag):
    refs = {}
    with open(ref_fp) as f:
        mask_translation = str.maketrans('abcdefghijklmnopqrstuvwxyz', '-' * 26)

        if prot_flag:
            invalid_translation = str.maketrans('BJOUXZ', '-' * 6)
        else:
            invalid_translation = str.maketrans('BDEFHIJKLMNOPQRSUVWXYZ', '-' * 22)

        def makeupper(s):
            if mask_flag:
                return s.translate(mask_translation)
            else:
                return s.upper()

        for name, seq, qual in readfq(f):
                refs[name] = np.frombuffer(seq.encode(), dtype='S1')
    return refs
