#!/usr/bin/env python3

"""
Write a python script that takes in input a list of identifiers and a fasta file and extract the sequences.
"""

import sys

def get_seqs(fileseq, pos=1):
    """Read a fasta file and create a dictionary where the keys are the ids and the values are the correspondent sequences"""
    dseq={}
    with open(fileseq, "r") as f:
        for line in f:  
            if line[0]==">":
                k = line[1:].rstrip().split("|")[pos] #line[1:] to exclude directly the > sign
                dseq[k]='' 
                continue #go to the next line
            dseq[k] = dseq[k] + line.rstrip()
    return dseq

if __name__=="__main__":
    fileseq = sys.argv[1]
    fileids = sys.argv[2]
    pos=1
    if len(sys.argv)>3: 
        pos=int(sys.argv[3])-1 #if there is a 3rd input use that to search the i
    dseq = get_seqs(fileseq, pos)
    #print(len(dseq.keys()))
    #Create a list of the identifiers from the file 
    ids=open(fileids, "r").read().rstrip().split('\n')
    for i in ids:
        if dseq.get(i,0) != 0: #if the id is not present the program will not crash
            print(">"+i+"\n"+dseq[i])
        else:
            print('WARNING: Sequence '+i+' not found', file=sys.stderr)