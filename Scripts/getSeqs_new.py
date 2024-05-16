#!/usr/bin/env python3

"""
Write a python script that takes in input a list of identifiers and a fasta file and extract the sequences.
"""
from getSeqs import get_seqs
import sys
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
       '--seq',
       action='store',
       help = 'File of sequences'
    )
    
    parser.add_argument(
       '--ids',
       action='store',
       help = 'List of ids that you want the sequences'
    )
    
    parser.add_argument(
       '--pos',
       action='store',
       help = 'In which position is the id respect to the pipe in the FASTA file'
    )
    
    args = parser.parse_args()
    
    
    fileseq = args.seq
    fileids = args.ids 
    try:
        pos = int(args.pos) - 1
    except:
        pos=1
        
    try:    
        dseq = get_seqs(fileseq, pos)
        #Create a list of the identifiers from the file 
        ids=open(fileids, "r").read().rstrip().split('\n')
        for i in ids:
            if dseq.get(i,0) != 0: #if the id is not present the program will not crash
                print(">"+i+"\n"+dseq[i])
            else:
                print('WARNING: Sequence '+i+' not found', file=sys.stderr)
    except:
        print('ERROR: wrong arguments', file=sys.stderr)


"""
user@LAPTOP:My_Proj$ py getSeqs_new.py --help
usage: getSeqs_new.py [-h] [--seq SEQ] [--ids IDS] [--pos POS]

options:
  -h, --help  show this help message and exit
  --seq SEQ   File of sequences
  --ids IDS   List of ids that you want the sequences
  --pos POS   In which position is the id respect to the pipe
"""