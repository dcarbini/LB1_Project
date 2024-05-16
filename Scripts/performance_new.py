#!/usr/bin/env python3

"""Given in input the preds file write a script that calculate the confusion matrix and all the the parameters to evaluate the performance given a threshold provided in input."""

import performance as pf
from getSeqs import get_seqs
import sys
import argparse
import numpy as np

def get_data(predfile, label, metr): 
    """Reading the info has to identify to each seq the e-value and the label"""
    preds = [] #output = list of lists
    with open(predfile, "r") as f:
        for line in f:
            v = line.rstrip().split() #v = [seqID, metric, label]
            v[1] = float(v[metr]) #e-value
            v[2] = int(v[label]) #label
            preds.append(v)
    return preds


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
      'P',
       type = str,
       help = 'Prediction File'
    )
    
    parser.add_argument(
       '--th',
       action='store',
       help = 'Threshold'
    )
    
    parser.add_argument(
       '--label',
       action='store',
       help = 'Column of the labels'
    )
    
    parser.add_argument(
       '--metric',
       action='store',
       help = 'Column of the metric to evaluate'
    )
    
    parser.add_argument(
       '--cm',
       action='store_true',
       help = 'Print the confusion matrix'
    )
    
    args = parser.parse_args()
    
    predfile = args.P
    try:
        metric = int(args.metric) - 1
        label = int(args.label) - 1
        preds = get_data(predfile, label, metric)
    except:
        preds = pf.get_data(predfile)
   
    try:
        th = float(args.th) #float because it is an e-value
        cm = pf.compute_cm(preds,th)
    except: 
        cm = pf.compute_cm(preds)  
        
    if args.cm:       
        print(cm)
        """
        OUTPUT:
        [[TN, FN]
        [FP, TP]]
        """
        print('TP=', cm[1][1], 'TN=', cm[0][0], 'FN=', cm[0][1], 'FP=', cm[1][0])
    
    q2 = pf.get_accuracy(cm)
    mcc = pf.get_mcc(cm)
    tpr = pf.get_tpr(cm)
    ppv = pf.get_ppv(cm)
    f1 = pf.get_F1_score(cm)
    
    print('TH=', th, 'Q2=', q2, 'MCC=', mcc, 'TPR=', tpr, 'PPV=', ppv, 'F1=', f1)


"""
user@LAPTOP:My_Proj$ py performance_new.py -h
usage: performance_new.py [-h] [--th TH] [--label LABEL] [--metric METRIC] [--cm] P

positional arguments:
  P                Prediction File

options:
  -h, --help       show this help message and exit
  --th TH          Threshold
  --label LABEL    Column of the labels
  --metric METRIC  Column of the metric to evaluate
  --cm             Print the confusion matrix
"""