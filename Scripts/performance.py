#!/usr/bin/env python3

"""Given in input the preds file write a script that calculate the confusion matrix and all the the parameters to evaluate the performance given a threshold provided in input."""

import sys
import numpy as np

def get_data(predfile): 
    """Reading the info has to identify to each seq the e-value and the label"""
    preds = [] #output = list of lists
    with open(predfile, "r") as f:
        for line in f:
            v = line.rstrip().split() #v = [seqID, e-value, label]
            v[1] = float(v[1]) #e-value
            v[2] = int(v[2]) #label
            preds.append(v)
    return preds
    
#CONFUSION MATRIX (CM)
def compute_cm(preds,th=0.5): 
    cm = np.zeros((2,2)) #matrix 2x2
    #Go line by line and assign the correct values
    for pred in preds: 
        p=0  #prediction is 0 by default
        if pred[1]<=th: p=1 #if e-value lower than the threshold the prediction is positive
        cm[p][pred[2]]+=1 #Identify the cell [prediction][label=reality] and add one case
    return cm #return the confusion matrix TN, FN, FP, TP
        
def _cells(cm):
    tp = cm[1][1]
    tn = cm[0][0]
    fn = cm[0][1]
    fp = cm[1][0]
    return tp,tn,fn,fp
        
#OVERALL ACCURACY        
def get_accuracy(cm): #(TN+TP)/(TN+TP+FN+FP) = (TN+TP)/tot
    q2 = float((cm[0][0]+cm[1][1])/np.sum(cm))
    #print('Q2=', q2, 'N=', np.sum(cm))
    return q2

#MATTHEW CORRELATION COEFFICIENT -> if the dataset is not balanced 
def get_mcc(cm): #(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    tp,tn,fn,fp = _cells(cm)
    mcc = (tp*tn-fp*fn)/np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return mcc

#SENSITIVITY OR RECALL OR TRUE POSITIVE RATE (TPR)
def get_tpr(cm): #TP/(TP+FN)
    tp,tn,fn,fp = _cells(cm)
    tpr = tp/(tp+fn)
    return tpr
    
#PRECISION OR POSITIVE PREDICTIVE VALUE (PPV)
def get_ppv(cm): #TP/(TP+FP)
    tp,tn,fn,fp = _cells(cm)
    ppv = tp/(tp+fp)
    return ppv
    
#F1
def get_F1_score(cm): #2*PPV*TPR/(PPV+TPR)
    ppv = get_ppv(cm)
    tpr = get_tpr(cm)
    f1 = 2*ppv*tpr/(ppv+tpr)
    return f1

if __name__=="__main__":
    predfile = sys.argv[1]
    preds = get_data(predfile)
    #print(preds)
    if len(sys.argv)<2: 
        cm = compute_cm(preds)
    else:
        th = float(sys.argv[2]) #float because it is an e-value
        cm = compute_cm(preds,th)
    print(cm)
    print('TP=', cm[1][1], 'TN=', cm[0][0], 'FN=', cm[0][1], 'FP=', cm[1][0])
    """
    OUTPUT:
    [[TN, FN]
    [FP, TP]]
    """
    q2 = get_accuracy(cm)
    mcc = get_mcc(cm)
    tpr = get_tpr(cm)
    ppv = get_ppv(cm)
    f1 = get_F1_score(cm)
    #print('Total sequences:', np.sum(cm))
    print('TH=', th, 'Q2=', q2, 'MCC=', mcc, 'TPR=', tpr, 'PPV=', ppv, 'F1=', f1)