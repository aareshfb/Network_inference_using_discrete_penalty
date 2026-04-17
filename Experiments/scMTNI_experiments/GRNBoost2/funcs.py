# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:21:42 2024

@author: aares
"""
import pandas as pd 
import numpy as np
import os
import time

import copy
from scipy.special import erfinv



def nearPSD(A,epsilon=1e-6):#Found this script online need to verify it
    n = A.shape[0]
    eigval, eigvec = np.linalg.eig(A)
    val = np.matrix(np.maximum(eigval,epsilon))
    vec = np.matrix(eigvec)
    T = 1/(np.multiply(vec,vec) * val.T)
    T = np.matrix(np.sqrt(np.diag(np.array(T).reshape((n)) )))
    B = T * vec * np.diag(np.array(np.sqrt(val)).reshape((n)))
    out = B*B.T
    return(out)

def Calculate_Bic(Theta,data,n):
    p,_,t=np.shape(Theta)
    bic=0
    
    for i in range(t):
        theta_psd=nearPSD(Theta[:,:,i]) #find nearest PSD matrix to estimate
        temp = copy.deepcopy(theta_psd)
        temp[abs(temp) < 1e-6] = 0
        df=np.count_nonzero(np.triu(temp,1))
        
        COV=np.cov(data[i+1][:,:p].T)
        
        bic+=n[i]*np.trace(COV@theta_psd) - n[i]*np.log(np.linalg.det(theta_psd)) + np.log(n[i])*df + 4*df*np.log(p)
    return(bic)

def errorMetrics(Theta,True_Precission):
    TP=(1*(Theta>0))*(1*(True_Precission>0))+(1*(Theta<0))*(1*(True_Precission<0))
    FP=(1*(Theta!=0))*(1*(True_Precission==0))
    FN=(1*(Theta==0))*(1*(True_Precission!=0))
    
    TP=np.sum(TP)
    FP=np.sum(FP)
    FN=np.sum(FN)
    
    Precision=TP/(TP+FP)
    Recall=TP/(TP+FN)
    F1=2*Precision*Recall/(Precision+Recall)
    
    print('Precision:',Precision,'\nRecall:',Recall,"\nF1:",F1,"\n")
    return(Precision,Recall,F1)
    
