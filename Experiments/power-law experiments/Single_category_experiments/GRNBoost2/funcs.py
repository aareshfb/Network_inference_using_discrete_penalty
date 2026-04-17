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

def Generate_Data(T,n,p,rd_seed):
    np.random.seed(rd_seed)
    fname='/nfs/turbo/umms-ukarvind/shared_data/l0-networks/simulations-hypergraphs/simulation_data/hypergraph_'+str(T)+'/covariance_matrix/clust_covariance_'
    #fname='../simulation_data/covariance-matrix/clust_covariance_'
    Data={}
    for i in range(1,T+1):
        cov=pd.read_csv(fname+str(i)+'.csv').to_numpy()
        cov=np.linalg.cholesky(cov)
        A=np.random.rand(n*p,1)
        A=np.sqrt(2)*erfinv(2*A-1)
        A=np.reshape(A, (n,p),'F')
        Data[i]=pd.DataFrame( A@cov)
        print("Population",i,"completed", end="\r")
    return(Data)
    
def Read_Data(file_name, file_ext='.txt', sep='\t'):
    #read expresion counts data and return list
    Data={}
    Name=['LE_counts','CT_counts','MVP_counts', 'PAN_counts']
    for i in range(len(Name)):
        Data[i]=pd.read_csv(file_name+Name[i]+file_ext,sep=sep).to_numpy()
        
    return(Data)
    
#functions for data input and output
def Read_Simulation_Data(file_name,t, file_ext='.txt', sep='\t'):
    #read expresion counts data and return list
    Data={}
    for i in range(1, t+1):
        Data[i-1]=pd.read_csv(file_name+str(i)+file_ext,sep=sep).to_numpy()
    return(Data)



    # t2=time.perf_counter()
    # time_mins=(t2-t1)/60
    # print("Total Time Taken",time_mins,"mins")
    return(theta,time_mins)

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
    TP=(1*(Theta!=0))*(1*(True_Precission!=0))
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
    
def errorMetrics2(Theta,True_Precission):
    TP=(1*(Theta!=0))*(1*(True_Precission!=0))
    FP=(1*(Theta!=0))*(1*(True_Precission==0))
    FN=(1*(Theta==0))*(1*(True_Precission!=0))
    
    TP=np.sum(TP)
    FP=np.sum(FP)
    FN=np.sum(FN)
    
    return(TP,FP,FN)
