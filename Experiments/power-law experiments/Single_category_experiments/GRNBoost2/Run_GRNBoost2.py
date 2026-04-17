import pandas as pd
import numpy as np
from funcs import *
import time
import multiprocessing as mp
from arboreto.algo import grnboost2
#conda activate conda1-env
import os
import sys
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

if __name__ == "__main__":

    file_path='../generated_data/'
    file_name = 'simulation_clust_'
    file_ext='.csv'
    t = int(sys.argv[1])
    p = int(sys.argv[2])
    n_by_p= float(sys.argv[3])
    rd_seed = int(sys.argv[4]) 
    
    print('generating expression data...')
    data = Generate_Data(t,int(n_by_p*p),p, rd_seed)
    # df = pd.DataFrame(data)

    #n =  [np.shape(cts)[0] for cts in data.values()] #get array of sample size for each population
    #p = np.shape(data[1])[1]
    #print(n, p)

    # Read true Precision matrix
    print('reading true precision matrices')
    True_Precission = Read_Simulation_Data(f'../simulation_data/true_networks/clust_precision_', t, file_ext, ',')
    True_Precission=np.stack([x[:p,:p] for x in True_Precission.values()],2)

    #Create outputfile for BIC
    folder_path = os.path.join('Output', str(rd_seed))
    os.makedirs(folder_path, exist_ok=True)
    BIC_filename = os.path.join(folder_path, f'GRNBoost_t{t}.csv')
    
    if not os.path.isfile(BIC_filename):
        with open(BIC_filename, 'w') as f:
            f.write('t,p,n,precision,recall,F1,AUPRC,Time(sec)')
            f.write('\n')
            

    network={}
    TP=0
    FP=0
    FN=0
    all_scores = []
    all_truth = []
    time_network=0
    for pop in range(t):
        print(pop)
        d1=data[pop+1]
        d1.columns=d1.columns.astype(str)	
        
        tstart=time.perf_counter()
        network[pop]= grnboost2(d1)
        tend=time.perf_counter()
        
        theta=network[pop].to_numpy()
        theta=theta.astype(float)
        
        theta_mat = np.zeros((p, p))
        theta_mat[theta[:,0].astype(int),theta[:,1].astype(int)] = theta[:,2]
        np.fill_diagonal(theta_mat, 0)
        
        true_pop = True_Precission[:, :, pop].copy()
        np.fill_diagonal(true_pop, 0)
        
        TP_pop, FP_pop, FN_pop  = errorMetrics2(theta_mat,true_pop)
        
        TP+= TP_pop
        FP+= FP_pop
        FN+= FN_pop
        
        #temp=True_Precission[theta[:,0].astype(int),theta[:,1].astype(int),pop]
        #TP+=np.sum((1*(theta[:,2]>0))*(1*(temp>0))+(1*(theta[:,2]<0))*(1*(temp<0)))
        #FP+=np.sum((1*(theta[:,2]!=0))*(1*(temp==0)))
        #FN+=np.sum((1*(theta[:,2]==0))*(1*(temp!=0)))
        time_network+=tend-tstart
        
        scores = theta_mat.flatten()
        truth  = (true_pop!= 0).astype(int).flatten()

        all_scores.append(scores)
        all_truth.append(truth)

    Precision=TP/(TP+FP)
    Recall=TP/(TP+FN)
    F1=2*Precision*Recall/(Precision+Recall)
    
    # one global AUPRC over all populations together
    all_scores = np.concatenate(all_scores)
    all_truth = np.concatenate(all_truth)
    
    precision_curve, recall_curve, _ = precision_recall_curve(all_truth, all_scores)
    AUPRC = auc(recall_curve, precision_curve)



    with open(BIC_filename, 'a') as f:
        f.write(str(rd_seed)+','+str(t)+','+str(p)+','+str(n_by_p*p)+','+str(Precision)+','+str(Recall)+','+str(F1)+','+str(AUPRC)+','+str(time_network))
        f.write('\n')



    
        
