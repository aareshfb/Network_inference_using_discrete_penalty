import pandas as pd
import numpy as np
from funcs import *
import time
import multiprocessing as mp
import sys


if __name__ == "__main__":

    file_path='../generated_data/'
    file_name = 'simulation_clust_'
    file_ext='.csv'
    
    t = int(sys.argv[1])
    p = int(sys.argv[2])
    n_by_p = float(sys.argv[3])
    nu0 = float(sys.argv[4])
    rand_seed = int(sys.argv[5]) 
    
    #t=20
    #p=2000
    #n_by_p=0.5
    #read expression data for all samples
    # print('reading expression data...')
    # data = Read_Simulation_Data(file_path+file_name, t,file_ext,',')
    print('generating expression data...')
    data = Generate_Data(t,int(n_by_p*p),p,rand_seed)

    n =  [np.shape(cts)[0] for cts in data.values()] #get array of sample size for each population
    p = np.shape(data[1])[1]
    print(n, p)

    # Read true Precision matrix
    print('reading true precision matrices')
    True_Precission = Read_Simulation_Data(f'../simulation_data/true_networks/clust_precision_', t, file_ext, ',')
    True_Precission=np.stack([x[:p,:p] for x in True_Precission.values()],2)

    #Create outputfile for BIC
    output_path = os.path.join("Output", "Vary_np", str(rand_seed))
    os.makedirs(output_path, exist_ok=True)
    
    n_by_p_str = str(n_by_p).replace('.', 'd')
    BIC_filename = os.path.join(output_path, f'ErrorMetrics_L0_np{n_by_p_str}.csv')
    print(f"BIC filename: {BIC_filename}")

    #BIC_filename='Output/Vary_np/ErrorMetrics_L0.csv' #/home/aareshfb/simulation-compare_l1_vs_l0/code_L0/Output/Vary_np

    if not os.path.isfile(BIC_filename):
        with open(BIC_filename, 'w') as f:
            f.write('rand_seed,t,p,n,nu_0,mu,gamma,bic value,precision,recall,F1,Time(min)')
            f.write('\n')
            
            
    #parameters for inference
    #nu0 = 1e-4
    S = backwardMap(data, nu0)

    for i in range(t):
        print(f'\n Number of non-diagonal edges in Pop-{i} : {np.count_nonzero(S[:,:,i]) - p}')

    gamma=1 #similarity
    mu=1e-2 #sparsity

    
    W=np.genfromtxt('../simulation_data/MST-adjacency_matrix.txt')
    W=np.triu(W)
    # W = np.array([[0,0,1,0,0,0,0],
    #          [0,0,0,0,1,0,0],
    #          [0,0,0,0,1,1,0],
    #          [0,0,0,0,0,0,1],
    #          [0,0,0,0,0,0,0],
    #          [0,0,0,0,0,0,1],
    #          [0,0,0,0,0,0,0]])
    #    W = np.diag(np.ones(t-1),1)

        
    muList = np.array([1e-4,1e-3])
#    muList = np.arange(1e-3, 1e-2, 1e-3)


    for mu in muList:
        print(nu0, mu, gamma)
        
        print('Estimating Theta ....')
        time0=time.time()
        Theta,time_network = estimateNetwork(S, W, mu, gamma);
        time1=time.time()
        print('Estimation complete ....')
        
        print('Number of errors:',np.count_nonzero(np.isnan(Theta)))
        #print number of non-zero edges in recovered network
        print('Number of non-zero edges in recovered network:')
        for j in range(t):
            y = Theta[:,:,j]
            print(len(y[abs(y) > 0]) - p)

        print('Evaluating bic of estimated model....')
        bic=Calculate_Bic(Theta, data, n)

        print('Error metrics for network estimation \n')
        Precision,Recall,F1=errorMetrics(Theta, True_Precission)

        with open(BIC_filename, 'a') as f:
            f.write(str(rand_seed)+','+str(t)+','+str(p)+','+str(n_by_p*p)+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(bic)+','+str(Precision)+','+str(Recall)+','+str(F1)+','+str(time_network))
            f.write('\n')



    
        
