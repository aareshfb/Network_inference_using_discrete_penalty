# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 12:39:34 2024

@author: aares
"""

import numpy as np
from funcs import *
# import pandas as pd
import time
import sys
# import multiprocessing as mp

# NOTE:
# The data path in this script is hard-coded.
# Please modify it (around line 35,36,67) to point to your local dataset.

if __name__ == "__main__":

    #t=int(sys.argv[1]) #[3,5,10,20,50,100]
    #p=250
    #n_by_p= float(sys.argv[2])
    #psx=0.25 #pertubation strength
    
    psx=float(sys.argv[1])
    t=int(sys.argv[2])
    p=int(sys.argv[3])
    n_by_p=float(sys.argv[4])
    rd_seed=int(sys.argv[5])
    
    
    
    n=int(n_by_p*p)
    file_path='path_to_data/vary_T/hypergraph_'+str(t)+'_P'+str(psx) #hard coded file path
    Precision_file_name='/true_networks/cat' #hard coded file path to ground truth data. see line 59,60
    # P_file_name = '/covariance_matrix/cat1_clust_covariance_'
    # R_file_path='../simulation_data/cat2/'
    # R_file_name = 'simulation_clust_'
    file_ext='.csv'
    #output_location='../Output/'#+file_name+'/'
    

    #Generate Data
    data_P=Generate_Data(t,1,n,p,psx,rd_seed)
    data_R=Generate_Data(t,2,n,p,psx,rd_seed)
    
    
    #read data
    #print('reading expression data...')
    #data_P = Read_Data(P_file_path+P_file_name,'.csv',',')
    #data_R = Read_Data(R_file_path+R_file_name,'.csv',',')
    
    print([np.shape(cts)[0] for cts in data_P.values()])
    print([np.shape(cts)[0] for cts in data_R.values()])
    #p = np.shape(data_P[1])[1]
    
    print('reading true precision matrices')
    True_Precission_Cat1 = Read_Simulation_Data(file_path+Precision_file_name+str(1)+'_clust_precision_', t, file_ext, ',')
    True_Precission_Cat2 = Read_Simulation_Data(file_path+Precision_file_name+str(2)+'_clust_precision_', t, file_ext, ',')
    True_Precission_Cat1=np.stack([x[:p,:p] for x in True_Precission_Cat1.values()],2)
    True_Precission_Cat2=np.stack([x[:p,:p] for x in True_Precission_Cat2.values()],2)
    
    #print(f'size of problems is: {n, p}')
    
    # W = np.diag(np.ones(T-1),1) #similarity matrix for sub-populations
    W=np.genfromtxt(file_path+'/MST.txt')
    W=np.triu(W)
    
    #parameters for inference
    muList = np.array([1e-4,1e-3])
    nu0List = np.array([5e-3,1e-3,5e-4,1e-4])#np.array([1e-3,5e-4,1e-4,5e-5,1e-5])
    gamma=1 #similarity
    
    psx_str = str(psx).replace('.', 'd')
    output_dir = os.path.join('Output', psx_str, str(rd_seed))
    os.makedirs(output_dir, exist_ok=True)
    BIC_filename = os.path.join(output_dir, f'perturb_T{t}.csv')
    if not os.path.isfile(BIC_filename):
        with open(BIC_filename, 'w') as f:
            f.write('rd_seed,Alg,t,p,n,Cat,nu_0,mu,gamma,precision,recall,F1,Time,BIC,BIC_mean')
            f.write('\n')
    
    for nu0 in nu0List:
        for mu in muList:
            S_P = backwardMap(data_P, nu0) #backward map of sample covariance matrix
            S_R = backwardMap(data_R, nu0) #backward map of sample covariance matrix              
            for i in range(t):
                print(f'\n Number of non-diagonal edges in Primary Pop-{i} : {np.count_nonzero(S_P[:,:,i]) - p}')
                print(f'\n Number of non-diagonal edges in Recurrent Pop-{i} : {np.count_nonzero(S_R[:,:,i]) - p}')
        
            # plt.plot(np.diag(S[:pp,:pp,0]))
            # print(np.diag(S[:pp,:pp, 0]))
            # Regulaizer=8.15 #this is to ensure Q is psd use >8.15
        
    
    
    
            print(nu0, mu,gamma)
    
    
            print('Estimating Theta (Cat_alg)....')
            time1_c=time.time()
            Theta = estimateNetwork_Category(S_P,S_R, W, mu, gamma,p)
            time2_c=time.time()
            time_c=time2_c-time1_c
            print('Estimation complete ....')
            
            print('Estimating Theta cat1 ....')
            time1_p=time.time()
            Theta_P = estimateNetwork(S_P, W, mu, gamma);
            time2_p=time.time()
            time_p=time2_p-time1_p
            
            print('Estimating Theta cat2 ....')
            time1_r=time.time()
            Theta_R = estimateNetwork(S_R, W, mu, gamma);
            time2_r=time.time()
            time_r=time2_r-time1_r
            print('Estimation complete ....')
            
            for i in range(t):
                print(f'\n Number of non-diagonal edges in Pop-{i} : {np.count_nonzero(Theta[:,:,i]) - p}')
    
            ## take sum with global network to get primary/recurrent networks
            
            
            #Save Data
            ## global networks, primary and recurrent
            #for i in range(t*3):
            #     np.savetxt(f'{output_location}/{i}.csv', Theta[:,:,i],delimiter=',')
    
            #print("Data saved in",output_location,"folder")
            Cat1_theta=np.zeros((p,p,t))
            Cat2_theta=np.zeros((p,p,t))
            for i in range(t):
                 Cat1_theta[:,:,i]=Theta[:,:,i]+Theta[:,:,i+t]
                 Cat2_theta[:,:,i]=Theta[:,:,i]+Theta[:,:,i+2*t]
                 #np.savetxt(f'{output_location}/Cat1_{i}.csv', Theta[:,:,i]+Theta[:,:,i+T],delimiter=',')
                 #np.savetxt(f'{output_location}/Cat2_{i}.csv', Theta[:,:,i]+Theta[:,:,i+2*T],delimiter=',')
         
            print('Error metrics for network estimation \n')
            print(nu0, mu,gamma)
            Precision_Cat1,Recall_Cat1,F1_Cat1=errorMetrics(Cat1_theta, True_Precission_Cat1)
            Precision_Cat2,Recall_Cat2,F1_Cat2=errorMetrics(Cat2_theta, True_Precission_Cat2)
    
            Precision_Cat1_direct,Recall_Cat1_direct,F1_Cat1_direct=errorMetrics(Theta_P, True_Precission_Cat1)
            Precision_Cat2_direct,Recall_Cat2_direct,F1_Cat2_direct=errorMetrics(Theta_R, True_Precission_Cat2)
            
            N=[n]*t
            bic_cat_1=Calculate_Bic(Cat1_theta, data_P, N)
            bic_cat_2=Calculate_Bic(Cat2_theta, data_R, N)

            bic_direct_1=Calculate_Bic(Theta_P, data_P, N)
            bic_direct_2=Calculate_Bic(Theta_R, data_R, N)


            
            with open(BIC_filename, 'a') as f:
                 f.write(str(rd_seed)+',Cat_alg,'+str(t)+','+str(p)+','+str(n_by_p*p)+','+'Cat1'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat1)+','+str(Recall_Cat1)+','+str(F1_Cat1)+','+str(time_c)+','+str(bic_cat_1)+','+str(0.5*(bic_cat_1+bic_cat_2)))
                 f.write('\n')
                 f.write(str(rd_seed)+',Cat_alg,'+str(t)+','+str(p)+','+str(n_by_p*p)+','+'Cat2'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat2)+','+str(Recall_Cat2)+','+str(F1_Cat2)+', ,'+str(bic_cat_2))
                 f.write('\n')
                 f.write(str(rd_seed)+',Direct_alg,'+str(t)+','+str(p)+','+str(n_by_p*p)+','+'Cat1'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat1_direct)+','+str(Recall_Cat1_direct)+','+str(F1_Cat1_direct)+','+str(time_p)+','+str(bic_direct_1)+','+str(0.5*(bic_direct_1+bic_direct_2)))
                 f.write('\n')
                 f.write(str(rd_seed)+',Direct_alg,'+str(t)+','+str(p)+','+str(n_by_p*p)+','+'Cat2'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat2_direct)+','+str(Recall_Cat2_direct)+','+str(F1_Cat2_direct)+','+str(time_r)+','+str(bic_direct_2))
                 f.write('\n')




    
