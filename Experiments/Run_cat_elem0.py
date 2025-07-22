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
import os
# import multiprocessin1g as mp

if __name__ == "__main__":

    # t=5 #[3,5,10,20]
    # p=250
    # n_by_p=20
    # n=n_by_p*p
    
    perturb_strength=float(sys.argv[1])
    t=int(sys.argv[2])
    p=int(sys.argv[3])
    n_by_p=float(sys.argv[4])
    rd_seed=int(sys.argv[5])
    
    n=int(n_by_p*p)
    file_path='/nfs/turbo/umms-ukarvind/shared_data/l0-networks/simulations-category/simulation_data/perturb_strength/perturb_'+str(perturb_strength)
    Precision_file_name='/true_networks/cat' #clust_precision_
    # P_file_name = '/covariance_matrix/cat1_clust_covariance_'
    # R_file_path='../simulation_data/cat2/'
    # R_file_name = 'simulation_clust_'
    file_ext='.csv'
    #output_location='../Output/'#+file_name+'/'
    

    #Generate Data
    data_P=Generate_Data(t,perturb_strength,1,n,p,rd_seed)
    data_R=Generate_Data(t,perturb_strength,2,n,p,rd_seed)
    
    
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
    W=np.genfromtxt('/nfs/turbo/umms-ukarvind/shared_data/l0-networks/simulations-category/simulation_data/perturb_strength/MST.txt')
    W=np.triu(W)
    
    #parameters for inference
    muList = np.array([1e-4,1e-3])
    nu0List =np.array([5e-3,1e-3,5e-4,1e-4])
    gamma=1 #similarity
    
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

            print('Evaluating bic of estimated model....')
            bic_cat_P=Calculate_Bic(Cat1_theta, data_P, [n]*t)
            bic_cat_R=Calculate_Bic(Cat2_theta, data_R, [n]*t)
            bic_direct_P=Calculate_Bic(Theta_P, data_P, [n]*t)
            bic_direct_R=Calculate_Bic(Theta_R, data_R, [n]*t)
            
            perturb_str = str(perturb_strength).replace('.', 'd')

            # Create the directory if it doesn't exist
            output_dir = os.path.join('Output', str(rd_seed))
            os.makedirs(output_dir, exist_ok=True)
            
            # Final filename
            BIC_filename = os.path.join(output_dir, f'perturb_{perturb_str}.csv')
            
            
            #BIC_filename=file_path+'/Error_Metrics_perturb'+'.csv'
            if not os.path.isfile(BIC_filename):
                with open(BIC_filename, 'w') as f:
                    f.write('rd_seed,Alg,perturb_strength,t,p,n_by_p,Cat,nu_0,mu,gamma,precision,recall,F1,Time,bic,bic_mean')
                    f.write('\n')
            
            with open(BIC_filename, 'a') as f:
                 f.write(str(rd_seed)+',Cat_alg,'+str(perturb_strength)+','+str(t)+','+str(p)+','+str(n_by_p)+','+'Cat1'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat1)+','+str(Recall_Cat1)+','+str(F1_Cat1)+','+str(time_c)+','+str(bic_cat_P)+','+str(0.5*(bic_cat_P+bic_cat_R)))
                 f.write('\n')
                 f.write(str(rd_seed)+',Cat_alg,'+str(perturb_strength)+','+str(t)+','+str(p)+','+str(n_by_p)+','+'Cat2'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat2)+','+str(Recall_Cat2)+','+str(F1_Cat2)+', ,'+str(bic_cat_R))
                 f.write('\n')
                 f.write(str(rd_seed)+',Direct_alg,'+str(perturb_strength)+','+str(t)+','+str(p)+','+str(n_by_p)+','+'Cat1'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat1_direct)+','+str(Recall_Cat1_direct)+','+str(F1_Cat1_direct)+','+str(time_p)+','+str(bic_direct_P)+','+str(0.5*(bic_direct_P+bic_direct_R)))
                 f.write('\n')
                 f.write(str(rd_seed)+',Direct_alg,'+str(perturb_strength)+','+str(t)+','+str(p)+','+str(n_by_p)+','+'Cat2'+','+str(nu0)+','+str(mu)+','+str(gamma)+','+str(Precision_Cat2_direct)+','+str(Recall_Cat2_direct)+','+str(F1_Cat2_direct)+','+str(time_r)+','+str(bic_direct_R))
                 f.write('\n')



#     heatmap(np.abs(Theta),output_location)

#    A = Theta[:,:,0]
#    print(np.count_nonzero(A,1))
#    plt.plot(np.diag(A))
    
#    MU = [3e-7,1e-6,3e-6,1e-5]
#    GAMMA = [1e-1]

#    BIC_filename='../Output/BIC-networks.txt'
#        
#    if not os.path.isfile(BIC_filename):
#        with open(BIC_filename, 'w') as f:
#            f.write('nu_0,mu,gamma,bic value')
#            
#    bic_list = []

#    for mu in MU:
#        for gamma in GAMMA:
#            
#            print(mu,gamma)
#                
#            print('Estimating Theta ....')
#            Theta = estimateNetwork(S, W, mu, gamma,p,output_location);
#            A = Theta[:,:,0]
#            print(np.count_nonzero(A,1))
#            plt.plot(np.diag(A))
#            
#            print('Estimation complete ....')
#    
#            print('Evaluating bic of estimated model....')
#            bic=Calculate_Bic(Theta, data, n)
#            print(bic)
#            bic_list.append(bic)
#    
#            with open(BIC_filename, 'a') as f:
#                f.write(str(nu0)+'\t'+str(mu)+'\t'+str(gamma)+'\t'+str(bic))
#                f.write(f'\n Number of non-diagonal edges in Pop-1 : {np.count_nonzero(Theta[:,:,0]) - p}')
#                f.write(f'\n Number of non-diagonal edges in Pop-2 : {np.count_nonzero(Theta[:,:,1]) - p}')
#                f.write(f'\n Number of non-diagonal edges in Pop-3 : {np.count_nonzero(Theta[:,:,2]) - p}')
#                f.write(f'\n Number of non-diagonal edges in Pop-4 : {np.count_nonzero(Theta[:,:,3]) - p}')
#                f.write('\n')
       

        
    # W=np.array([[ 0., 1., 1.],
    #             [0., 0., 0.],
    #             [0., 0., 0.]])
    
