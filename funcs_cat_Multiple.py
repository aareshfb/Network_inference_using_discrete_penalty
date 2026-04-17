# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:21:42 2024

@author: aares
"""

import numpy as np
import os
import time
import multiprocessing as mp
from Para_algo2 import Para_Algo
import matplotlib.pyplot as plt
import copy



def ST(x,nu):
    #perform soft-thresholding of sample covariance matrix (backward mapping)
    
    # y=1/n*(x.T@x)
    y=np.cov(x.T)
    d=np.diag(y)
    # y=(y > nu)*(y - nu) + (y < -nu)*(y + nu)
    y=y-np.sign(y)*np.minimum(np.abs(y),nu)
    np.fill_diagonal(y,d)
    print('Lowest eig after thresholding:',min((np.linalg.eigvals(y))))
    y=np.linalg.inv(y)
    return(y)

def backwardMap(ctsList, nu0):
    '''
    Parameters
    ----------
    Data : list of expression matrices
    nu : int
        Thresholding Coefficent.
    Returns
    -------
    Backward Mapping:(N,N,t) np.array
    '''
    
    Data = {}
    for k in ctsList.keys():
        cts = ctsList[k]
        (n,p) = np.shape(cts)
        nu = nu0 * np.sqrt(np.log(p) / n) * np.linalg.norm(cts)
        Data[k]= ST(cts,nu)

    S=np.stack([x for x in Data.values()],2)
    return(S)
        
    
def Save_Theta(Theta,n,file_location='../Output/',file_name='theta'):
    isExist=os.path.exists(file_location)
    if not isExist:
        os.makedirs(file_location)
        print('Folder',file_location,' has been created. Data in this folder.')
    T=len(Theta[0][1])
    #n=int(np.sqrt(len(Theta)))
    np_theta=np.zeros((n,n,T))
    for i in Theta:
        np_theta[i[0]]=i[1]
    # for t in range(T):
    #     np.savetxt(file_location+file_name+str(t)+'.csv', np_theta[:,:,t],delimiter=',')
    # print("Data saved in",file_location,"folder")
    #np_theta=np_theta+np.transpose(np.triu(np_theta,1),(1,0,2))
    for t in range(T):
       np_theta[:,:,t]=np_theta[:,:,t]+np.triu(np_theta[:,:,t],1).T
    return(np_theta)

def Get_Q(W,gamma):
    Q=-2*gamma*(W+W.T)
    Q=Q+2*np.diag(1+gamma*(np.sum(W,1)+np.sum(W,0)))#multiply the diag elements by 2 because 1/2xTQx
    return(Q) 


    return(Q)



def Get_Q_l2_allM(W,gamma,no_cat,Reg=0.01):
    n=len(W[0,:])
    Q=Get_Q(W, gamma)
    I=np.eye(n) 
    O=np.zeros((n,n))
    # Q=Q+np.diag(np.diag(Q)) #Scale the diag elements by 2 for the global variables
    Q=Q+2*(no_cat-1)*I #You add 2 to the diagonal element becasue you need to add 1 for the category (i) which gets scales by 2 because of 1/2 Q factor
    # Q=np.block([[Q+Reg*I,2*I,2*I],
    #             [2*I,(2+Reg)*I,O],
    #             [2*I,O,(2+Reg)*I]])

    blocks = [[Q + Reg*I] + [2*I]*no_cat]

    for i in range(no_cat):
        row = [2*I]
        for j in range(no_cat):
            if i == j:
                row.append((2 + Reg)*I)
            else:
                row.append(O)
        blocks.append(row)
    
    Q_big = np.block(blocks)
    return(Q_big)



def Parallel_MP_row(Q, C, lam, i, Theta, p):
    # print(LAM)
    for j in range(i,p):
        if i==j:
            temp=-np.linalg.inv(Q)@C[j,:]
        else:
            
            _,temp=Para_Algo(Q, C[j,:], lam*np.ones(len(Q)))
            # _,temp=Para_Algo(Q, C[j,:], LAM)
        Theta.append([(i,j),temp])
    print("Gene",i,"completed", end="\r")



def estimateNetwork_CategoryM(S_list, W, mu, gamma,p=None,output_location='../Output/'):
    '''
    Parameters
    ----------
    S : soft-thresholded sample precision matrix
    W : adjacency matrix for subpopulations
    mu : sparsity parameter
    gamma : similarity strength between connected subpopulations
    
    Returns
    -------
    Network:(N,N,t) np.array
    '''
    # Q=Get_Q_l2(W,gamma)
    no_cat=len(S_list)
    Q=Get_Q_l2_allM(W,gamma,no_cat,0.01)
    if p==None:
        p = np.shape(S_list[0])[0]
    # p = np.shape(S)[0]
    
    # LAM=mu2*np.ones(len(Q))
    # LAM[:int(len(Q)/3)]=mu
    t1=time.time()
    with mp.Manager() as manager:
        Theta=manager.list()
        #create pool object
        nprocs = mp.cpu_count()
        pool = mp.Pool(nprocs - 0)
        
        # pool.starmap_async(Parallel_MP_row, [(Q,np.hstack((-2*S_P[i,:,:]-2*S_R[i,:,:],
        #                                                    -2*S_P[i,:,:],
        #                                                    -2*S_R[i,:,:])),mu,i,Theta,p) for i in range(p)], chunksize = 5)
        
        pool.starmap_async(Parallel_MP_row, [(Q,np.hstack([-2*sum(S_cat[i, :, :] for S_cat in S_list)]
                                                          +[-2 * S_cat[i, :, :] for S_cat in S_list]),mu,i,Theta,p) for i in range(p)], chunksize=5)

        pool.close()
        pool.join()
        theta = Save_Theta(Theta,p,output_location)

    t2=time.time()
    print("Total Time Taken",(t2-t1)/60,"mins")
    return(theta)

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
        
        bic+=np.trace(theta_psd@COV) - n[i]*np.log(np.linalg.det(theta_psd)) + np.log(n[i])*df + 4*df*np.log(p)
    return(bic)

# def errorMetrics(Theta,True_Precission):
#     TP=(1*(Theta!=0))*(1*(True_Precission!=0))
#     FP=(1*(Theta!=0))*(1*(True_Precission==0))
#     FN=(1*(Theta==0))*(1*(True_Precission!=0))
    
#     TP=np.sum(TP)
#     FP=np.sum(FP)
#     FN=np.sum(FN)
    
#     Precision=TP/(TP+FP)
#     Recall=TP/(TP+FN)
#     F1=2*Precision*Recall/(Precision+Recall)
    
#     print('Precision:',Precision,'\nRecall:',Recall,"\nF1:",F1,"\n")
#     return(Precision,Recall,F1)

# def Get_Q_category(W,gamma,a):
#     #Unused
#     n=len(W[0,:])
    
#     Q=-2*gamma*(W+W.T)
#     Q=Q+2*np.diag(2*np.square(a)+gamma*(np.sum(W,1)+np.sum(W,0)))#multiply the diag elements by 2 because 1/2xTQx
    
#     I=np.eye(n) 
#     O=np.zeros((n,n))
#     # Q=Q+np.diag(np.diag(Q)) #Scale the diag elements by 2 for the global variables
#     # Q=Q #You add 2 to the diagonal element becasue you need to add 1 for the recurrect term which gets scales by 2 because of 1/2Q factor
#     Q=np.block([[Q,2*a*I,2*a*I],
#                 [2*a*I,2*I,O],
#                 [2*a*I,O,2*I]])
#     return(Q)

# def Get_Q_l2(W,gamma,Reg=0.01):
#     n=len(W[0,:])
#     Q=Get_Q(W, gamma)
#     I=np.eye(n) 
#     O=np.zeros((n,n))
#     # Q=Q+np.diag(np.diag(Q)) #Scale the diag elements by 2 for the global variables
#     Q=Q+2*I #You add 2 to the diagonal element becasue you need to add 1 for the recurrect term which gets scales by 2 because of 1/2 Q factor
#     Q=np.block([[Q+Reg*I,2*I,2*I],
#                 [2*I,(2)*I,O],
#                 [2*I,O,(2)*I]])

# def Get_Q_l2_all(W,gamma,Reg=0.01):
#     n=len(W[0,:])
#     Q=Get_Q(W, gamma)
#     I=np.eye(n) 
#     O=np.zeros((n,n))
#     # Q=Q+np.diag(np.diag(Q)) #Scale the diag elements by 2 for the global variables
#     Q=Q+2*I #You add 2 to the diagonal element becasue you need to add 1 for the recurrect term which gets scales by 2 because of 1/2 Q factor
#     Q=np.block([[Q+Reg*I,2*I,2*I],
#                 [2*I,(2+Reg)*I,O],
#                 [2*I,O,(2+Reg)*I]])
#     return(Q)

# def estimateNetwork_Category(S_P,S_R, W, mu, gamma,p=None,output_location='../Output/'):
#     '''
#     Parameters
#     ----------
#     S : soft-thresholded sample precision matrix
#     W : adjacency matrix for subpopulations
#     mu : sparsity parameter
#     gamma : similarity strength between connected subpopulations
    
#     Returns
#     -------
#     Network:(N,N,t) np.array
#     '''
#     # Q=Get_Q_l2(W,gamma)
#     Q=Get_Q_l2_all(W,gamma,0.01)
#     if p==None:
#         p = np.shape(S_P)[0]
#     # p = np.shape(S)[0]
    
#     # LAM=mu2*np.ones(len(Q))
#     # LAM[:int(len(Q)/3)]=mu
#     t1=time.time()
#     with mp.Manager() as manager:
#         Theta=manager.list()
#         #create pool object
#         nprocs = mp.cpu_count()
#         pool = mp.Pool(nprocs - 0)
        
#         pool.starmap_async(Parallel_MP_row, [(Q,np.hstack((-2*S_P[i,:,:]-2*S_R[i,:,:],
#                                                            -2*S_P[i,:,:],
#                                                            -2*S_R[i,:,:])),mu,i,Theta,p) for i in range(p)], chunksize = 5)

#         pool.close()
#         pool.join()
#         theta = Save_Theta(Theta,p,output_location)

#     t2=time.time()
#     print("Total Time Taken",(t2-t1)/60,"mins")
#     return(theta)