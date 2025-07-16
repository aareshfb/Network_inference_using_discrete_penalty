# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 17:35:27 2023

@author: aares
"""

import numpy as np
from Break_points_fun2 import *
# import time
# from Verify import Solve_using_GP,Verification
# import matplotlib.pyplot as plt

def argmin_g(Node,xr): 
    C_xr=Node.Cx+np.array([[0],[Node.q*xr],[0]])
    with np.errstate(divide='ignore',invalid='ignore'): 
        opt_array=np.nan_to_num(-0.25*np.square(C_xr[1])/C_xr[0])+C_xr[2]
        opt_idx=np.argmin(opt_array)
        if opt_array[opt_idx]>Node.f0:
            return(0)

        xstar=-0.5*C_xr[1]/C_xr[0]
        return(xstar[opt_idx])

def Para_Algo(Q,c,lam,M=5e3):
    '''
    

    Parameters
    ----------
    Q : numpy nxn dim array
        Matrix of the quadratic term.
    c : numpy n dim array
        Vector for the linear term.
    lam : numpy n dim array
        Sparsity regularizer.
    M : number, optional
        Bounds on solution x. The default is 5e3.

    Returns
    -------
    Objective value, Array.

    '''
    # if (lam==np.zeros(len(c))).all():
        # return(-0.5*c.T@np.linalg.inv(Q)@c,-np.linalg.inv(Q)@c)
    # print('Setting Variable Bounds as:',M)
    class N2:
        def __init__(self,name,Cx,Dx,f0):
            #object attributes
            self.name=name
            self.Cx=Cx
            self.Dx=Dx
            self.f0=f0
        
        def show_f(self):
            print('Cx:',self.Cx)
            print('Dx:',self.Dx)
            
        def get_values(self):
            return(self.Dx,self.Cx,self.f0)
            
        def opt(self):
            Cx=self.Cx
            f0=self.f0
            opt_array=-0.25*np.square(Cx[1])/Cx[0]+Cx[2]
            opt_idx=np.argmin(opt_array)
            opt=np.min(opt_array)
            xstar=-0.5*Cx[1]/Cx[0]
            if f0<opt:
                # print(opt)
                return(f0,0)
            else:
                return(opt,xstar[opt_idx])
            
        def plot(self):
            print('Cx plotted for node',self.name)
            x=np.linspace(-75, 75,100)
            plt.figure(self.name)
            plt.plot(np.outer(Cx[0], np.square(x)).T+np.outer(Cx[1], x).T+np.reshape(Cx[2], (1,-1)))
            plt.grid()
    
    root=0
    Leaves=np.where(np.count_nonzero(Q,1)-1==1)[-1].tolist()
    if root in Leaves: Leaves.remove(root)


    ccl_Dict,rank=Graph_Hierarchy(Q,root)
    rankv=list(rank.values())
    rankk=list(rank.keys())
    
    Nodes={}
    for leaf in Leaves:
        Cx=np.array([[Q[leaf,leaf]/2],[c[leaf]],[lam[leaf]]])
        Dx=[-np.inf,np.inf]
        f0=0
        Nodes[leaf]=N2(leaf,Cx,Dx,f0)

    Levels=max(rankv)
    for level in range(Levels+1):
        # level=1
        # print(level)
        level_idx=np.where(np.array(rankv)==level)[0]
        
        for l_elm in range(len(level_idx)):
            # l_elm=0
            # ndx=rankk[level_idx[l_elm]]
            ndx=rankk[level_idx[l_elm]]
            Cndx=[]
            Dndx=[]
            for branch in ccl_Dict[ndx]:
                # print('Branch:',branch)
                for i in range(len(branch)-1):
                    q=Q[branch[i],branch[i+1]] 
                    Dx,Cx,f0=Nodes[branch[i]].get_values()
                    Dx,Cx=Truncate(Dx, Cx,M)
                    Dy,Cy=fi_to_gj(Dx, Cx, f0, q)
                    Nodes[branch[i]].Cy=copy.deepcopy(Cy)
                    Nodes[branch[i]].q=q
                    if branch[i+1]!=ndx:
                        # f0_next=np.min(Cy[2,:]) #lowest coefficent note fj(0)=gj(0)
                        f0_next=Nodes[branch[i]].opt()[0]
                        h_x=np.array([[Q[branch[i+1],branch[i+1]]/2],[c[branch[i+1]]],[lam[branch[i+1]]]]) #g_sum_to_fj(Cy, Q[branch[i+1],branch[i+1]], c[branch[i+1]], lam[branch[i+1]])
                        Cy_next=Cy+h_x
                        Nodes[branch[i+1]]=N2(branch[i+1],Cy_next,Dy,f0_next)
                        # print('Node',i+1,'calculated')
                    else:
                        Cndx.append(Cy)
                        Dndx.append(Dy)
            Dy,Cy=Sum_g(Dndx,Cndx)
            f0_next=np.min(Cy[2,:])
            h_x=np.array([[ Q[ndx,ndx]/2],[c[ndx]],[lam[ndx]]])
            Cy_next= Cy+h_x
            Nodes[ndx]=N2(ndx,Cy_next,Dy,f0_next)
            
            '''
            #### Temp code
            N_idx=[]
            for branch in ccl_Dict[ndx]:
                N_idx=N_idx+branch[:-1]
            N_idx=N_idx +[ndx]
            # print(Solve_using_GP(N_idx, Q, c, lam))
            # print('Parametric sol:',Nodes[ndx].opt())
            '''
    x={ndx:Nodes[ndx].opt()[1]}
    for level in reversed(range(Levels+1)):
        # print(level)
        level_idx=np.where(np.array(rankv)==level)[0]
        for l_elm in range(len(level_idx)):
            ndx=rankk[level_idx[l_elm]]
            # print(ndx)
            for branch in ccl_Dict[ndx]:
                # print('Branch:',branch)
                branch_rev=list(reversed(branch))
                for i in range(len(branch)-1):
                    x[branch_rev[i+1]]=argmin_g(Nodes[branch_rev[i+1]],x[branch_rev[i]])
                    
    xnp=np.array([x[i] for i in range(len(x.keys()))])
    return(Nodes[root].opt()[0],xnp)