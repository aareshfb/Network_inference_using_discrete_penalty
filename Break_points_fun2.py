# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 18:50:10 2023

@author: aareshfb

Algorithm to find the breakpoints
"""

import numpy as np
import copy
import networkx as nx

def Project_Breakpoints_Robust(C,Dx,epsilon=1e-6):
    """
    Parameters
    ----------
    C : np.array size(3,n)
        Coefficents of f(x)
    Dx : List
        Breakpoints in f(x)
    epsilon: float
        Error Margin
    Returns
    -------
    Breakpoints in g(y), index of coefficents of g(y)
    """
    n=np.shape(C)[1]
    Dy=[-np.inf]
    Discard=[]
    for j in range(1,n):
        # print('j',j,'Dy:',Dy,'Discard',Discard)
        i=j-1
        while i in Discard:
            i-=1

        while i>=0:
            # print('i',i)
            try:
                y1,y2=Get_Slopes(C[:,i], C[:,j])
            except:
                y1=Get_Slopes(C[:,i], C[:,j])[0]
                y2=np.inf
                
            yli,yui=Slope_bounds(C[:,i], Dx[i:i+2])
            ylj,yuj=Slope_bounds(C[:,j], Dx[j:j+2])
            
            # x1i=Get_intercept(C[:,i], y1)
            # x1j=Get_intercept(C[:,j], y1)
            
            # x2i=Get_intercept(C[:,i], y2)
            # x2j=Get_intercept(C[:,j], y2)
            
            # if (x1j>=Dx[j] and x1j<=Dx[j+1]) and (x1i>=Dx[i] and x1i<=Dx[i+1]) and not isinstance(y1, complex):
            if ((yli-y1)/np.abs(np.maximum(yli,-1e300))<=epsilon and (y1-yui)/np.abs(yui)<= epsilon) and ((ylj-y1)/np.abs(ylj)<=epsilon and (y1-yuj)/np.abs(np.minimum(yuj,1e300))<= epsilon) and not isinstance(y1, complex):
                if y1>=max(Dy):
                    Dy.append(y1)
                    i=-1
                else:
                    del Dy[-1]
                    Discard.append(i)
                    while i in Discard:
                        i-=1
                        
            # elif (x2j>=Dx[j] and x2j<=Dx[j+1]) and (x2i>=Dx[i] and x2i<=Dx[i+1])and not isinstance(y2, complex):
            elif ((yli-y2)/np.abs(np.maximum(yli,-1e300))<=epsilon and (y2-yui)/np.abs(yui)<= epsilon) and ((ylj-y2)/np.abs(ylj)<=epsilon and (y2-yuj)/np.abs(np.minimum(yuj,1e300))<=epsilon) and not isinstance(y2, complex):
                if y2>=max(Dy):
                    Dy.append(y2)
                    i=-1
                else:
                    del Dy[-1]
                    Discard.append(i)
                    while i in Discard:
                        i-=1 
                '''edits
            else:
                if Dy[-1] != -np.inf: del Dy[-1]
                Discard.append(i)
                while i in Discard:
                    i-=1 
                '''
            else:
                

                # if (x1i<Dx[itemp] or x1i>Dx[itemp+1]) or (x2i<Dx[itemp] or x1i>Dx[itemp+1]):# or isinstance(y1, complex):
                if (y1<yli or y1>yui) and (y2<yli or y2>yui):    
                    Discard.append(i)
                    if Dy[-1] != -np.inf: del Dy[-1]
                    while i in Discard:
                        i-=1 
                else:
                    i=-1
                # if (x1j<Dx[j] or x1j>Dx[j+1]) or (x2j<Dx[j] or x1j>Dx[j+1]):# or isinstance(y1, complex):
                if (y1<ylj or y1>yuj) and (y2<ylj or y2>yuj):
                    Discard.append(j)
                    # itemp=i
                    # i=-1

            '''
            else:  
                Discard.append(j)
                i=-1
            '''
    Dy.append(np.inf)
    Cy_idx=np.setdiff1d(np.arange(0,n),Discard)
    return(Dy,Cy_idx)




def Project_Breakpoints(C,Dx):
    """
    Parameters
    ----------
    C : np.array size(3,n)
        Coefficents of f(x)
    Dx : List
        Breakpoints in f(x)

    Returns
    -------
    Breakpoints in g(y), index of coefficents of g(y)
    """
    n=np.shape(C)[1]
    Dy=[-np.inf]
    Discard=[]
    for j in range(1,n):
        # print('j',j,'Dy:',Dy,'Discard',Discard)
        i=j-1
        while i in Discard:
            i-=1

        while i>=0:
            # print('i',i)
            y1,y2=Get_Slopes(C[:,i], C[:,j])
            
            yli,yui=Slope_bounds(C[:,i], Dx[i:i+2])
            ylj,yuj=Slope_bounds(C[:,j], Dx[j:j+2])
            
            # x1i=Get_intercept(C[:,i], y1)
            # x1j=Get_intercept(C[:,j], y1)
            
            # x2i=Get_intercept(C[:,i], y2)
            # x2j=Get_intercept(C[:,j], y2)
            
            # if (x1j>=Dx[j] and x1j<=Dx[j+1]) and (x1i>=Dx[i] and x1i<=Dx[i+1]) and not isinstance(y1, complex):
            if (yli<=y1 and y1<= yui) and (ylj<=y1 and y1<= yuj) and not isinstance(y1, complex):
                if y1>=max(Dy):
                    Dy.append(y1)
                    i=-1
                else:
                    del Dy[-1]
                    Discard.append(i)
                    while i in Discard:
                        i-=1
                        
            # elif (x2j>=Dx[j] and x2j<=Dx[j+1]) and (x2i>=Dx[i] and x2i<=Dx[i+1])and not isinstance(y2, complex):
            elif (yli<=y2 and y2<= yui) and (ylj<=y2 and y2<= yuj) and not isinstance(y2, complex):
                if y2>=max(Dy):
                    Dy.append(y2)
                    i=-1
                else:
                    del Dy[-1]
                    Discard.append(i)
                    while i in Discard:
                        i-=1 
                '''edits
            else:
                if Dy[-1] != -np.inf: del Dy[-1]
                Discard.append(i)
                while i in Discard:
                    i-=1 
                '''
            else:
                

                # if (x1i<Dx[itemp] or x1i>Dx[itemp+1]) or (x2i<Dx[itemp] or x1i>Dx[itemp+1]):# or isinstance(y1, complex):
                if (y1<yli or y1>yui) and (y2<yli or y2>yui):    
                    Discard.append(i)
                    if Dy[-1] != -np.inf: del Dy[-1]
                    while i in Discard:
                        i-=1 
                else:
                    i=-1
                # if (x1j<Dx[j] or x1j>Dx[j+1]) or (x2j<Dx[j] or x1j>Dx[j+1]):# or isinstance(y1, complex):
                if (y1<ylj or y1>yuj) and (y2<ylj or y2>yuj):
                    Discard.append(j)
                    # itemp=i
                    # i=-1

            '''
            else:  
                Discard.append(j)
                i=-1
            '''
    Dy.append(np.inf)
    Cy_idx=np.setdiff1d(np.arange(0,n),Discard)
    return(Dy,Cy_idx)

def Add_indicator2(Dy,Cy,f0): #This has a bug
    yl=-np.inf
    yu=np.inf
    if np.shape(Cy)[1]==0:
        return(Dy,Cy)
    for i in range(np.shape(Cy)[1]):
        temp1,temp2=np.roots(Cy[:,i]-np.array([0,0,f0]))
        y1=np.minimum(temp1,temp2)
        y2=np.maximum(temp1,temp2)
        if isinstance(y1, complex):
            return(Dy,Cy)
        else:
            yl=np.maximum(yl,y1)
            yu=np.minimum(yu,y2)
    temp1=list(np.where(Dy<yl)[0])
    temp2=list(np.where(Dy>yu)[0])
    D2y=[Dy[i] for i in temp1]+[yl,yu]+[Dy[i] for i in temp2]
    C2y=np.concatenate((Cy[:,:np.max(temp1)+1], np.array([[0,0,f0]]).T, Cy[:,np.min(temp2)-1:]),1)#note, Dy has +-inf and Cy does not
    return(D2y,C2y)

def Add_indicator3(Dy,Cy,f0):
    yl=-np.inf
    yu=np.inf
    if np.shape(Cy)[1]==0:
        return(Dy,Cy)
    for i in range(np.shape(Cy)[1]):
        temp1,temp2=np.roots(Cy[:,i]-np.array([0,0,f0]))
        y1=np.minimum(temp1,temp2)
        y2=np.maximum(temp1,temp2)
        if isinstance(y1, complex):
            return(Dy,Cy)
        else:
            yl=np.maximum(yl,y1)
            yu=np.minimum(yu,y2)
    if yl>=yu:
        return(Dy,Cy)
    temp1=list(np.where(Dy<yl)[0])
    temp2=list(np.where(Dy>yu)[0])
    D2y=[Dy[i] for i in temp1]+[yl,yu]+[Dy[i] for i in temp2]
    C2y=np.concatenate((Cy[:,:np.max(temp1)+1], np.array([[0,0,f0]]).T, Cy[:,np.min(temp2)-1:]),1)#note, Dy has +-inf and Cy does not
    return(D2y,C2y)

    
def Truncate(Dx,Cx,lim=5e3):
    Dt=np.array(copy.deepcopy(Dx))
    Ct=copy.deepcopy(Cx)
    idx=np.where(Dt>=-lim)[0]
    # print(idx)
    Dt=Dt[idx]
    Ct=Ct[:,np.min(idx)-1:]
    Dt=np.append([-np.inf], Dt)
    idx=np.where(Dt<=lim)[0]
    # print(idx)
    Dt=Dt[idx]
    Ct=Ct[:,:np.max(idx)+1] #this is because you delete -np.inf so +1, and then because of : need one more+1
    return(list(Dt)+[np.inf],Ct)

def fi_to_gj(Dx,Cx,f0,q):
    '''
    Parameters
    ----------
    Dx : List
        Breakpoints in f(x).
    Cx : np.array(3,)
        Coefficents of x.

    Returns
    -------
    Dy :List
        Breakpoints in g(y)
    Cy : np.array(3,)
        Coefficents of y.
    '''
    # f0=np.min(Cx[2,:])
    Dy,Index=Project_Breakpoints_Robust(copy.deepcopy(Cx),copy.deepcopy(Dx))
    # Dy,Index=Project_Breakpoints(Cx,Dx)
    Dy2=list(np.array(Dy)/-q) #because the above function assume q=-1
    Cy=Conjugate_coeff(Cx[:,Index],q)
    Dy3,Cy=Add_indicator3(Dy2, Cy, f0)
    return(Dy3,Cy)


def solve_quad_array(C,Y):
    #Y is the set of breakpoints
    return(C[0]*np.square(Y)+C[1]*Y+C[2])

'''
Dy,Index=Project_Breakpoints(C,Dx) #Get Breakpoints
Cy=Conjugate_coeff(C[:,Index],1) 
gC=Cy[:,1:] # Coefficents of gy without first quadratic
# Add - to gy?? Answer: 
gy=solve_quad_array(gC,Dy[1:-1]) # Get all the intercepts p^*(y) 
gy=np.concatenate(([-np.infty],gy,[-np.infty])) #add infty
'''

def Add_indicator(gy,Dy,f0,Cy):
    '''
    Parameters
    ----------
    gy : np.array
        g(y) for y in Dy.
    Dy : List
        Breakpoints of the function g(y).
    f0 : real number
        Value of f(0).
    Cy : np.array(3,n)
        Coefficents of quadratics in g(y).

    Returns
    -------
    D2y : Updated list of breakpoints Dy
    '''
    temp=np.where(gy<f0)[0]
    D2y=[Dy[i] for i in temp]
    temp2=np.where(gy>f0)[0] # Speed up this line
    if len(temp2)==0:
        return(Dy,Cy)
    else:
        idx,jdx=min(temp2)-1,max(temp2) #delete ith index so choose i-1 index. Delete jth index but still choose j
        # for idx
        ai,bi,di=Cy[:,idx]
        y1,y2=np.roots([ai,bi,di-f0])
        y=y1*(y1>=Dy[idx] and y1<=Dy[idx+1])+y2*(y2>=Dy[idx] and y2<=Dy[idx+1])
        D2y.insert(idx+1,y)
        # for jdx
        ai,bi,di=Cy[:,jdx]
        y1,y2=np.roots([ai,bi,di-f0])
        y=y1*(y1>=Dy[jdx] and y1<=Dy[jdx+1])+y2*(y2>=Dy[jdx] and y2<=Dy[jdx+1])
        D2y.insert(idx+2,y) #inset at idx+2 and not jdx because we removed elements between idx and jdx from Dy
        Cy=np.concatenate((Cy[:,:idx+1],np.array([[0,0,f0]]).T,Cy[:,jdx:]),1)
        return(D2y,Cy)

def Sum_g(D,C_list):
    #This function is to add g(x1),g(x2).. when x1,x2 are children of some breanching node
    D=copy.deepcopy(D)
    C_list=copy.deepcopy(C_list)
    n=len(D)
    Dl=[]
    l=[]
    t=[0]*n
    Ctemp=Add_quadratic_equations(t, C_list)
    for i in D:
        Dl=Dl+i[1:-1]
        l.append(len(i[1:-1]))
    u=np.cumsum(l)
    Dls=sorted(set(Dl)) # set is to remove duplicates

    # for t,i in enumerate(Dl):
    for i in Dls:
        # idx=Dl.index(i)
        idxes=np.where(Dl==i)[0]
        for idx in idxes:
            tdx=min(np.where(u>idx)[0])
            t[tdx]+=1
        # print(t)
        Ctemp=np.column_stack((Ctemp,Add_quadratic_equations(t,C_list)))
    Dls.append(np.inf)
    Dls.insert(0,-np.inf)
    return(Dls,Ctemp)

def Add_quadratic_equations(t,C_list):
    # t is the index array of length number of branches. 
    q=np.zeros(3)
    for idx,C in enumerate(C_list):
        q+=C[:,t[idx]]
    return(q)

def Conjugate_coeff(coeff_I,q=-1):
    a_i,b_i,d_i=coeff_I
    """
    try:
        a_star=(q**2)/(-4*a_i)
        b_star=q*b_i/(2*a_i)
        d_star=d_i+(b_i**2)/(-4*a_i)
        return(a_i,b_i,d_i)
    except:
        return(a_i,b_i,d_i)
    """
    a_star=(np.square(q))/(-4*a_i)
    b_star=q*b_i/(-2*a_i)
    d_star=d_i+(b_i**2)/(-4*a_i)
    if len(np.shape(coeff_I))==1:
        return(a_star,b_star,d_star)
    else:
        return(np.stack((a_star,b_star,d_star)))

def Get_Slopes(coeff_I,coeff_J):
    a_i,b_i,d_i=coeff_I
    a_j,b_j,d_j=coeff_J
    '''
    try:
        A=0.25*(1/a_i-1/a_j)
        B=-0.5*(b_i/a_i-b_j/a_j)
        C=0.25*(b_i**2/a_i-b_j**2/a_j)-d_i+d_j    
        return(np.roots([A,B,C]))      
    except:
        return(np.roots([a_i,b_i,d_i-d_j]))
    '''
    A=0.25*(1/a_i-1/a_j)
    B=-0.5*(b_i/a_i-b_j/a_j)
    D=0.25*(np.square(b_i)/a_i-np.square(b_j)/a_j)-d_i+d_j
    
    return(np.roots([A,B,D]))
    
def Slope_bounds(coeff_I,Dx):
    l,u=Dx
    a_i,b_i,d_i=coeff_I
    yl,yu=2*a_i*l+b_i,2*a_i*u+b_i
    return(np.minimum(yl,yu),np.maximum(yl,yu))

def Get_intercept(coeff_I,y1):

    a_i,b_i,d_i=coeff_I
    a_s,b_s,d_s=Conjugate_coeff(coeff_I,-1)
    
    g_star1=a_s*y1**2+b_s*y1+d_s
    x1=np.roots([a_i,b_i-y1,d_i-g_star1])[0]
    
    return(x1)

def solve_quad(coeff_I,x):
    a_i,b_i,d_i=coeff_I
    return(a_i*x**2+b_i*x+d_i)


def g_sum_to_fj(Cy,Qjj,cj,lamj):
    Cy=Cy+np.array([[1/2*Qjj],[cj],[lamj]])
    return(Cy)


def Graph_Hierarchy(Q,r=0):

    G1=nx.Graph(Q)
    G1.remove_edges_from(nx.selfloop_edges(G1))
    a=np.array(G1.degree)
    # a=np.count_nonzero(Q,1)-1
    # d_nodes=np.where(a>2)[-1]
    d_nodes=a[np.where(a[:,1]>2)[-1],0]
    rank = {key: -1 for key in d_nodes}
    rank[r]=-1
    
    level=0
    breakflag=0
    ccl_Dict={key: [] for key in d_nodes}
    ccl_Dict[r]=[]
    while list(G1.nodes)!=[r]:
        a=np.array(G1.degree)
        leaves=a[np.where(a[:,1]==1)[-1],0]
        leaves=np.delete(leaves,np.where(leaves == r))
        trim=[]
        for i in leaves:
            parent=i
            son=None #this to get rid of an error
            trim_temp=[]
            #trim.append(parent)
            #set the base level for parent
            if parent not in rank.keys():
                base_level=0
            else:
                base_level=rank[parent]
            while parent != r:
                nbhd=list(G1.neighbors(parent)) #del parent from np.where(Q[parent,:]!=0)[0]
                if len(nbhd)!=1:
                    nbhd.remove(son)
                for n in nbhd:
                    if n in rank.keys():
                        temp=level
                        rank[n]=np.max(temp)
                        #print("Node ",n," has level ", level)
                        breakflag=1
                        #trim_temp.append(n)
                        ccl_Dict[n].append(trim_temp)
                        break          
    
                son=parent
                parent=nbhd[0]
                #print("Will trim ",son)
                trim.append(son)
                trim_temp.append(son)
                if breakflag==1:
                    breakflag=0
                    break
    
        #print("\n\n\n",trim)
        G1.remove_nodes_from(trim)
        #plt.figure()
        #nx.draw(G1, pos=pos, with_labels=True)
        level+=1
    
    #Add root nodes
    for i in ccl_Dict.keys():
        for j in ccl_Dict[i]:
            j.append(i)
    return(ccl_Dict,rank)
    
# def Add_indicator(gy,Dy,f0,C):
#     temp=np.where(gy<f0)[0]
#     D2y=[Dy[i] for i in temp]
#     idx,jdx=min(temp),max(temp)
#     # for idx
#     ai,bi,di=C[:,idx]
#     y1,y2=np.roots(ai,bi,di-f0)
#     y=y1*(y1>=Dy[idx-1] and y1<=Dy[idx])+y2*(y2>=Dy[idx-1] and y2<=Dy[idx])
#     D2y.insert(idx,y)
#     # for jdx
#     ai,bi,di=C[:,jdx]
#     y1,y2=np.roots(ai,bi,di-f0)
#     y=y1*(y1>=Dy[idx-1] and y1<=Dy[idx])+y2*(y2>=Dy[idx-1] and y2<=Dy[idx])
#     D2y.insert(idx+1,y) #inset at idx+1 and not jdx because we removed elements between idx and jdx from Dy
#     return(D2y)