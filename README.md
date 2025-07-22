# Efficient inference of dynamic gene regulatory networks using discrete penalty

This repository contains the code of the algorithm in the paper "Efficient inference of dynamic gene regulatory networks using l_0 penalty".


Use the line below to import the function:
```
from elem_0 import elem_0
```

The function *Para_Algo* requires 5 arguments.

*data*: dict.
        A dictionary where each key corresponds to a population index (0, 1, ..., K-1),  
        and each value is an (n_k,p) np.array representing the gene expression data for that category.

        
   *W*: np.ndarray size (K,K).
        The symmetric matrix representing the hypergraph structure
        
 *nu0*: float.
        Soft-thresholding parameter used in the approximate backward mapping.
        
  *mu*: float.
        Penalty weight for the sparsity-inducing \(\ell_0\) norm. Larger values encourage sparser solutions.
        
*gamma*: float.
        Penalty weight for enforcing similarity between the precision matrices of neighboring categories
        as defined by the hypergraph W.

