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



        # Efficient Inference of Dynamic Gene Regulatory Networks Using Discrete Penalty

This repository contains the code for the algorithm described in the paper  
**"Efficient Inference of Dynamic Gene Regulatory Networks Using ℓ₀ Penalty."**

---

## Usage

To run the method, import the main function:

```python
from elem_0 import elem_0
```

This returns a handle to the core optimization routine `Para_Algo`.

---

## Function: `Para_Algo`

```python
Theta, history = Para_Algo(data, W, nu0, mu, gamma)
```

This function performs joint estimation of multiple precision matrices  
under an ℓ₀ sparsity penalty and a similarity constraint.

### Arguments:

- **`data`**: `dict`  
  A dictionary where each key is a population index (`0, 1, ..., K-1`),  
  and each value is an `(n_k, p)` NumPy array of gene expression data.

- **`W`**: `np.ndarray` of shape `(K, K)`  
  A symmetric matrix representing the hypergraph (adjacency structure) over the populations.

- **`nu0`**: `float`  
  Soft-thresholding parameter used in the approximate backward mapping step.

- **`mu`**: `float`  
  Penalty weight for the ℓ₀ norm. Larger values encourage sparser solutions.

- **`gamma`**: `float`  
  Penalty weight that enforces similarity between the precision matrices  
  of neighboring populations, as defined by `W`.

### Returns:

- **`Theta`**: `np.ndarray` of shape `(p, p, K)`  
  The estimated precision matrices for each population.

