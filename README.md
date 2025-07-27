# Efficient Inference of Dynamic Gene Regulatory Networks Using Discrete Penalty

This repository contains the code for the algorithm described in the paper. 


---

## Usage

To run the method, import the main function:

```python
from elem_0 import elem_0
```

---

## Function: `elem_0`

```python
Theta = elem_0(data, W, nu0, mu, gamma)
```

This function performs joint estimation of multiple precision matrices (1 category) 
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

