# Efficient Inference of Dynamic Gene Regulatory Networks Using Discrete Penalty

This repository contains the code for the algorithm described in the paper. The method supports both single-category and two-category joint network inference using a discrete ℓ₀ penalty and structural similarity constraints. It enables efficient joint estimation of sparse precision matrices across multiple populations, with support for separate inference in two categories (e.g., primary vs. recurrent disease).

---

## Usage

To run the method, import the appropriate function depending on the number of categories:

```python
from ELEM_0 import elem_0
from Cat_ELEM_0 import elem_0_cat
```

---

## Function: `elem_0`

```python
Theta = elem_0(data, W, nu0, mu, gamma)
```

This function performs joint estimation of multiple precision matrices (from a single category) under an ℓ₀ sparsity penalty and a similarity constraint.

### Arguments:

- **`data`**: `dict`  
  Gene expression data 
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

---

## Function: `elem_0_cat`

```python
Theta = elem_0_cat(data1,data2, W, nu0, mu, gamma)
```

This function performs joint estimation of networks for two distinct categories of populations (e.g., primary vs. recurrent disease), enforcing sparsity within each category and structural similarity across populations within each category.

### Arguments:

- **`data1`**: `dict`  
  Gene expression data for **category 1**.    
  A dictionary where each key is a population index (`0, 1, ..., K-1`),  
  and each value is an `(n_k, p)` NumPy array of gene expression data.

- **`data2`**: `dict`  
  Gene expression data for **category 2**.  
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

- **`Theta`**: `np.ndarray` of shape `(p, p, 3K)`  
  A stacked set of estimated precision matrices:
  
  - `Theta[:, :, 0:K]`: estimated **global (shared)** networks  
  - `Theta[:, :, K:2K]`: estimated networks for **category 1**  
  - `Theta[:, :, 2K:3K]`: estimated networks for **category 2**
