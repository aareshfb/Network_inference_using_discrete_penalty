# Efficient Inference of Dynamic Gene Regulatory Networks Using Discrete Penalty

This repository contains the code for the algorithm described in the following paper:

**Efficient inference of dynamic gene regulatory networks using discrete penalty.**  
arXiv preprint: https://arxiv.org/pdf/2507.23106

 The method supports both single-category and multiple-category joint network inference using a discrete ℓ₀ penalty and structural similarity constraints. It enables efficient joint estimation of sparse precision matrices across multiple populations, with support for separate inference in two or more categories (e.g., primary vs. recurrent disease).


---

## Usage

To run the method, import the appropriate function depending on the number of categories:

```python
from ELEM_0 import elem_0
from Mult_Cat_ELEM_0 import elem_0_catM
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
Theta = elem_0_catM(data_list, W, nu0, mu, gamma)
```

This function performs joint estimation of networks for two distinct categories of populations (e.g., primary vs. recurrent disease), enforcing sparsity within each category and structural similarity across populations within each category.

### Arguments:

- **`data_list`**: `list`  
  A list of gene expression datasets for different categories. The length of the list corresponds to the number of categories.

  Each element of `data_list` is a dictionary containing the gene expression data for a single category. Let `data_c` be an element of `data_list`, then:

    - **`data_c`**: `dict`  
      Gene expression data for **category c**.    
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

- **`Theta`**: `np.ndarray` of shape `(p, p, (C+1)*K)`  
  Here `C` denotes the number of categories. `Theta` is stacked set of estimated precision matrices:
  
  - `Theta[:, :, 0:K]`: estimated **global (shared)** networks  
  - `Theta[:, :, K:2*K]`: estimated networks for **category 1**  
  - `Theta[:, :, 2*K:3*K]`: estimated networks for **category 2**
  - `Theta[:, :, C*K:(C+1)*K]`: estimated networks for **category C**
 
--

The scripts used for the experiments were originally implemented for the two-category setting (primary vs. recurrent disease). For completeness, we include them here. They can be executed as follows:

```python

from Cat_ELEM_0 import elem_0_cat
Theta = elem_0_catM(data_1,data_2, W, nu0, mu, gamma)
```
Here, `data_1` and `data_2` are dictionaries with the same structure as `data_c` described above.


---


## References

The following methods are used in this repository:

- **ELEM-1**  
  Ravikumar, Visweswaran, et al. "Efficient inference of spatially-varying Gaussian Markov random fields with applications in gene regulatory networks." IEEE/ACM transactions on computational biology and bioinformatics 20.5 (2023): 2920-2932.

- **FASJEM**  
  Wang, Beilun, Ji Gao, and Yanjun Qi. "A fast and scalable joint estimator for learning multiple related sparse Gaussian graphical models." Artificial Intelligence and Statistics. PMLR, 2017.
  GitHub: https://github.com/QData/FASJEM

- **FGL / GGL**  
  Danaher, Patrick, Pei Wang, and Daniela M. Witten. "The joint graphical lasso for inverse covariance estimation across multiple classes." Journal of the Royal Statistical Society Series B: Statistical Methodology 76.2 (2014): 373-397.

- **GRNBoost2**  
  Moerman, Thomas, et al. "GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks." Bioinformatics 35.12 (2019): 2159-2161.

- **SILGGM**  
  Zhang, Rong, Zhao Ren, and Wei Chen. "SILGGM: An extensive R package for efficient statistical inference in large-scale gene networks." PLoS computational biology 14.8 (2018): e1006369.

- **jointGHS**  
  Lingjærde, Camilla, et al. "Scalable multiple network inference with the joint graphical horseshoe." The Annals of Applied Statistics 18.3 (2024): 1899-1923.
  GitHub: https://github.com/Camiling/jointGHS

- **scMTNI benchmark / datasets**  
  Zhang, Shilu, et al. "Inference of cell type-specific gene regulatory networks on cell lineages from single cell omic datasets." Nature Communications 14.1 (2023): 3064.
  GitHub: https://github.com/Roy-lab/scMTNI
