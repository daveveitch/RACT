# RACT
 
**Rank-adaptive covariance testing**

R package to apply RACT method for two-sample covariance testing.

### References ### 

> Veitch, D., He, Y., & Park, J. Y. (2024). Rank-adaptive covariance testing with applications to genomics and neuroimaging. arXiv preprint arXiv:2309.10284. [link](https://arxiv.org/abs/2309.10284)

## Contents

1. [Installation](#id-installation)
2. [Usage](#id-usage)

<div id='id-installation'/>

---

### Installation
To install the latest development builds directly from GitHub, please run the following:

```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("daveveitch/RACT")
```

<div id='id-usage'/>

---

### Usage

This package contains one main function `RACT()`, and one helper function `K_calculate()`.

`RACT()` is the main function which takes two data matrices and tests the null hypotheses $H_0: \Sigma_1 = \Sigma_2$. This function takes the following arguments:

* `X_1`: a ($p$ x $n_1$) data matrix, where $p$ is the number of features and $n_1$ is the number of subjects in group 1.
* `X_2`: a ($p$ x $n_2$) data matrix, where $p$ is the number of features and $n_2$ is the number of subjects in group 2.
* `n_perm`: (optional, default is `n_perm=1000`) Number of permutations to use in the permutation procedure.
* `K`: (optional, default is `K=NULL`) If `K=NULL`, Ky-Fan($k$) norms from 1 to $K$ are included where $K$ the smallest integer such that sum of top $K$ singular values of the pooled covariance are greater than 80% of the sum of all of the singular values. Otherwise this argument takes a vector of integers, and these integers represent the Ky-Fan(k) norms that are included. For example, `K=c(1,5)` means RACT is based on the Ky-Fan(1) and Ky-Fan(5) norms only.
* `min_P`: (optional, default is `min_P=FALSE`) Whether a minimum $p$-value statistic is used.
* `cov`: (optional, default is `cov=TRUE`) If `cov=TRUE` RACT tests for the equality of the covariance of `X_1`, `X_2`. If `cov=FALSE` RACT tests for the equality of the correlation of `X_1`, `X_2`.

`RACT()` will return a list with two entries, 'RACT p value' the p-value from the adaptive testing procedure, and 'Individual Ky-Fan(k) p values' the p-values from each individual Ky-Fan($k$) norm.

```R
RACT_p_values = RACT(X_1,X_2,n_perm=1000,K=NULL,min_P=FALSE,cov=TRUE)
``` 

---

`K_calculate()` is a helper function which calculates $K$ (a vector representing the Ky-Fan($k$) to include in RACT) based on the pooled covariance or pooled correlation. The function takes the following arguments

* `X_1`: a ($p$ x $n_1$) data matrix, where $p$ is the number of features and $n_1$ is the number of subjects in group 1.
* `X_2`: a ($p$ x $n_2$) data matrix, where $p$ is the number of features and $n_2$ is the number of subjects in group 2.
* `K_pct`: (optional, default is `K_pct = 0.8`) 0 < `K_pct` <= 1, where `K_pct` determines how large $K$ should be so that sum of top $K$ singular values of pooled covariance/correlation at least `K_pct` of the sum of all singular values.
* `cov`: (optional, default is `cov=TRUE`) If `cov=TRUE` RACT tests for the equality of the covariance of `X_1`, `X_2`. If `cov=FALSE` RACT tests for the equality of the correlation of `X_1`, `X_2`.

```R
K = K_calculate(X_1,X_2,K_pct = 0.8,cov = TRUE)
``` 

`K_calculate()` will return a vector of values c(1,...,K) which can then be used as an argument for `RACT()`.

---

