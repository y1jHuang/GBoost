
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GBoost

<!-- badges: start -->

<!-- badges: end -->

Group boosting is an algorithm developed to perform variable selection
in a scalar-on-network regression setting, where the network has
underlying group structure.

## Installation

To install the updated package `GBoost`, please enter R's console, and type in:

```R
install.packages('#path/GBoost', repos=NULL, type='source')
```

replace `#path` with the folder path where the downloaded source code of `GBoost` exists.

Installation error might occur for `GroupBoosting` package. It requires Intel TBB 2020 to build up the package. If your operating system is centOS, then you might have to update the version of TBB and edit the `Makevar` file, whose path generally is `GroupBoosting >> src >> Makevars`. So do other packages. Here is my command, and you need to replace the TBB path with your own in `Makevars` file:

`PKG_LIBS=$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -L/data/app/intel2020.4/tbb/lib/intel64/gcc4.8 -ltbb`

However, if your operating system is Windows, then nothing should be done before installing the package. And the installation might be successful without any errors.

## Example

To avoid masking functions from `GroupBoosting` and for convenience while debugging, I renamed the package with `GBoost` and modified names of its functions. There are also some slight change to call them. Here I would like to use an example for further illustration.

``` r
library(GBoost)
data1 = simul_group_data(nodes = 6, n = 100, 
                         num.groups = 2, q.groups = 1, 
                         sparse_g = 0, dense_g = 1, 
                         effect_size = 5)
test1 = GBoost_fit(data1$X, data1$Y, data1$group, 
                   total_steps=5000, step_size=c(1e-2,1e-2),
                   adj_var = 999, stop_tol=-1e-5, gamma = 0, 
                   lasso_lambda = 0.0314, weighted = 'n')
cbind(test1$beta, data1$beta)

          [,1] [,2]
 [1,] 4.833239    5
 [2,] 5.305087    5
 [3,] 0.000000    0
 [4,] 0.000000    0
 [5,] 0.000000    0
 [6,] 5.173285    5
 [7,] 0.000000    0
 [8,] 0.000000    0
 [9,] 0.000000    0
[10,] 0.000000    0
[11,] 0.000000    0
[12,] 0.000000    0
[13,] 0.000000    0
[14,] 0.000000    0
[15,] 0.000000    0
```

First, `step_size` should be entered with a vector, which is the step size $v$ for group and edge selection, respectively.

Second, If not use weighted version to summarize within groups for the first stage, set `weighted = 'n'`. Otherwise, set `weighted = 'lasso'` or `weighted = 'boosting'` to specify the method for a weighted average of the edges.