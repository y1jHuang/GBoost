# `GBoost`

## Introduction

`GBoost` is a package implements variable screening via $L_2$ boosting algorithm in a scalar-on-network regression setting, where the network has underlying group structure. Compared to other linear methods, the algorithm is relatively fast with high performance, and can exploit sparsity in the input matrix `X`.

The authors of `GBoost` are Emily L. Morris and Yingjie Huang, who are maintainers as well.

## Installation

#### prerequisites

Imports: `Rcpp`, `RcppArmadillo`, `RcppParallel`

To successfully compile `.cpp` files in the package, Intel&reg; `Threading Building Blocks(TBB)` is required. Generally, for `Windows` users, it has already been installed in the system. However, `GNU/Linux` users might have to install Intel&reg; `TBB` on their own, and manually edit `Makevar` file, whose path is `GBoost >> src >> Makevar`:

```bash
PKG_LIBS=$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -L/data/app/intel/tbb/lib/intel64/gcc4.8 -ltbb
```

Replace `/data/app/intel/tbb/lib/intel64/gcc4.8` with your own Intel&reg; `TBB` path.

#### Installation

Please type the following command in R console:

```R
install.packages('#path/GBoost', repos=NULL, type='source')
```

Replace `#path` with the folder path where downloaded `GBoost` package located.

## Quick Start

```R
library(GBoost)
library(glmnet)
# Simulate data with the assigned effect
data = simul_group_data(nodes = 6, n = 100, 
                        num.groups = 2, q.groups = 1,
                        sparse_g = 0, dense_g = 1, 
                        effect_size = 5)
# Estimation by GBoost
mdl.GBoost = GBoost_fit(data$X, data$Y, data$group, 
                        total_steps=5000, step_size=c(1e-2,1e-2), 
                        adj_var = 999, stop_tol=-1e-5, gamma = 1, 
                        lasso_lambda = 0.0314, weighted = 'n')
beta.GBoost = mdl.GBoost$beta

# Estimation by Lasso
mdl.lasso <- cv.glmnet(data$X, data$Y, alpha=1, parallel = FALSE)
beta.lasso = as.vector(stats::coef(mdl.lasso, s="lambda.min"))[-1]

table = cbind(beta.GBoost, beta.lasso, data$beta, data$group)
colnames(table) = c('GBoost', 'Lasso', 'effect','group')

> table
        GBoost       Lasso     effect     group
 [1,] 4.829109  4.53929582          5         1
 [2,] 5.302042  5.05273894          5         1
 [3,] 0.000000  0.00000000          0         2
 [4,] 0.000000 -0.09168638          0         2
 [5,] 0.000000 -0.32576087          0         2
 [6,] 5.168462  4.97995221          5         1
 [7,] 0.000000 -0.18823506          0         2
 [8,] 0.000000  0.00000000          0         2
 [9,] 0.000000  0.13367415          0         2
[10,] 0.000000 -0.07816475          0         2
[11,] 0.000000  0.00000000          0         2
[12,] 0.000000  0.00000000          0         2
[13,] 0.000000 -0.35384280          0         3
[14,] 0.000000  0.00000000          0         3
[15,] 0.000000  0.00000000          0         3
```

Here is an example of `GBoost`. Using `simul_group_data()` function, we first simulate the data with assigned effect:


$$
\pmb y = \mathbf X \pmb \omega + \pmb \epsilon,~
\pmb \omega = 
\begin{pmatrix}
5&5&0&0&0&5&0&\cdots&0
\end{pmatrix}^\intercal
$$


where $\epsilon_i \sim \mathcal{N}(0, 1), \ x_{i,j} \sim \mathcal{N}(\pmb{\mu},1), \ \mu_{i,j} \sim \rm{Ber}(\lbrace -1, 1 \rbrace; \theta)$

Then we use `GBoost_fit` for estimation. 

`data$X` is the input `X` which contains both adjustment variables and predictors.   
`data$Y` is the response variable `Y`.   
`data$group` is a factor indicates the group structure of predictors.  
`totoal_steps` is the maximum number for iteration.  
`step_size` is step size $v$. Generally, it is better to use small step size, while takes longer time to process.  
`adj_var` is a vector indexes the column of adjustment variables. If there are no adjustment variables, set it to `999 `  
`stop_tol` is a minimum number of AIC/BIC's derivative. Iteration will not stop until reaches it.  
`gamma` is   
`weighted` is an option to process $\mathbf{X}$. If `'n'`, we normally average edges within each group; if `'lasso'`, we use Lasso algorithm to identify informative edges, and perform weighted average. `'boosting'` is similar.

## Algorithm

**Data:** $\lbrace \mathbf x_i, y_i \rbrace ^n_{i=1};$ number of iterations $M_1$ and $M_2$; updating rate $v$  
**Result:** $\lbrace\lambda_g\rbrace$ and $\lbrace\omega_{j,k}\rbrace$  
**begin**  
&emsp; Initialize $\lambda_g = 0 $, $g = 1,...,G$  
&emsp; Calculate the average of edges in sub-network $g$  
&emsp; **for** $g=1,...,G$ **do**  
&emsp; &emsp; $\overline{\mathbf X }^{(g)} = m_g^{-1} \sum_{(j,k)\in I_g} \mathbf X^{(j,k)}$  
&emsp; **end**  
&emsp; Stage 1 of $L_2$ boosting algorithm:  
&emsp; **for** $m = 1,...,M_1$ **do**  
&emsp;&emsp; **for** $g = 1,...,G$ **do**  
&emsp;&emsp;&emsp; Compute the first partial derivative with respect to $\lambda_g$  
&emsp;&emsp;&emsp; $L_1(g) = (\pmb{y}-\overline{\mathbf{X}}\pmb{\lambda} )^T \overline{\mathbf{X}}^{(g)} $  
&emsp;&emsp; **end**  
&emsp;&emsp; Find $\hat{g} = \underset{1 \leq g \leq G}{\arg\max} \left\vert L_1(g) \right\vert $  
&emsp;&emsp; Calculate the second partial derivative with respect to $\hat{g}$  
&emsp;&emsp; $L_2(\hat{g}) = \left\Vert \overline{\mathbf{X}}^{(\hat{g})} \right\Vert ^2$  
&emsp;&emsp; Update $\lambda_{\hat{g}} = \lambda_{\hat{g}} + v L_2(\hat{g})^{-1}L_1(\hat{g})$  
&emsp; **end**  
&emsp; Initialize $\omega_{j,k} = 0, (j, k) \in S = \bigcup_{g:\lambda_g \neq 0} I_g$  
&emsp; Filter out edges belong to irrelevant networks  
&emsp; $\mathbf{X} = \lbrace\mathbf{X}^{(j,k)}\rbrace, (j, k) \in S = \bigcup_{g':\lambda_{g'} \neq 0} I_g$  
&emsp; Stage 2 of $L_2$ boosting algorithm:  
&emsp; **for** $m = 1,...,M_2$ **do**  
&emsp;&emsp; **for** $(j,k)\in S$ **do**  
&emsp;&emsp;&emsp; Compute the first partial derivative with respect to ???????  
&emsp;&emsp;&emsp; $L_1(j, k) = (\pmb{y} - \mathbf{X}\pmb{\omega})^T \mathbf{X}$  
&emsp;&emsp; **end**  
&emsp;&emsp; Find $(\hat{j}, \hat{k}) = \underset{(j,k)\in S}{\arg\max} \left\vert L_1(j,k) \right\vert$  
&emsp;&emsp; Calculate the second partial derivative with respect to $(\hat{j}, \hat{k})$  
&emsp;&emsp; $L_2(\hat{j}, \hat{k}) = \left\Vert \mathbf{X}^{(\hat{j}, \hat{k})} \right\Vert^2$  
&emsp;&emsp; Update $\omega_{\hat{j}, \hat{k}} = \omega_{\hat{j}, \hat{k}} + v L_2(\hat{j}, \hat{k})^{-1}L_1(\hat{j}, \hat{k})$  
&emsp; **end**  
**end**

## Functions

