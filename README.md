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
data = simul_group_data(nodes = 6, n = 100, 
                        num.groups = 2, q.groups = 1,
                        sparse_g = 0, dense_g = 1, 
                        effect_size = 5)
mdl = GBoost_fit(data$X, data$Y, data$group, 
                 total_steps=5000, step_size=c(1e-2,1e-2), 
                 adj_var = 999, stop_tol=-1e-5, gamma = 1, 
                 lasso_lambda = 0.0314, weighted = 'n')
table = cbind(mdl$beta, data$beta, data$group)
colnames(table) = c('GBoost','effect','group')

> table
        GBoost    effect    group
 [1,] 4.829109         5        1
 [2,] 5.302042         5        1
 [3,] 0.000000         0        2
 [4,] 0.000000         0        2
 [5,] 0.000000         0        2
 [6,] 5.168462         5        1
 [7,] 0.000000         0        2
 [8,] 0.000000         0        2
 [9,] 0.000000         0        2
[10,] 0.000000         0        2
[11,] 0.000000         0        2
[12,] 0.000000         0        2
[13,] 0.000000         0        3
[14,] 0.000000         0        3
[15,] 0.000000         0        3

```

Here is an example of `GBoost`. Using `simul_group_data()` function, we first simulate the data with assigned $\begin{matrix}
\pmb{\omega} = (5&5&0&0&0&5&0&...&0)^T
\end{matrix}$ 
$$
\pmb{y} = \pmb{X}\pmb{\omega} + \pmb{\epsilon}
$$
where $\epsilon_i \sim \mathcal{N}(0, 1), x_{i,j} \sim \mathcal{N}(\pmb{\mu},1), \ \mu_{i,j} \sim \rm{Ber}(x; 0.3),\
\rm{Ber}(x;\theta) = \left\{ \begin{array}{rcl} 
\theta & \mbox{if} & x = 1\\
1-\theta & \mbox{if} & x = -1\\
\end{array}\right.$

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

**Data:** $\{ \mathbf{x}_i, y_i \} ^n_{i=1};$ number of iterations $M_1$ and $M_2$; updating rate $v$  
**Result:** $\{\lambda_g\}$ and $\{\omega_{j,k}\}$  
**begin**  
	Initialize $\lambda_g = 0 $, $g = 1,...,G$  
	Calculate the average of edges in sub-network $g$  
​	**for** $g=1,...,G$ **do**  
​		$\overline{\mathbf{X}}^{(g)} = m_g^{-1} \sum_{(j,k)\in I_g}\mathbf{X}^{(j,k)}$  
​	**end**  
​	Stage 1 of $L_2$ boosting algorithm:  
​	**for** $m = 1,...,M_1$ **do**  
​		**for** $g = 1,...,G$ **do**  
​			Compute the first partial derivative with respect to $\lambda_g$  
​			$L_1(g) = (\pmb{y}-\overline{\mathbf{X}}\pmb{\lambda} )^T \overline{\mathbf{X}}^{(g)} $  
​		**end**  
​		Find $g^* = \underset{1 \leq g \leq G}{\arg\max} \vert L_1(g) \vert $  
​		Calculate the second partial derivative with respect to $g^*$  
​		$L_2(g^*) = \parallel \overline{\mathbf{X}}^{(g^*)} \parallel ^2$  
​		Update $\lambda_{g^*} = \lambda_{g^*} + v L_2(g^*)^{-1}L_1(g^*)$  
​	**end**  
​	Initialize $\omega_{j,k} = 0, (j, k) \in S = \bigcup_{g:\lambda_g \neq 0} I_g$  
​	Filter out edges belong to irrelevant networks  
​	$\mathbf{X} = \{\mathbf{X}^{(j,k)}\}, (j, k) \in S = \bigcup_{g':\lambda_{g'} \neq 0} I_g$  
​	Stage 2 of $L_2$ boosting algorithm:  
​	**for** $m = 1,...,M_2$ **do**  
​		**for** $(j,k)\in S$ **do**  
​			Compute the first partial derivative with respect to ???????  
​			$L_1(j, k) = (\pmb{y} - \mathbf{X}\pmb{\omega})^T \mathbf{X}$  
​		**end**  
​		Find $(j^*, k^*) = \underset{(j,k)\in S}{\arg\max} \vert L_1(j,k) \vert$  
​		Calculate the second partial derivative with respect to $(j^*,k^*)$  
​		$L_2(j^*,k^*) = \parallel \mathbf{X}^{(j^*,k^*)} \parallel ^2$  
​		Update $\omega_{j^*, k^*} = \omega_{j^*, k^*} + v L_2(j^*, k^*)^{-1}L_1(j^*,k^*)$  
​	**end**  
**end**  

