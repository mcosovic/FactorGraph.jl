# [Output Data](@id outputdata)

The function `gbp()` returns a struct variable `results` with fields `means`, `variances` and `iterations` containing the GBP algorithm results, and the additional fields `rmse`, `mae` and `wrss` if used command `out = "error"`. 

Also, if used command `out = "wls"`, the struct variable `results` contains additional fields `meansWLS`, `rmseWLS`, `maeWLS`, `wrssWLS`, `rmseBPWLS`, `maeBPWLS`.  

Finally, a struct variable `system` contains the input data.

---

#### The GBP algorithm results
The vectors `means` and `variances` contains means and variances of the marginals obtained using the GBP algorithm, while `iterations` save the maximum number of iterations. Using the keyword `out = "iterate"`, fields `means` and `variances` are given as matrices, where each column contains mean and variance values in the corresponding iteration according to the vector `iterations`.     

```@raw html
&nbsp;
```

#### The GBP error metrics
Using the command `out = "error"`, we obtained root mean square error `rmse`, mean absolute error `mae`, and weighted residual sum of squares `wrss` of the GBP algorithm at the converged point. Further, fields `rmse`, `mae` and `wrss` become vectors if `out = ["iterate", "error"]` command is used, where each error value corresponds to the iteration according to the vector `iterations`.

The root mean square error is obtained as:
```math
  \begin{aligned}
    \text{rmse} = \sqrt {\cfrac{\sum_{i=1}^m \left[z_i - h_i(\hat{\mathbf x}) \right]^2}{m}},
  \end{aligned}
```  
where ``m`` denotes the number of observations, ``z_i`` is observation value, and corresponding equation ``h_i(\hat{\mathbf x})`` is evaluated at the point ``\hat{\mathbf x}`` obtained using the GBP throughout iterations or at the final iteration.  

The mean absolute error is obtained as:
```math
  \begin{aligned}
    \text{mae} = \cfrac{\sum_{i=1}^m \left|z_i - h_i(\hat{\mathbf x}) \right|}{m}.
  \end{aligned}
```  

The weighted residual sum of squares is obtained as:
```math
  \begin{aligned}
    \text{wrss} = \sum_{i=1}^m \cfrac{\left[z_i - h_i(\hat{\mathbf x}) \right]^2}{v_i},
  \end{aligned}
```  
where ``v_i`` is observation variance. Note that `wrss` is the value of the objective function of the optimization problem we are solving.

```@raw html
&nbsp;
```

#### The WLS and GBP error metrics
Using command `out = "wls"`, we obtained error metrics `rmseWLS`, `maeWLS` and `wrssWLS`, evaluated according to above equations, where ``\hat{\mathbf x}`` is obtained by solving WLS problem. Fields `rmseBPWLS` and `maeBPWLS` determine distance beetwen the GBP estimate ``\hat{x}_{i,{\text{gbp}}}`` and WLS estimate ``\hat{x}_{i,{\text{wls}}}``, where root mean square is obtainde using:
```math
  \begin{aligned}
    \text{rmse} = \sqrt {\cfrac{\sum_{i=1}^n \left[\hat{x}_{i,\text{wls}} - \hat{x}_{i,{\text{gbp}}}) \right]^2}{n}},
  \end{aligned}
```  
where ``n`` is the number of state variabels. The mean absolute error is computed as:
```math
  \begin{aligned}
    \text{mae} = {\cfrac{\sum_{i=1}^n \left|\hat{x}_{i,\text{wls}} - \hat{x}_{i,{\text{gbp}}}) \right|}{n}}.
  \end{aligned}
```  
Note that error metrics `rmseBPWLS` and `maeBPWLS` can be given also for each GBP iteration, if `out = ["iterate", "wls"]` command is used.

