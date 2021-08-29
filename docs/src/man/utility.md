# [Utility Functions](@id utilityfunction)

The GaussBP provides several utility functions to evaluate and compare obtained results.

---
#### The WLS results
The function provides the estimate obtained by the WLS method and root mean square error (RMSE), the mean absolute error (MAE) and the weighted residual sum of squares (WRSS) error metrics evaluated according to the WLS solutions. These results can be used to compare results obtained by the GBP algorithm.
```julia-repl
exact = wls(gbp)
```
The function returns the composite type `WeightedLeastSquares` with fields `estimate`, `rmse`, `mae`, `wrss`. Note that results are obtained according to variables `SystemModel.jacobian`, `SystemModel.observation` and `SystemModel.variance`.

----
#### The GBP error metrics
The package provides the function to obtain RMSE, MAE, and WRSS error metrics of the GBP algorithm.
```julia-repl
evaluation = errorMetric(gbp)
```
The function returns the composite type `ErrorMetric` with fields `rmse`, `mae`, `wrss`. Further, passing the composite type `WeightedLeastSquares`, we obtained additional fields `rmseGBPWLS` and `maeGBPWLS` that determine the distance between the GBP estimate and WLS estimate.
```julia-repl
evaluation = errorMetric(gbp, exact)
```
The function returns the composite type `ErrorMetricWiden` with fields `rmse`, `mae`, `wrss`, `rmseGBPWLS`, `maeGBPWLS`.

----
#### Display results
The function shows data display in the Julia REPL, and can provide different views depending on the input variables.

The following function can be used to show GBP results:
```julia-repl
displayData(gbp)
```

To show the GBP results and error metric use:
```julia-repl
displayData(gbp, evaluation)
```

To show the GBP and WLS results use:
```julia-repl
displayData(gbp, exact)
```

To show the GBP and WLS results, and error metrics use:
```julia-repl
displayData(gbp, exact, evaluation)
```

---

#### Error metrics

The root mean square error, the mean absolute error and the weighted residual sum of squares are evaluated according to:
```math
  \begin{aligned}
    \text{rmse} = \sqrt {\cfrac{\sum_{i=1}^m \left[z_i - h_i(\hat{\mathbf x}) \right]^2}{m}}; \quad
    \text{mae} = \cfrac{\sum_{i=1}^m \left|z_i - h_i(\hat{\mathbf x}) \right|}{m}; \quad
    \text{wrss} = \sum_{i=1}^m \cfrac{\left[z_i - h_i(\hat{\mathbf x}) \right]^2}{v_i},
  \end{aligned}
```
where ``m`` denotes the number of observations, ``z_i`` is observation value, ``v_i`` is observation variance, and corresponding equation ``h_i(\hat{\mathbf x})`` is evaluated at the point ``\hat{\mathbf x}`` obtained using the GBP or WLS algorithm. Note, `wrss` is the value of the objective function of the optimization problem we are solving.

```@raw html
&nbsp;
```

Fields `rmseGBPWLS` and `maeGBPWLS` determine distance beetwen the GBP estimate ``\hat{x}_{\text{gbp},i}`` and WLS estimate ``\hat{x}_{\text{wls},i}``, where root mean square error and mean absolute error are obtained using:
```math
  \begin{aligned}
    \text{rmse} = \sqrt {\cfrac{\sum_{i=1}^n \left[\hat{x}_{\text{wls},i} - \hat{x}_{{\text{gbp}},i}) \right]^2}{n}}; \quad
    \text{mae} = {\cfrac{\sum_{i=1}^n \left|\hat{x}_{\text{wls},i} - \hat{x}_{{\text{gbp}},i}) \right|}{n}},
  \end{aligned}
```
where ``n`` is the number of state variables.
