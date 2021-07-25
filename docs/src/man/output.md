# [Output Data](@id outputdata)

The function `gbp()` returns a struct variable `results` with fields `gbp` and `wls`. 

The variable `results.gbp` contains the following fields:
- `mean`, `variance`, `iteration`, `rmse`, `mae`, `wrss` 
These variables save the GBP algorithm results. Values of the error metrics `rmse`, `mae` and `wrss` are available only if we use the command `out = evaluation`. 

The variable `results.wls` contains the following fields:
- `mean` `rmse`, `mae`, `wrss`, `rmseGBPWLS`, `maeGBPWLS` 
Those variables save the WLS results, which can be used to compare the results obtained by the GBP algorithm. However, values are only available if we use the command `out = wls`.

Finally, a struct variable `system` contains the input data.

---

## [The GBP Results] (@id gbpresults)
The struct variable `results.gbp` contains the GBP algorithm results. To describe the outputs, we will use the example shown below. 
```julia-repl
using GaussBP

H = [1.5 2.0; 
     3.1 4.6]
z = [0.8; 4.1]
v = [1.0; 1.0]  

d = [5 2 2.4  1.5;
     8 1 0.85 0.9]
```

```@raw html
&nbsp;
```

#### Compute means and variances 
The `mean` and `variance` vectors contain mean and variance values of the marginals obtained using the GBP algorithm after the algorithm reaches the maximum number of iterations, while variable `iteration` save the maximum number of iterations.

```julia-repl
results, ~ = gbp(H, z, v; algorithm = vanilla)

julia> results.gbp.mean
2-element Vector{Float64}:
 -5.109180124020112
  4.148396436639841

julia> results.gbp.variance
2-element Vector{Float64}:
 4035.0842444764703
 4033.696884882772

julia> results.gbp.iteration
30 
```

```@raw html
&nbsp;
```

#### Compute means, variances and error metrics
In addition, using the command `out = evaluation`, except `mean`, `variance` and `iteration` variables, we obtained root mean square error `rmse`, mean absolute error `mae`, and weighted residual sum of squares `wrss` of the GBP algorithm after the algorithm reaches the maximum number of iterations. 

```julia-repl
results, ~ = gbp(H, z, v; algorithm = vanilla, out = evaluation)

julia> results.gbp.rmse
1-element Vector{Float64}:
 0.616577078168593

julia> results.gbp.mae
1-element Vector{Float64}:
 0.511406044334784

julia> results.gbp.wrss
1-element Vector{Float64}:
 1.022812088669568
```

```@raw html
&nbsp;
```


#### Compute means and variances through GBP iterations
Using the keyword `out = iteration`, fields `mean` and `variance` are given as matrices, where each column contains mean and variance values in the corresponding iteration according to the vector `iteration`.   

```julia-repl
results, ~ = gbp(H, z, v; algorithm = vanilla, max = 5, out = iteration)

julia> results.gbp.mean
2×5 Matrix{Float64}:
 0.885904  -0.108862   0.643748  -0.413415  0.206589
 0.671831   0.0883903  0.822573   0.335673  1.0947

julia> results.gbp.variance
2×5 Matrix{Float64}:
 98361.6  48877.4  94041.0  45748.4  86286.4
 25127.9  48877.3  24024.2  45748.0  22043.2

julia> results.gbp.iteration
5-element Vector{Int64}:
 1
 2
 3
 4
 5 
```

```@raw html
&nbsp;
```

#### Compute means, variances and error metrics through GBP iterations
Further, fields `rmse`, `mae` and `wrss` become vectors if we use the command `out = [iteration, evaluation]`, where each error value corresponds to the iteration according to the vector `iteration`.

```julia-repl
results, ~ = gbp(H, z, v; algorithm = vanilla, max = 5, out = [iteration, evaluation])

julia> results.gbp.rmse
5-element Vector{Float64}:
 1.8058973401950722
 2.904011755263201
 1.7463457881882882
 2.764690391451094
 1.6388328035960575

julia> results.gbp.mae
5-element Vector{Float64}:
 1.8046205188319902
 2.4086945100967365
 1.745111115124653
 2.293135757520099
 1.6376741235455516

julia> results.gbp.wrss
5-element Vector{Float64}:
 -3.6092410376639803
  4.817389020193473
 -3.490222230249306
  4.586271515040198
 -3.2753482470911033
```

```@raw html
&nbsp;
```

#### Compute means, variances and error metrics in the dynamic framework
The dynamic framework allows computing only means and variances, and error metrics if we use the command `out = evaluation`.
In the dynamic framework, means, variances and/or error metrics are evaluated before each new measurement update. Thus, fields `mean` and `variance` are given as matrices, where each column contains mean and variance values in the corresponding iteration according to the vector `iteration`. Note that according to the variable `d`, updates occur in the fifth and eighth iteration. 

```julia-repl
results, ~ = gbp(H, z, v, d; algorithm = vanillaDynamic)

julia> results.gbp.mean
2×3 Matrix{Float64}:
 -0.413415  -0.161671  -1.25493
  0.335673   0.944561   1.2827

julia> results.gbp.variance
2×3 Matrix{Float64}:
 45748.4  76467.6  4035.08
 45748.0  19535.0  4033.72

julia> results.gbp.iteration
3-element Vector{Int64}:
  4
  7
 30
``` 

```@raw html
&nbsp;
```

#### Compute means, variances and error metrics in the dynamic framework through iteration
Similar to before, using the keyword `out = iteration`, fields `mean` and `variance` are given as matrices, where each column contains mean and variance values in the corresponding iteration according to the vector `iteration`. In addition, using the keyword `out = [iteration evaluation]` error metrics can be also evaluated through iterations. 

```julia-repl
results, ~ = gbp(H, z, v, d; algorithm = vanillaDynamic, max = 5, out = iteration)

ulia> results.gbp.mean
2×5 Matrix{Float64}:
 0.885904  -0.108862   0.643748  -0.413415  0.0664998
 0.671831   0.0883903  0.822573   0.335673  0.819545

julia> results.gbp.variance
2×5 Matrix{Float64}:
 98361.6  48877.4  94041.0  45748.4  86286.4
 25127.9  48877.3  24024.2  45748.0  22043.2

julia> results.gbp.iteration
5-element Vector{Int64}:
 1
 2
 3
 4
 5
``` 
---

## The WLS Results
The struct variable `results.wls` contains results obtained using the WLS, which can be used to compare results obtained by the GBP algorithm. To describe the outputs, we will use the example given in [the GBP algorithm results] (@ref gbpresults) section. 

```@raw html
&nbsp;
```

#### Compute WLS solution and error metrics
Using the command `out = wls`, we obtained error metrics `rmse`, `mae` and `wrss` according to the WLS solution `mean`. Fields `rmseGBPWLS` and `maeGBPWLS` determine distances between the GBP estimate and WLS estimate after the GBP algorithm reaches the maximum number of iterations.

```julia-repl
results, ~ = gbp(H, z, v, d; algorithm = vanilla, out = wls)

julia> results.wls.mean
2-element Vector{Float64}:
 -6.457142857142001
  5.242857142856556

julia> results.wls.rmse
1-element Vector{Float64}:
 8.537275777709492e-14

julia> results.wls.mae
1-element Vector{Float64}:
 7.938094626069869e-14

julia> results.wls.wrss
1-element Vector{Float64}:
 -6.283862319378386e-14

julia> results.wls.rmseGBPWLS
1-element Vector{Float64}:
 0.17925300226918328

julia> results.wls.maeGBPWLS
1-element Vector{Float64}:
 1.2212117196693022
```

```@raw html
&nbsp;
```

#### Compute WLS solution and error metrics through GBP iterations
Using command `out = [iteration, wls]`, except variables `mean`, `rmse`, `mae` and `wrss`, we obtained fields `rmseGBPWLS` and `maeGBPWLS` for the each GBP iteration, where each element is related to the elements of the vector `results.gbp.iteration`

```julia-repl
results, ~ = gbp(H, z, v, d; algorithm = vanilla, max = 5, out = [iteration, wls])

julia> results.wls.rmseGBPWLS
5-element Vector{Float64}:
 1.9601144215518003
 0.8441540294119615
 1.895475083793785
 0.8036572271155497
 1.7787820268143777

julia> results.wls.maeGBPWLS
5-element Vector{Float64}:
 5.957036597700661
 5.751373926050387
 5.760587282671439
 5.475456017760832
 5.405943343568386

julia> results.gbp.iteration
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```

```@raw html
&nbsp;
```

#### Compute WLS solution and error metrics in the dynamic framework
In the dynamic framework, means and error metrics are evaluated before each new measurement update. Thus, field `mean` is given as matrix, where each column contains mean values in the corresponding iteration according to the vector `results.gbp.iteration`. Note that according to the variable `d`, updates occur in the fifth and eighth iteration. Also, error metrics are evaluated at the same iterations. 


```julia-repl
results, ~ = gbp(H, z, v, d; algorithm = vanillaDynamic, out = wls)

julia> results.wls.mean
2×3 Matrix{Float64}:
 -6.45714  -1.6  -1.27143
  5.24286   1.6   1.37857

julia> results.wls.rmse
3-element Vector{Float64}:
 8.537275777709492e-14
 1.3892435022089475e-14
 2.0123079445926824e-14

julia> results.wls.mae
3-element Vector{Float64}:
 7.938094626069869e-14
 1.3655743202889425e-14
 2.0039525594484076e-14

julia> results.wls.wrss
3-element Vector{Float64}:
 -6.283862319378386e-14
 -8.807769328692908e-15
  1.2163110025337826e-14

julia> results.wls.rmseGBPWLS
3-element Vector{Float64}:
 0.8036572271155497
 0.553587197227041
 0.05612434414868057

julia> results.wls.maeGBPWLS
3-element Vector{Float64}:
 5.475456017760832
 1.0468839767712772
 0.05618486344011808

julia> results.gbp.iteration
3-element Vector{Int64}:
  4
  7
 30
```

```@raw html
&nbsp;
```

#### Compute means, variances and error metrics in the dynamic framework through iteration
Using command `out = [iteration, wls]`, variables `mean`, `rmse`, `mae` and `wrss` are still evaluated before each new measurement update. Further, variables `rmseGBPWLS` and `maeGBPWLS` are obtained for each GBP iteration, where each element is related to the elements of the vector `results.gbp.iteration`.

```julia-repl
results, ~ = gbp(H, z, v, d; algorithm = vanillaDynamic, max = 5, out = [iteration, wls])

julia> results.wls.rmseGBPWLS
5-element Vector{Float64}:
 1.1014847586826113
 0.014475633457227476
 1.0368454209245965
 0.8036572271155497
 0.6265283007716769

julia> results.wls.maeGBPWLS
5-element Vector{Float64}:
 0.7788673422381023
 0.09862607394889142
 0.7331604281780024
 5.475456017760832
 1.2234773564418984

julia> results.gbp.iteration
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```
---

## Error Metrics 

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

#### The WLS and GBP error metrics
Fields `rmseGBPWLS` and `maeGBPWLS` determine distance beetwen the GBP estimate ``\hat{x}_{\text{gbp},i}`` and WLS estimate ``\hat{x}_{\text{wls},i}``, where root mean square error and mean absolute error are obtained using:
```math
  \begin{aligned}
    \text{rmse} = \sqrt {\cfrac{\sum_{i=1}^n \left[\hat{x}_{\text{wls},i} - \hat{x}_{{\text{gbp}},i}) \right]^2}{n}}; \quad
    \text{mae} = {\cfrac{\sum_{i=1}^n \left|\hat{x}_{\text{wls},i} - \hat{x}_{{\text{gbp}},i}) \right|}{n}},
  \end{aligned}
```  
where ``n`` is the number of state variabels.
