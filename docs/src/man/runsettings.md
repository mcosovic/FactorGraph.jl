# [Run Settings](@id runpf)

Input arguments of the function `gbp()` describe the GBP algorithm settings. The order of inputs and their appearance is arbitrary, with only DATA input required. Still, for methodological reasons, the syntax examples follow a certain order.

#### Syntax
```julia-repl
gbp(DATA)
gbp(DATA; ALGORITHM)
gbp(DATA; ALGORITHM, ITERATION)
gbp(DATA; ALGORITHM, ITERATION, VIRTUAL)
gbp(DATA; ALGORITHM, ITERATION, VIRTUAL, OUT)
```
```@raw html
&nbsp;
```
#### Description
```julia-repl
gbp(DATA) runs the GBP algorithm using DATA input 
gbp(DATA; ALGORITHM) selects the type of the GBP algorithm 
gbp(DATA; ALGORITHM, ITERATION) sets the iteration scheme
gbp(DATA; ALGORITHM, ITERATION, VIRTUAL) sets virtual factor nodes
gbp(DATA; ALGORITHM, ITERATION, VIRTUAL, OUT) controls output variable structure and display 
```
```@raw html
&nbsp;
```
#### Output
```julia-repl
results, system = gbp() returns results and input system data
```
```@raw html
&nbsp;
```
####  Examples
```julia-repl
results, system = gbp("data897_300.h5"; algorithm = "efficient", out = "display")
```
```julia-repl
results, system = gbp("data897_300.h5"; algorithm = "kahan", out = ["error", "display"])
```
```julia-repl
results, system = gbp("data33_14.xlsx"; variance = 1e60, out = ["iterate", "display"])
```
---


## Variable Arguments
The GBP function `gbp()` receives the variable argument DATA. 


| DATA            | input system data                                                           |
|:----------------|:----------------------------------------------------------------------------|
|                 |                                                                             |
| **Example**     | `"data33_14.h5"`                                                            |
| **Description** | loads the system data using h5-file from the package                        |
|                 |                                                                             |
| **Example**     | `"data33_14.xlsx"`                                                          |
| **Description** |  loads the system data using xlsx-file from the package                     |
|                 |                                                                             |
| **Example**     | `"C:/name.h5"`                                                              |
| **Description** | loads the system data using h5-file from a custom path                      |
|                 |                                                                             |
| **Example**     | `"C:/name.xlsx"`                                                            |
| **Description** | loads the system data using xlsx-file from a custom path                    |
|                 |                                                                             |
| **Example**     | `jacobian, observation, variances, dynamic`                                 |
| **Description** | loads the system data passing arguments directly                            |


---

## Keyword Arguments
The GBP function `gbp()` receives a group of arguments by keyword: ALGORITHM, ITERATION, INITIAL, OUT.

| ALGORITHM       | selects the type of the algorithm                                                                   |
|:----------------|:----------------------------------------------------------------------------------------------------|
|                 |                                                                                                     |
| **Command**     | `algorithm = "vanilla"`                                                                             |
| **Description** |  runs the solver using the native GBP algorithm, `default ALGORITHM setting`                        |
|                 |                                                                                                     |
| **Command**     | `algorithm = "efficient"`                                                                           |
| **Description** |  runs the solver using the computation-efficient GBP algorithm                                      |
|                 |                                                                                                     |
| **Command**     | `algorithm = "kahan"`                                                                               |
| **Description** |  runs the solver using the computation-efficient GBP algorithm with compensated summation           |
|                 |                                                                                                     |
| **Command**     | `algorithm = "vanillaDynamic"`                                                                      |
| **Description** |  runs the dynamic solver using the native GBP algorithm                                             |
|                 |                                                                                                     |
| **Command**     | `algorithm = "efficientDynamic"`                                                                    |
| **Description** |  runs the dynamic solver using the computation-efficient GBP algorithm                              |
|                 |                                                                                                     |
| **Command**     | `algorithm = "kahanDynamic"`                                                                        |
| **Description** |  runs the dynamic solver using the computation-efficient GBP algorithm with compensated summation   |
|                 |                                                                                                     |
| **Command**     | `algorithm = "vanillaAgeing"`                                                                       |
| **Description** |  runs the ageing solver using the native GBP algorithm                                              |
|                 |                                                                                                     |
| **Command**     | `algorithm = "efficientAgeing"`                                                                     |
| **Description** |  runs the ageing solver using the computation-efficient GBP algorithm                               |
|                 |                                                                                                     |
| **Command**     | `algorithm = "kahanAgeing"`                                                                         |
| **Description** |  runs the ageing solver using the computation-efficient GBP algorithm with compensated summation    |

```@raw html
&nbsp;
```

| ITERATION       | sets the iteration scheme                                                   |
|:----------------|:----------------------------------------------------------------------------|
|                 |                                                                             |
| **Command**     | `max = value`                                                               |
| **Description** |  specifies the maximum number of iterations, `default setting: max = 30`    |
|                 |                                                                             |
| **Command**     | `bump = value`                                                              |
| **Description** |  specifies the iteration when to suspend the computation of variances (in a usual scenario, variances converge much faster than means) `default setting: bump = max` |
|                 |                                                                             |
| **Command**     | `damp = value`                                                              |
| **Description** |  specifies the iteration where applied randomised damping, `default setting: damp = max` |
|                 |                                                                             |
| **Command**     | `prob = value`                                                              |   
| **Description** |  a Bernoulli random variable with probability `prob = value` independently sampled for each mean value message from a factor node to a variable node, applied for randomised damping iteration scheme with `value` between 0 and 1, `default setting: prob = 0.6`     |
|                 |                                                                             |
| **Command**     | `alpha = value`                                                             |
| **Description** |  the damped message is evaluated as a linear combination of the message from the previous and the current iteration, with weights `alpha = value` and `1 - alpha`, applied for randomised damping iteration scheme where `alpha` is between 0 and 1, `default setting: alpha = 0.4` |

!!! note "Randomised Damping"
    We provide an improved GBP algorithm that applies synchronous scheduling with randomised damping. The randomised damping parameter pairs lead to a trade-off between the number of non-converging simulations and the rate of convergence. In general, for the selection of `prob` and `alpha` for which only a small fraction of messages are combined with their values in a previous iteration, and that is a case for `prob` close to 0 or `alpha` close to 1, we observe a large number of non-converging simulations.


```@raw html
&nbsp;
```

| VIRTUAL         | sets virtual factor nodes                                                           |
|:----------------|:------------------------------------------------------------------------------------|
|                 |                                                                                     |
| **Command**     | `mean = value`                                                                      |
| **Description** |  the mean value of the virtual factor nodes, `default setting: mean = 0.0`          |
|                 |                                                                                     |
| **Command**     | `variance = value`                                                                  |
| **Description** |  the variance value of the virtual factor nodes, `default setting: variance = 1e10` |

!!! note "Virtual Factor Nodes"
    The virtual factor node is a singly-connected factor node used if the variable node is not directly observed. In a usual scenario, without prior knowledge, the variance of virtual factor nodes tend to infinity. Using virtual factor nodes, we also improved convergence performance.
    
```@raw html
&nbsp;
```

| OUT             | controls output variable structure and display                  |
|:----------------|:----------------------------------------------------------------|
|                 |                                                                 |
| **Command**     | `out = "iterate"`                                               |
| **Description** |  saves means and variances of the marginals through iterations  |
|                 |                                                                 |
| **Command**     | `out = "error"`                                                 |
| **Description** |  computes error metrics of the GBP algorithm, note that if combined `"error"` with `"iterate"` (`out = ["error", "iterate"]`) computes error metrics through iterations                                                                  |
|                 |                                                                 |
| **Command**     | `out = "wls"`                                                   |
| **Description** |  computes the solution using the WLS method and error metrics of the GBP algorithm according to the WLS solution, for the combination of `"wls"` with `"iterate"` (`out = ["iterate", "wls"]`) error metrics are computed through iterations |
|                 |                                                                  |
| **Command**     | `out = "display"`                                                |
| **Description** |  shows data display in the Julia REPL                            |


!!! note "OUT"
    The keyword `out` accepts any subset of commands `"iterate"`, `"error"`, `"wls"`, `"display"`. 