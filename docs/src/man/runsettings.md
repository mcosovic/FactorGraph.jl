# [Run Settings](@id runpf)

Input arguments of the function `gbp()` describe the GBP algorithm settings. The order of inputs and their appearance is arbitrary, with only DATA input required. Still, for methodological reasons, the syntax examples follow a certain order.

#### Syntax
```julia-repl
gbp(DATA)
gbp(DATA; ALGORITHM)
gbp(DATA; ALGORITHM, ITERATIONS)
gbp(DATA; ALGORITHM, ITERATIONS, VIRTUAL)
gbp(DATA; ALGORITHM, ITERATIONS, VIRTUAL, OUT)
```
```@raw html
&nbsp;
```
#### Description
```julia-repl
gbp(DATA) runs the GBP algorithm using DATA input 
gbp(DATA; ALGORITHM) selects the type of the GBP algorithm 
gbp(DATA; ALGORITHM, ITERATIONS) sets the iteration scheme
gbp(DATA; ALGORITHM, ITERATIONS, VIRTUAL) sets virtual factor nodes
gbp(DATA; ALGORITHM, ITERATIONS, VIRTUAL, OUT) controls output variable structure and display 
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
results, system = gbp("data33_14.csv"; variance = 1e60, out = ["iterate" "error"])
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
| **Example**     | `"data33_14.csv"`                                                           |
| **Description** |  loads the system data using csv-file from the package                      |
|                 |                                                                             |
| **Example**     | `"C:/name.h5"`                                                              |
| **Description** | loads the system data using h5-file from a custom path                      |
|                 |                                                                             |
| **Example**     | `"C:/name.csv"`                                                             |
| **Description** | loads the system data using csv-file from a custom path                     |
|                 |                                                                             |
| **Example**     | `Jacobian, means, variances`                                                |
| **Description** | loads the system data passing arguments directly                            |


---

## Keyword Arguments
The GBP function `gbp()` receives a group of arguments by keyword: ALGORITHM, ITERATIONS, VIRTUAL, OUT.

| ALGORITHM       | selects the type of the algorithm                                                         |
|:----------------|:------------------------------------------------------------------------------------------|
|                 |                                                                                           |
| **Command**     | `algorithm = "gbp"`                                                                       |
| **Description** |  runs the solver using native GBP algorithm, `default ALGORITHM setting`                  |
|                 |                                                                                           |
| **Command**     | `algorithm = "efficient"`                                                                 |
| **Description** |  runs the solver using computation-efficient GBP algorithm                                |
|                 |                                                                                           |
| **Command**     | `algorithm = "kahan"`                                                                     |
| **Description** |  runs the solver using computation-efficient GBP algorithm with compensated summation     |

```@raw html
&nbsp;
```

| ITERATIONS      | sets the iteration scheme                                                   |
|:----------------|:----------------------------------------------------------------------------|
|                 |                                                                             |
| **Command**     | `max = value`                                                               |
| **Description** |  specifies the maximum number of iterations, `default setting: max = 30`    |
|                 |                                                                             |
| **Command**     | `bump = value`                                                              |
| **Description** |  specifies the iteration where suspend the computation of variances (in a usual scenario, variances converge much faster than means) `default setting: bump = max` |
|                 |                                                                             |
| **Command**     | `damp = value`                                                              |
| **Description** |  specifies the iteration where applied randomized damping, `default setting: damp = max` |
|                 |                                                                             |
| **Command**     | `prob = value`                                                              |   
| **Description** |  a Bernoulli random variable with probability `prob = value` independently sampled for each mean value message from factor node to a variable node, applied for randomized damping iteration scheme with `value` between 0 and 1, `default setting: prob = 0.6`     |
|                 |                                                                             |
| **Command**     | `alpha = value`                                                             |
| **Description** |  the damped message is evaluated as a linear combination of the message from the previous and the current iteration, with weights `alpha = value` and `1 - alpha`, applied for randomized damping iteration scheme where `alpha` is between 0 and 1, `default setting: alpha = 0.4` |

!!! note "Randomized Damping"
    We provide an improved GBP algorithm that applies synchronous scheduling with randomized damping. The randomized damping parameter pairs lead to a trade-off between the number of non-converging simulations and the rate of convergence. In general, for the selection of `prob` and `alpha` for which only a small fraction of messages are combined with their values in a previous iteration, and that is a case for `prob` close to 0 or `alpha` close to 1, we observe a large number of non-converging simulations.


```@raw html
&nbsp;
```

| VIRTUAL         | sets virtual factor nodes                        |
|:----------------|:-------------------------------------------------|
|                 |                                                  |
| **Command**     | `mean = value`                                   |
| **Description** |  the mean value of virtual factor nodes          |
|                 |                                                  |
| **Command**     | `variance = value`                               |
| **Description** |  the variance value of the virtual factor nodes  |

!!! note "Virtual Factor Nodes"
    The virtual factor node is a singly-connected factor node used if the variable node is not directly observed. In a usual scenario, without prior knowledge, the variance of virtual factor nodes tend to infinity.

```@raw html
&nbsp;
```

| OUT             | controls output variable structure and display   |
|:----------------|:-------------------------------------------------|
|                 |                                                  |
| **Command**     | `out = "iterate"`                                |
| **Description** |  saves means and variances of the marginals through iterations |
|                 |                                                  |
| **Command**     | `out = "error"`                                  |
| **Description** |  computes error metrics of the GBP algorithm, note that if combined `"error"` with `iterate` (`out = ["iterate", "error"]`) computes error metrics through iterations                 |
|                 |                                                  |
| **Command**     | `out = "wls"`                                    |
| **Description** |  computes the solution using WLS method and error metrics of the GBP algorithm according to WLS solution, note that if combined `"wls"` with `iterate` (`out = ["iterate", "wls"]`) computes error metrics through iterations  |
|                 |                                                  |
| **Command**     | `out = "display"`                                |
| **Description** |  shows data display in the Julia REPL            |


!!! note "OUT"
    The keyword `out` accept any subset of commands `"iterate"`, `"error"`, `"wls"`, `"display"`. 