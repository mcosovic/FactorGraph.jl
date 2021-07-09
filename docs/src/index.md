GaussBP
=============

The GaussBP solver provides the solution of the linear system of equations with/without Gaussian noise using the Gaussian belief propagation (GBP) algorithm applied over the factor graph.

The software package includes:
 - [Native GBP algorithm] (@ref nativeGBP);
 - [Computation-efficient GBP algorithm] (@ref efficientGBP);
 - [Computation-efficient GBP algorithm with Kahan–Babuška algorithm] (@ref kahanGBP).
---


### Installation
GaussBP requires Julia 1.6 and higher. To install the GaussBP package, run the following command:
```julia-repl
pkg> add https://github.com/mcosovic/GaussBP.jl
```

To load the package, use the command:
```julia-repl
julia> using GaussBP
```
---


###  Quick Start
```julia-repl
using GaussBP

results, system = gbp("data33_14.h5"; max = 50, out = "display")
```
```julia-repl
using GaussBP

results, system = gbp("data33_14.h5"; algorithm = "kahan", out = ["error", "display"])
```

```julia-repl
using GaussBP
using Plots

results, system = gbp("data33_14.csv"; variance = 1e60, out = ["iterate" "error", "display"])
plot(results.iterations, results.rmse)
```
```julia-repl
using GaussBP

H = [1.5 0.0 2.0; 0.0 3.1 4.6; 2.6 8.1 0.5]
z = [0.8; 4.1; 2.2]
v = [1.0; 1.0; 1.0]     

results, settings = gbp(H, z, v; out = ["error", "wls", "display"], algorithm = "kahan")
```
---

### More Information
- M. Cosovic and D. Vukobratovic, "Distributed Gauss-Newton Method for State Estimation Using Belief Propagation," in IEEE Transactions on  Power Systems, vol. 34, no. 1, pp. 648-658, Jan. 2019. [arxiv.org](https://arxiv.org/pdf/1702.05781.pdf)
- M. Cosovic, "Design and Analysis of Distributed State Estimation Algorithms Based on Belief Propagation and Applications in Smart Grids." arXiv preprint arXiv:1811.08355 (2018). [arxiv.org](https://arxiv.org/pdf/1811.08355.pdf)

