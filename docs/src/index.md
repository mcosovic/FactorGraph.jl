GaussBP
=============

The GaussBP solver provides the solution of the linear system of equations with/without Gaussian noise using the Gaussian belief propagation (GBP) algorithm applied over the factor graph in a static or dynamic framework.

The software package includes:
 - [Vanilla GBP algorithm] (@ref vanillaGBP);
 - [Computation-efficient GBP algorithm] (@ref efficientGBP);
 - [Computation-efficient GBP algorithm with Kahan–Babuška algorithm] (@ref kahanGBP);
 - [Dynamic GBP algorithm] (@ref dynamicGBP).
---


### Requirements 
GaussBP requires Julia 1.6 and higher. 

---

### Installation  
To install the GaussBP package, run the following command:
```julia-repl
pkg> add GaussBP
```

To load the package, use the command:
```julia-repl
julia> using GaussBP
```
---


###  Quick Start
```julia-repl
using GaussBP

results, system = gbp("data33_14.h5"; max = 20, out = terminal)
```

```julia-repl
using GaussBP

results, system = gbp("data33_14.h5"; algorithm = efficient, out = [evaluation, terminal])
```

```julia-repl
using GaussBP
using Plots

results, system = gbp("data33_14.xlsx"; variance = 1e60, out = [iteration, evaluation])
plot(results.gbp.iteration, results.gbp.rmse)
```

```julia-repl
using GaussBP

H = [1.5 0.0 2.0; 0.0 3.1 4.6; 2.6 8.1 0.4]
z = [0.8; 4.1; 2.2]
v = [1.0; 1.0; 1.0]     

results, settings = gbp(H, z, v; algorithm = kahan, out = [wls, terminal])
```

```julia-repl
using GaussBP

H = [1.5 0.0 2.0; 0.0 3.1 4.6; 2.6 8.1 0.4]
z = [0.8; 4.1; 2.2]
v = [1.0; 1.0; 1.0]  
d = [2 3 2.4 1.5; 15 1 0.85 0.9]

results, settings = gbp(H, z, v, d; algorithm = vanillaDynamic, out = terminal)
```
---

### More Information
- M. Cosovic and D. Vukobratovic, "Distributed Gauss-Newton Method for State Estimation Using Belief Propagation," in IEEE Transactions on  Power Systems, vol. 34, no. 1, pp. 648-658, Jan. 2019. [arxiv.org](https://arxiv.org/pdf/1702.05781.pdf)
- M. Cosovic, "Design and Analysis of Distributed State Estimation Algorithms Based on Belief Propagation and Applications in Smart Grids." arXiv preprint arXiv:1811.08355 (2018). [arxiv.org](https://arxiv.org/pdf/1811.08355.pdf)

