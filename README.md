# GaussBP

[![Documentation][documentation-badge]][documentation] ![Build][build-badge]


<a href="https://mcosovic.github.io/GaussBP.jl/stable/"><img align="left" width="320" src="/docs/src/assets/logo2.svg" /></a>

GaussBP is an open-source, easy-to-use simulation tool/solver for researchers and educators provided as a Julia package, with source code released under MIT License. The solver provides the solution of the linear system of equations with/without Gaussian noise using the Gaussian belief propagation (GBP) algorithm applied over the factor graph in a static or dynamic framework.

We have tested and verified simulation tool using different scenarios to the best of our ability. As a user of this simulation tool, you can help us to improve future versions, we highly appreciate your feedback about any errors, inaccuracies, and bugs. For more information, please visit [documentation][documentation] site.

---

### Requirement
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


### Quick Start: The GBP Algorithm
```julia-repl
using GaussBP

results, system = gbp("data33_14.h5"; max = 20, out = "display")
```
```julia-repl
using GaussBP

results, system = gbp("data33_14.h5"; algorithm = "efficient", out = ["error", "display"])
```

```julia-repl
using GaussBP
using Plots

results, system = gbp("data33_14.xlsx"; variance = 1e60, out = ["iterate", "error", "display"])
plot(results.gbp.iteration, results.gbp.rmse)
```
```julia-repl
using GaussBP

jacobian = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6; 2.6 8.1 0.4]
observation = [0.5; 0.8; 4.1; 2.2]
variance = [0.1; 1.0; 1.0; 1.0]    

results, ~ = gbp(jacobian, observation, variance; algorithm = "kahan", out = ["wls", "display"])
```
---

### Quick Start: The GBP Algorithm with Dynamic Updates 
```julia-repl
using GaussBP

jacobian = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6; 2.6 8.1 0.4]
observation = [0.5; 0.8; 4.1; 2.2]
variance = [0.1; 1.0; 1.0; 1.0]  
dynamic = [2 3 2.4 1.5; 15 1 0.85 0.9]

results, ~ = gbp(jacobian, observation, variance, dynamic; algorithm = "vanillaDynamic", out = "display")
```
---

### Quick Start: The GBP Algorithm with Ageing Information
```julia-repl
using GaussBP

jacobian = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6; 2.6 8.1 0.4]
observation = [0.5; 0.8; 4.1; 2.2]
variance = [0.1; 1.0; 1.0; 1.0]  
dynamic = [2 3 2.45 1.5 1 4 0.5 0 50; 8 1 0.85 0.9 2 12 0.3 0.1 30]

results, ~ = gbp(jacobian, observation, variance, dynamic; algorithm = "vanillaAgeing", out = "display")
```

[documentation-badge]: https://github.com/mcosovic/GaussBP.jl/workflows/Documentation/badge.svg
[build-badge]: https://github.com/mcosovic/GaussBP.jl/workflows/Build/badge.svg
[documentation]: https://mcosovic.github.io/GaussBP.jl/stable/