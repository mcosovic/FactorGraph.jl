# GaussBP

[![Documentation][documentation-badge]][documentation] ![Build][build-badge]


<a href="https://mcosovic.github.io/GaussBP.jl/stable/"><img align="right" width="180" src="/docs/src/assets/logo2.png" /></a>

GaussBP is an open-source, easy-to-use simulation tool/solver for researchers and educators provided as a Julia package, with source code released under MIT License. The solver provides the solution of the linear system of equations with/without Gaussian noise using belief propagation algorithm applied over the factor graph.

We have tested and verified simulation tool using different scenarios to the best of our ability. As a user of this simulation tool, you can help us to improve future versions, we highly appreciate your feedback about any errors, inaccuracies, and bugs. For more information, please visit [documentation][documentation] site.

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

results, system = gbp("data33_14.csv"; variance = 1e60, out = ["iterate", "error", "display"])
plot(results.iterations, results.rmse)
```
```julia-repl
using GaussBP

H = [1.5 0.0 2.0; 0.0 3.1 4.6; 2.6 8.1 0.5]
z = [0.8; 4.1; 2.2]
v = [1.0; 1.0; 1.0]     

results, settings = gbp(H, z, v; out = ["error", "wls", "display"], algorithm = "kahan")
```

[documentation-badge]: https://github.com/mcosovic/GaussBP.jl/workflows/Documentation/badge.svg
[build-badge]: https://github.com/mcosovic/GaussBP.jl/workflows/Build/badge.svg
[documentation]: https://mcosovic.github.io/GaussBP.jl/stable/