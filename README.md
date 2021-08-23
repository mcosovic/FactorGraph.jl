# GaussBP

[![Documentation][documentation-badge]][documentation] ![Build][build-badge]


<a href="https://mcosovic.github.io/GaussBP.jl/stable/"><img align="left" width="320" src="/docs/src/assets/logo2.svg" /></a>

GaussBP is an open-source, easy-to-use simulation tool/solver for researchers and educators provided as a Julia package, with source code released under MIT License. The GaussBP package provides the set of different functions to perform inference over the factor graph in a static or dynamic framework using the linear Gaussian belief propagation (GBP) algorithm. The linear GBP model requires the set of linear equations and provides the minimum mean squared error (MMSE) estimate of the state variables.

We have tested and verified simulation tool using different scenarios to the best of our ability. As a user of this simulation tool, you can help us to improve future versions, we highly appreciate your feedback about any errors, inaccuracies, and bugs. For more information, please visit [documentation][documentation] site.

---

#### Requirement
GaussBP requires Julia 1.6 and higher.

---

#### Installation
To install the GaussBP package, run the following command:
```julia-repl
pkg> add GaussBP
```

To use GaussBP package, add the following code to your script, or alternatively run the same command in command prompt:
```julia-repl
using GaussBP
```
---


#### Quick Start: MMSE estimator
```julia-repl
using GaussBP

gbp = graphicalModel("data33_14.h5")        # initialize the graphical model using HDF5 input
for iteration = 1:200                       # the GBP inference
    messageFactorVariableVanilla(gbp)       # compute message using the native GBP
    messageVariableFactorVanilla(gbp)       # compute message using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```

```julia-repl
using GaussBP

gbp = graphicalModel("data33_14.xlsx")      # initialize the graphical model using XLSX input
for iteration = 1:200                       # the GBP inference
    messageFactorVariableEfficient(gbp)     # compute message using the efficient GBP
    messageVariableFactorEfficient(gbp)     # compute message using the efficient GBP
    marginal(gbp)                           # compute marginals in each iteration
end
exact = wls(gbp)                            # compute WLS estimate
displayData(gbp, exact)                     # show results
```

```julia-repl
using GaussBP

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = graphicalModel(H, z, v)               # initialize the graphical model via arguments
for iteration = 1:50                        # the GBP inference
    messageFactorVariableEfficient(gbp)     # compute message using the efficient GBP
    messageVariableFactorEfficient(gbp)     # compute message using the efficient GBP
end
marginal(gbp)                               # compute marginals
```

```julia-repl
using GaussBP

gbp = graphicalModel("data33_14.h5")        # initialize the graphical model
for iteration = 1:30                        # the GBP inference
    messageFactorVariableEfficient(gbp)     # compute message using the efficient GBP
    messageVariableFactorEfficient(gbp)     # compute message using the efficient GBP
end
for iteration = 31:200                      # continues the GBP inference
    meanFactorVariableVanilla(gbp)          # compute only means using the native GBP
    meanVariableFactorVanilla(gbp)          # compute only means using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```

```julia-repl
using GaussBP
using Plots

gbp = graphicalModel("data33_14.h5")        # initialize the graphical model
x6 = []                                     # save the state variable marginal
for iteration = 1:50                        # the GBP inference
    messageFactorVariableKahan(gbp)         # compute message using the GBP with Kahan-Babuska
    messageVariableFactorKahan(gbp)         # compute message using the GBP with Kahan-Babuska
    marginal(gbp)                           # compute marginals
    push!(x6, gbp.inference.mean[6])        # save state variable marginal
end
plot(collect(1:50), x6)                     # show plot
```

---

#### Quick Start: MMSE estimator in the dynamic framework
```julia-repl
using GaussBP

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = graphicalModel(H, z, v)               # initialize the graphical model
for iteration = 1:200                       # the GBP inference
    messageFactorVariableVanilla(gbp)       # compute message using the native GBP
    messageVariableFactorVanilla(gbp)       # compute message using the native GBP
end

dynamic = [1 0.85 1e-10; 3 2.4 1e-10]       # factor nodes change the mean and variance values
dynamicInference(gbp, dynamic)              # integrate changes in the running GBP
for iteration = 201:400                     # continues the GBP inference
    messageFactorVariableVanilla(gbp)       # compute message using the native GBP
    messageVariableFactorVanilla(gbp)       # compute message using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```
---

#### Quick Start: MMSE estimator in the dynamic ageing framework
```julia-repl
using GaussBP

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = graphicalModel(H, z, v)               # initialize the graphical mode
for iteration = 1:200                       # the GBP inference
    messageFactorVariableVanilla(gbp)       # compute message using the native GBP
    messageVariableFactorVanilla(gbp)       # compute message using the native GBP
end

dynamic = [1 0.85 1.0 1 0.05 0 50]          # factor node with the linear ageing
dynamicInference(gbp, dynamic)              # integrate changes in the running GBP
for iteration = 201:400                     # continues the GBP inference
    ageingInference(gbp, dynamic)           # integrate variance ageing
    messageFactorVariableVanilla(gbp)       # compute message using the native GBP
    messageVariableFactorVanilla(gbp)       # compute message using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```


[documentation-badge]: https://github.com/mcosovic/GaussBP.jl/workflows/Documentation/badge.svg
[build-badge]: https://github.com/mcosovic/GaussBP.jl/workflows/Build/badge.svg
[documentation]: https://mcosovic.github.io/GaussBP.jl/stable/