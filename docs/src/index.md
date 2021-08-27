GaussBP
=============

The GaussBP package provides the set of different functions to perform inference over the factor graph in a static or dynamic framework using the linear Gaussian belief propagation (GBP) algorithm. The linear GBP model requires the set of linear equations and provides the minimum mean squared error (MMSE) estimate of the state variables.

The software package includes:
 - [Vanilla GBP algorithm] (@ref vanillaGBP);
 - [Computation-efficient GBP algorithm] (@ref efficientGBP);
 - [Computation-efficient GBP algorithm with Kahan–Babuška algorithm] (@ref kahanGBP);
 - [Dynamic GBP algorithm] (@ref dynamicGBP);
 - [Ageing GBP algorithm] (@ref ageingGBP).
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

#### Quick start: MMSE estimator
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

#### Quick start: MMSE estimator in the dynamic framework
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
dynamicInference!(gbp, dynamic)             # integrate changes in the running GBP
for iteration = 201:400                     # continues the GBP inference
    messageFactorVariableVanilla(gbp)       # compute message using the native GBP
    messageVariableFactorVanilla(gbp)       # compute message using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```
---

#### Quick start: MMSE estimator in the dynamic ageing framework
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
dynamicInference!(gbp, dynamic)             # integrate changes in the running GBP
for iteration = 201:400                     # continues the GBP inference
    ageingInference!(gbp, dynamic)          # integrate variance ageing
    messageFactorVariableVanilla(gbp)       # compute message using the native GBP
    messageVariableFactorVanilla(gbp)       # compute message using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```

---

#### More information
- M. Cosovic and D. Vukobratovic, "Distributed Gauss-Newton Method for State Estimation Using Belief Propagation," in IEEE Transactions on  Power Systems, vol. 34, no. 1, pp. 648-658, Jan. 2019. [arxiv.org](https://arxiv.org/pdf/1702.05781.pdf)
- M. Cosovic, "Design and Analysis of Distributed State Estimation Algorithms Based on Belief Propagation and Applications in Smart Grids." arXiv preprint arXiv:1811.08355 (2018). [arxiv.org](https://arxiv.org/pdf/1811.08355.pdf)

