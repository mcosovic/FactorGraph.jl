GaussBP
=============

The GaussBP package provides the set of different functions to perform inference over the factor graph in a static or dynamic framework using the linear Gaussian belief propagation (GBP) algorithm. The linear GBP model requires the set of linear equations and provides the minimum mean squared error (MMSE) estimate of the state variables.

The software package includes algorithms based on the [synchronous message passing schedule] (@ref synchronous):
 - [vanilla GBP algorithm] (@ref vanillaGBP);
 - [computation-efficient GBP algorithm] (@ref efficientGBP);
 - [computation-efficient GBP algorithm with Kahan–Babuška algorithm] (@ref kahanGBP);
 - [dynamic GBP algorithm] (@ref dynamicGBP);
 - [ageing GBP algorithm] (@ref ageingGBP).

The software package also includes a message passing algorithm that allows exact inference in tree factor graph:
- [forward–backward algorithm] (@ref treeGBP).

---

#### Requirement
GaussBP requires Julia 1.6 and higher.

---

#### Installation
To install the GaussBP package, run the following command:
```julia-repl
pkg> add GaussBP
```

To use GaussBP package, add the following code to your script, or alternatively run the same command in Julia REPL:
```julia-repl
using GaussBP
```
---

#### Quick start
Following examples are intended for a quick introduction to GaussBP package.

- Synchronous message passing schedule using the native GBP algorithm.
```julia-repl
using GaussBP

gbp = graphicalModel("data33_14.h5")        # initialize the graphical model using HDF5 input
for iteration = 1:200                       # the GBP inference
    messageFactorVariableVanilla(gbp)       # compute messages using the native GBP
    messageVariableFactorVanilla(gbp)       # compute messages using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```

- Synchronous message passing schedule using the efficient GBP algorithm.
```julia-repl
using GaussBP

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = graphicalModel(H, z, v)               # initialize the graphical model via arguments
for iteration = 1:50                        # the GBP inference
    messageFactorVariableEfficient(gbp)     # compute messages using the efficient GBP
    messageVariableFactorEfficient(gbp)     # compute messages using the efficient GBP
end
marginal(gbp)                               # compute marginals
```

- Synchronous message passing schedule using the GBP and Kahan-Babuska algorithm with the plotting of the marginal mean through iteration.
```julia-repl
using GaussBP
using Plots

gbp = graphicalModel("data33_14.h5")        # initialize the graphical model
x6 = []                                     # save the state variable marginal
for iteration = 1:50                        # the GBP inference
    messageFactorVariableKahan(gbp)         # compute messages using the GBP with Kahan-Babuska
    messageVariableFactorKahan(gbp)         # compute messages using the GBP with Kahan-Babuska
    marginal(gbp)                           # compute marginals
    push!(x6, gbp.inference.mean[6])        # save state variable marginal
end
plot(collect(1:50), x6)                     # show plot
```

- Synchronous message passing schedule using the native GBP algorithm in the dynamic framework.
```julia-repl
using GaussBP

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = graphicalModel(H, z, v)               # initialize the graphical model
for iteration = 1:200                       # the GBP inference
    messageFactorVariableVanilla(gbp)       # compute messages using the native GBP
    messageVariableFactorVanilla(gbp)       # compute messages using the native GBP
end

dynamicFactor!(gbp;                         # integrate changes in the running GBP
    factor = 1,
    mean = 0.85,
    variance = 1e-10)
for iteration = 201:400                     # continues the GBP inference
    messageFactorVariableVanilla(gbp)       # compute messages using the native GBP
    messageVariableFactorVanilla(gbp)       # compute messages using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```

- Synchronous message passing schedule using the native GBP algorithm in the dynamic ageing framework.
```julia-repl
using GaussBP

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = graphicalModel(H, z, v)               # initialize the graphical model
for iteration = 1:200                       # the GBP inference
    messageFactorVariableVanilla(gbp)       # compute messages using the native GBP
    messageVariableFactorVanilla(gbp)       # compute messages using the native GBP
end

for iteration = 1:400                       # continues the GBP inference
    ageingVariance!(gbp;                    # integrate changes in the running GBP
        factor = 4,
        initial = 1,
        limit = 50,
        model = 1,
        a = 0.05,
        tau = iteration)
    messageFactorVariableVanilla(gbp)       # compute messages using the native GBP
    messageVariableFactorVanilla(gbp)       # compute messages using the native GBP
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results
```

 - Forward–backward algorithm over the tree factor graph.
```julia-repl
using GaussBP

H = [1 0 0 0 0; 6 8 2 0 0; 0 5 0 0 0;       # jacobian matrix
     0 0 2 0 0; 0 0 3 8 2]
z = [1; 2; 3; 4; 5]                         # observation vector
v = [3; 4; 2; 5; 1]                         # variance vector

gbp = graphicalModelTree(H, z, v; root = 1) # initialize the tree graphical model
while gbp.graph.forward                     # inference from leaves to the root
     forwardVariableFactor(gbp)             # compute forward messages
     forwardFactorVariable(gbp)             # compute forward messages
end
while gbp.graph.backward                    # inference from the root to leaves
     backwardVariableFactor(gbp)            # compute backward messages
     backwardFactorVariable(gbp)            # compute backward messages
end
marginal(gbp)                               # compute marginals
displayData(gbp)                            # show results

```

---

#### More information
- M. Cosovic and D. Vukobratovic, "Distributed Gauss-Newton Method for State Estimation Using Belief Propagation," in IEEE Transactions on  Power Systems, vol. 34, no. 1, pp. 648-658, Jan. 2019. [arxiv.org](https://arxiv.org/pdf/1702.05781.pdf)
- M. Cosovic, "Design and Analysis of Distributed State Estimation Algorithms Based on Belief Propagation and Applications in Smart Grids." arXiv preprint arXiv:1811.08355 (2018). [arxiv.org](https://arxiv.org/pdf/1811.08355.pdf)

