FactorGraph
=============

The FactorGraph package provides the set of different functions to perform inference over the factor graph with continuous or discrete random variables using the belief propagation (BP) algorithm, also known as the sum-product algorithm.

----

#### Continuous framework
In the case of continuous random variables described by Gaussian distributions, we are using the linear Gaussian belief propagation (GBP) algorithm to solve the inference problem. The linear GBP model requires the set of linear equations and provides the minimum mean squared error (MMSE) estimate of the state variables. To perform inference the FactorGraph package uses several algorithms based on the [synchronous message passing schedule] (@ref synchronousSchedule):
 - [vanilla GBP algorithm] (@ref vanillaGBP);
 - [broadcast GBP algorithm] (@ref broadcastGBP);
 - [broadcast GBP with Kahan–Babuška algorithm] (@ref kahanGBP).
Within these algorithms, the packege provides several routines to allow dynamic GBP framework:
 - [dynamic GBP algorithm] (@ref dynamicGBP);
 - [ageing GBP algorithm] (@ref ageingGBP).
Finally, the package also includes a message passing algorithm that allows inference in the tree factor graph:
- [forward–backward GBP algorithm] (@ref treeSchedule).

---

#### Discrete framework
In the case of discrete random variables the package currently provides only the BP algorithm that allows exact inference in the tree factor graph:
- [forward–backward BP algorithm] (@ref treeSchedule).

---

#### Requirement
FactorGraph requires Julia 1.6 and higher.

---

#### Installation
To install the FactorGraph package, run the following command:
```julia-repl
pkg> add FactorGraph
```

To use FactorGraph package, add the following code to your script, or alternatively run the same command in Julia REPL:
```julia-repl
using FactorGraph
```
---

#### Quick start whitin continuous framework
The following examples are intended for a quick introduction to FactorGraph package within the continuous framework.

- The broadcast GBP algorithm:
```julia-repl
using FactorGraph

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = continuousModel(H, z, v)              # initialize the graphical model
for iteration = 1:50                        # the GBP inference
    messageFactorVariableBroadcast(gbp)     # compute messages using the broadcast GBP
    messageVariableFactorBroadcast(gbp)     # compute messages using the broadcast GBP
end
marginal(gbp)                               # compute marginals
```

- The vanilla GBP algorithm in the dynamic framework:
```julia-repl
using FactorGraph

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = continuousModel(H, z, v)              # initialize the graphical model
for iteration = 1:200                       # the GBP inference
    messageFactorVariable(gbp)              # compute messages using the vanilla GBP
    messageVariableFactor(gbp)              # compute messages using the vanilla GBP
end

dynamicFactor!(gbp;                         # integrate changes in the running GBP
    factor = 1,
    mean = 0.85,
    variance = 1e-10)
for iteration = 201:400                     # continues the GBP inference
    messageFactorVariable(gbp)              # compute messages using the vanilla GBP
    messageVariableFactor(gbp)              # compute messages using the vanilla GBP
end
marginal(gbp)                               # compute marginals
```

- The vanilla GBP algorithm in the dynamic ageing framework:
```julia-repl
using FactorGraph

H = [1.0 0.0 0.0; 1.5 0.0 2.0; 0.0 3.1 4.6] # jacobian matrix
z = [0.5; 0.8; 4.1]                         # observation vector
v = [0.1; 1.0; 1.0]                         # variance vector

gbp = continuousModel(H, z, v)              # initialize the graphical model
for iteration = 1:200                       # the GBP inference
    messageFactorVariable(gbp)              # compute messages using the vanilla GBP
    messageVariableFactor(gbp)              # compute messages using the vanilla GBP
end

for iteration = 1:400                       # continues the GBP inference
    ageingVariance!(gbp;                    # integrate changes in the running GBP
        factor = 3,
        initial = 1,
        limit = 50,
        model = 1,
        a = 0.05,
        tau = iteration)
    messageFactorVariable(gbp)              # compute messages using the vanilla GBP
    messageVariableFactor(gbp)              # compute messages using the vanilla GBP
end
marginal(gbp)                               # compute marginals
```

 - The forward–backward GBP algorithm over the tree factor graph:
```julia-repl
using FactorGraph

H = [1 0 0 0 0; 6 8 2 0 0; 0 5 0 0 0;       # jacobian matrix
     0 0 2 0 0; 0 0 3 8 2]
z = [1; 2; 3; 4; 5]                         # observation vector
v = [3; 4; 2; 5; 1]                         # variance vector

gbp = continuousTreeModel(H, z, v)          # initialize the tree graphical model
while gbp.graph.forward                     # inference from leaves to the root
     forwardVariableFactor(gbp)             # compute forward messages
     forwardFactorVariable(gbp)             # compute forward messages
end
while gbp.graph.backward                    # inference from the root to leaves
     backwardVariableFactor(gbp)            # compute backward messages
     backwardFactorVariable(gbp)            # compute backward messages
end
marginal(gbp)                               # compute marginals
```

---

#### Quick start whitin discrete framework
Following example is intended for a quick introduction to FactorGraph package within the discrete framework.

 - The forward–backward BP algorithm over the tree factor graph:
```julia-repl
using FactorGraph

probability1 = [1]
table1 = [0.2; 0.3; 0.4; 0.1]

probability2 = [1; 2; 3]
table2 = zeros(4, 3, 1)
table2[1, 1, 1] = 0.2; table2[2, 1, 1] = 0.5; table2[3, 1, 1] = 0.3; table2[4, 1, 1] = 0.0
table2[1, 2, 1] = 0.1; table2[2, 2, 1] = 0.1; table2[3, 2, 1] = 0.7; table2[4, 2, 1] = 0.1
table2[1, 3, 1] = 0.5; table2[2, 3, 1] = 0.2; table2[3, 3, 1] = 0.1; table2[4, 3, 1] = 0.1

probability = [probability1, probability2]
table = [table1, table2]

bp = discreteTreeModel(probability, table)  # initialize the tree graphical model
while bp.graph.forward                      # inference from leaves to the root
    forwardVariableFactor(bp)               # compute forward messages
    forwardFactorVariable(bp)               # compute forward messages
end
while bp.graph.backward                     # inference from the root to leaves
    backwardVariableFactor(bp)              # compute backward messages
    backwardFactorVariable(bp)              # compute backward messages
end
marginal(bp)                                # compute normalized marginals
```

---

#### More information
- M. Cosovic and D. Vukobratovic, "Distributed Gauss-Newton Method for State Estimation Using Belief Propagation," in IEEE Transactions on  Power Systems, vol. 34, no. 1, pp. 648-658, Jan. 2019. [arxiv.org](https://arxiv.org/pdf/1702.05781.pdf)
- M. Cosovic, "Design and Analysis of Distributed State Estimation Algorithms Based on Belief Propagation and Applications in Smart Grids." arXiv preprint arXiv:1811.08355 (2018). [arxiv.org](https://arxiv.org/pdf/1811.08355.pdf)