# [Inference](@id inferenceTreeDiscrete)

To exchange information over the tree factor graph, the FactorGraph provides forward–backward BP algorithm. We advise the reader to read the Section [forward–backward message passing schedule] (@ref treeSchedule) which provides a detailed description of the inference algorithm.

Each of the inference functions accepts only the composite type `DiscreteTreeModel`, i.e., an output variable of the function `bp = discreteTreeModel()`.

---

#### Forward inference
The set of functions that can be used to preform forward message inference:
```julia-repl
forwardVariableFactor(bp)
forwardFactorVariable(bp)
```
---

#### Backward inference
The set of functions that can be used to preform backward message inference:
```julia-repl
backwardVariableFactor(bp)
backwardFactorVariable(bp)
```
---

#### Marginal inference
To compute normalized marginals the FactorGraph provides the function:
```julia-repl
marginal(bp)
```
To compute unnormalized marginals the FactorGraph provides the function:
```julia-repl
marginalUnnormalized(bp)
```
Same as before, functions accept the composite type `DiscreteTreeModel`.