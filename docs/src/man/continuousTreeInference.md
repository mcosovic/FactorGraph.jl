# [Inference](@id inferenceTreeContinuous)

To exchange information over the tree factor graph, the FactorGraph provides forward–backward algorithm. We advise the reader to read the Section [forward–backward message passing] (@ref treeSchedule) which provides a detailed description of the inference algorithm.

Each of the inference functions accepts only the composite type `ContinuousTreeModel`, i.e., an output variable of the function `gbp = continuousTreeModel()`.

---

#### Forward inference
The set of functions that can be used to preform forward message inference:
```julia-repl
forwardVariableFactor(gbp)
forwardFactorVariable(gbp)
```
---

#### Backward inference
The set of functions that can be used to preform backward message inference:
```julia-repl
backwardVariableFactor(gbp)
backwardFactorVariable(gbp)
```
---

#### Marginal inference
To compute marginals the FactorGraph provides the function:
```julia-repl
marginal(gbp)
```
Same as before, the function accepts the composite type `ContinuousTreeModel`.