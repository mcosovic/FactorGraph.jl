# [Inference](@id vanilla)

To exchange information over the tree factor graph, the GaussBP provides forward–backward algorithm. We advise the reader to read the [forward–backward algorithm] (@ref treeGBP) which provides a detailed description of the inference algorithm.

Each of the inference functions accepts only the composite type `GraphicalModelTree`, i.e., an output variable of the function `gbp = graphicalModelTree()`.

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
To compute marginals the GaussBP provides the function:
```julia-repl
marginalTree(gbp)
```
Same as before, the function accepts only the composite type `GraphicalModelTree`.