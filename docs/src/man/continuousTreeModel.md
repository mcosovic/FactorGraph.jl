# [Tree Graphical Model](@id graphicalTreeModelContinuous)

The FactorGraph supports the composite type `ContinuousTreeModel` based on the [forward–backward message passing schedule] (@ref synchronousSchedule), with three fields:
- `ContinuousTreeGraph`;
- `ContinuousInference`;
- `ContinuousSystem`.

The subtype `ContinuousTreeGraph` describes the tree factor graph obtained based on the input data. The GBP inference and marginal values are kept in the subtype `ContinuousInference`. The system of the linear equations being solved is preserved in the subtype `ContinuousSystem`. Note that the function `continuousTreeModel()` returns the main composite type `ContinuousTreeModel` with all subtypes.

---

#### Build graphical model

Input arguments of the function `continuousTreeModel()` describe the tree graphical model, while the function returns `ContinuousTreeModel` type.

Loads the system data passing arguments:
```julia-repl
gbp = continuousTreeModel(coefficient, observation, variances)
```

---

#### Virtual factor nodes

The function `continuousTreeModel()` receives arguments by keyword to set the mean and variance of the virtual factor nodes to initiate messages from leaf variable nodes if the corresponding variable node does not have a singly connected factor node.

```julia-repl
gbp = continuousTreeModel(DATA; mean = value, variance = value)
```
Default setting of the mean value is `mean = 0.0`, while the default variance is equal to `variance = 1e10`.

---

#### Root variable node

The function `continuousTreeModel()` receives argument by keyword to set the root variable node.
```julia-repl
gbp = continuousTreeModel(DATA; root = index)
```
Default setting of the root variable node is `root = 1`.

---

#### Tree factor graph
Function checks whether the factor graph has a tree structure.
```julia-repl
tree = isTree(gbp)
```
The tree structure of tha factor graph is marked as `tree = true`, the opposite is `tree = false`. The function accepts the composite type `ContinuousTreeModel`, as well as the type `ContinuousModel`.