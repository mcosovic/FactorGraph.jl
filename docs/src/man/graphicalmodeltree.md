# [Graphical Model Tree](@id graphicalModelTree)

The FactorGraph supports the composite type `GraphicalModelTree` related with the [forward–backward message passing] (@ref treeGBP), with three fields:
- `FactorGraphTree`;
- `Inference`;
- `SystemModel`.

The subtype `FactorGraphTree` describes the tree factor graph obtained based on the input data. The GBP inference and marginal values are kept in the subtype `Inference`. The system of the linear equations being solved is preserved in the subtype `SystemModel`. Note that the function `graphicalModelTree()` returns the main composite type `GraphicalModelTree` with all subtypes.

---

#### Build graphical model

Input arguments DATA of the function `graphicalModelTree()` describe the tree graphical model, while the function returns `GraphicalModelTree` type.

Loads the system data using h5-file from the package:
```julia-repl
gbp = graphicalModelTree("data12_11.h5")
```

Loads the system data using xlsx-file from the package:
```julia-repl
gbp = graphicalModelTree("data13_11.xlsx")
```

Loads the system data from a custom path:
```julia-repl
gbp = graphicalModelTree("C:/name.h5")
```

Loads the system data passing arguments directly:
```julia-repl
gbp = graphicalModelTree(jacobian, observation, variances)
```

---

#### Virtual factor nodes

The GBP function `graphicalModelTree()` receives arguments by keyword to set the mean and variance of the virtual factor nodes to initiate messages from leaf variable nodes if the corresponding variable node does not have a singly connected factor node.

```julia-repl
gbp = graphicalModelTree(DATA; mean = value, variance = value)
```
Default setting of the mean value is `mean = 0.0`, while the default variance is equal to `variance = 1e10`.

---

#### Root variable node

The GBP function `graphicalModelTree()` receives argument by keyword to set the root variable node.
```julia-repl
gbp = graphicalModelTree(DATA; root = index)
```
Default setting of the root variable node is `root = 1`.

---

#### Tree factor graph
Function checks whether the factor graph has a tree structure.
```julia-repl
tree = isTree(gbp)
```
The tree structure of tha factor graph is marked as `tree = true`, the opposite is `tree = false`. The function accepts the composite type `GraphicalModelTree`, as well as the type `GraphicalModel`.