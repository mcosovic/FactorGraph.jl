# [Tree Graphical Model](@id graphicalTreeModelDiscrete)

The FactorGraph supports the composite type `DiscreteTreeModel` based on the [forward–backward message passing schedule] (@ref treeSchedule), with three fields:
- `DiscreteTreeGraph`;
- `DiscreteInference`;
- `DiscreteSystem`.

The subtype `DiscreteTreeGraph` describes the tree factor graph obtained based on the input data. The BP inference and marginal values are kept in the subtype `DiscreteInference`. The input system being solved is preserved in the subtype `DiscreteSystem`. Note that the function `discreteTreeModel()` returns the main composite type `DiscreteTreeModel` with all subtypes.

---

#### Build graphical model

Input arguments of the function `discreteTreeModel()` describe the tree graphical model, while the function returns `DiscreteTreeModel` type.

Loads the system data passing arguments:
```julia-repl
bp = discreteTreeModel(probability, table)
```

---

#### Virtual factor nodes

The function `discreteTreeModel()` receives arguments by keyword to set the message of the virtual factor nodes to initiate messages from leaf variable nodes.  This value is applied to all variable nodes and to all their possible states if the corresponding variable node does not have a singly connected factor node.

```julia-repl
bp = discreteTreeModel(DATA; message = value)
```
Default setting of the mean value is `message = 1.0`.

---

#### Root variable node

The function `discreteTreeModel()` receives argument by keyword to set the root variable node.
```julia-repl
bp = discreteTreeModel(DATA; root = index)
```
Default setting of the root variable node is `root = 1`.

---

#### Tree factor graph
Function checks whether the factor graph has a tree structure.
```julia-repl
tree = isTree(gbp)
```
The tree structure of tha factor graph is marked as `tree = true`, the opposite is `tree = false`. The function accepts the composite type `DiscreteTreeModel`.