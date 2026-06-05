# [Discrete Factor Graph](@id discrete-factor-graph)

```@meta
CurrentModule = FactorGraph
```

A discrete factor graph represents a model with finite-state variables and
nonnegative factor tables. It is a bipartite graph in which variable nodes store
finite state sets, factor nodes store local potentials, and each edge connects a
factor to a variable that appears in that factor table.

The full model factorizes as a product of local potentials:

```math
p(\mathbf{x}) \propto
\prod_i \psi_i(\mathbf{x}_{\mathcal{V}_i}).
```

Each local potential is a nonnegative function:

```math
\psi_i :
\mathcal{X}_{\mathcal{V}_i}
\rightarrow
\mathbb{R}_{\ge 0}.
```

Here, ``\mathcal{V}_i`` is the ordered collection of variables connected to
factor ``i``, and ``\mathcal{X}_{\mathcal{V}_i}`` is the Cartesian product of
their state sets. In the implementation, ``\psi_i`` is stored as a factor table
whose axes follow the variable order used in the factor constructor. Factor
tables are potentials and do not need to be normalized.

---


## Variable Nodes

Variable nodes are created with [`DiscreteVariable`](@ref):

```@example discrete_factor_graph
using FactorGraph # hide

x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on], probability = [0.7, 0.3])
x2 = DiscreteVariable(:x2, 3; label = "x2", states = ["low", "mid", "high"])
x3 = DiscreteVariable(:x3, 2; label = "x3")

nothing # hide
```

---

##### Identifiers and Labels

The identifier, such as `:x1`, is the compact programmatic reference used by
factors. It can be either a `Symbol` or an `Integer`.

The `label` is a `String` used for printed output and lookup functions. Factor
definitions may also refer to variables by label.

---

##### States

The second argument defines the number of states. State references can be named
with the `states` keyword and can be integers, symbols, or strings. The state
order determines the axis order in connected factor tables.

In this example, factors connected to `:x1` use the order `[:off, :on]`, factors
connected to `:x2` use the order `["low", "mid", "high"]`, and factors connected
to `:x3` use the default order `[1, 2]`.

---

##### Initial Belief

A variable may optionally define its initial belief using the `probability`
keyword. This vector must have the same length as the variable cardinality, must
be finite and nonnegative, and must contain at least one positive entry.

If `probability` is omitted, [`sumproduct`](@ref) starts from a uniform initial
message over the variable states, and [`minsum`](@ref) uses the corresponding
uniform initial cost.

This initial belief is used only to initialize [`sumproduct`](@ref) messages and
[`minsum`](@ref) costs. It is not a substitute for a unary factor when the prior
should be part of the probabilistic model.

---

## Factor Nodes

Factor nodes are created with [`DiscreteFactor`](@ref):

```@example discrete_factor_graph
ψ1 = [0.8, 0.2]
f1 = DiscreteFactor(x1, ψ1; label = "f1", initialize = true)

ψ2 = [1.0 0.2 0.4; 0.1 1.1 0.3]
f2 = DiscreteFactor(x1, x2, ψ2; label = "f2")

ψ3 = cat([0.9 0.5 0.2; 0.1 0.5 0.8], [0.7 0.4 0.3; 0.3 0.6 0.7]; dims = 3)
f3 = DiscreteFactor(x1, x2, x3, ψ3; label = "f3")

nothing # hide
```

---

##### Variable References

The first arguments define the variable nodes connected to the factor node. Their
order is important and must match the dimensions of the factor table. A variable
can be referenced by identifier, such as `:x1` or `1`, by label, such as `"x1"`,
or by passing the [`DiscreteVariable`](@ref) node itself.

---

##### Factor Table

The factor table stores a nonnegative potential over the connected variables
``\psi_i(\mathbf{x}_{\mathcal{V}_i})``. The table dimensions follow the order of
the variables passed to [`DiscreteFactor`](@ref). For a factor connected to
variables with cardinalities `2` and `3`, the table must have size `(2, 3)`.

In the example above, the rows of `f2` are indexed by the `:x1` states
`[:off, :on]`, and the columns are indexed by the `:x2` states
`["low", "mid", "high"]`.

For a factor connected to three variables with cardinalities `2`, `3`, and `2`,
the table must have size `(2, 3, 2)`. The two matrix slices of `f3` along the
third dimension correspond to the `:x3` states `[1, 2]`.

Tables must be finite, nonnegative, and contain at least one positive entry.
They do not need to sum to one.

---

##### Initialization from a Unary Factor

A unary factor can be marked with `initialize = true`. The factor remains an
ordinary factor in the model, but its table is also used to initialize messages
on all edges incident to the connected variable.

When this option is used, the initializing unary factor overrides the variable's
own `probability` initial belief for that variable.

Only unary factors can be used for initialization, and each variable can have at
most one initializing unary factor.

---

## Constructing the Factor Graph

Use [`factorGraph`](@ref) to validate and build the factor graph:

```@example discrete_factor_graph
graph = factorGraph([x1, x2, x3], [f1, f2, f3])

nothing # hide
```

The resulting [`DiscreteFactorGraph`](@ref) stores the validated model, edge
list, and adjacency lists used for message passing.

For quick debugging, the graph structure can be rendered as an SVG with
[`saveGraphFigure`](@ref). This is useful for checking that factor-table
dimensions match the intended variable connections before running inference:

```@example discrete_factor_graph
saveGraphFigure("../dfg.svg", graph)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../dfg.svg"
    type="image/svg+xml"
    aria-label="Discrete factor graph"
    style="width: 47%; height: auto;">
    <a href="../../dfg.svg">Discrete factor graph</a>
  </object>
</div>
```

The rendered SVG is interactive. Hover over variables, factors, and edges to
inspect summary metadata while checking the graph structure. Full discrete
probability, state, and table tooltips are available through the graph figure
API; see [`graphFigure`](@ref) and [`saveGraphFigure`](@ref).

---

## Incremental Construction

The graph can also be built incrementally. Start with an empty
[`DiscreteFactorGraph`](@ref), add variables with [`addVariable!`](@ref), and add
factors with [`addFactor!`](@ref):

```@example discrete_factor_graph
graph = DiscreteFactorGraph()

addVariable!(graph, :x1, 2; label = "x1", states = [:off, :on])
addVariable!(graph, :x2, 3; label = "x2", states = ["low", "mid", "high"])

addFactor!(graph, :x1, [0.8, 0.2]; label = "f1")
addFactor!(graph, :x1, :x2, [1.0 0.2 0.4; 0.1 1.1 0.3]; label = "f2")

nothing # hide
```

Both functions also accept pre-built node objects:

```@example discrete_factor_graph
addVariable!(graph, DiscreteVariable(:x3, 2; label = "x3"))
addFactor!(graph, DiscreteFactor(:x3, [0.55, 0.45]; label = "f3"))

nothing # hide
```

Graph-only construction calls can be interleaved while the model is being
assembled. If inference has already started, graph-only topology changes make
existing inference objects stale. In that case, create a new inference object
before running message passing again.

!!! note "Info"
    This style is useful when a model is assembled from data or loops. The batch
    style remains useful when variables and factors are naturally available as
    collections.

---

## [Tree Factor Graph](@id tree-discrete-factor-graph)

A tree-structured discrete factor graph is a connected discrete factor graph with
no cycles. When a root variable is provided to [`factorGraph`](@ref), the result
is a [`TreeFactorGraph`](@ref) object that supports exact forward-backward
inference on trees:

```@example discrete_factor_graph
f1 = DiscreteFactor(:x1, [0.6, 0.4]; label = "f1")
f2 = DiscreteFactor(:x1, :x2, [1.0 0.2 0.4; 0.2 1.0 0.3]; label = "f2")
f3 = DiscreteFactor(:x2, :x3, [0.8 0.1; 0.1 0.8; 0.4 0.6]; label = "f3")

tree = factorGraph([x1, x2, x3], [f1, f2, f3]; root = :x1)

nothing # hide
```

The returned [`TreeFactorGraph`](@ref) stores the selected root variable,
parent-edge orientation, and the forward/backward edge orders used by tree
inference. The root controls the message order, but not the final exact
marginals after a complete sum-product sweep. The `root` keyword accepts the same
variable references used elsewhere, such as identifiers like `:x3` or labels
like `"x3"`.

If a regular [`DiscreteFactorGraph`](@ref) already exists, it can be converted
explicitly:

```@example discrete_factor_graph
graph = factorGraph([x1, x2, x3], [f1, f2, f3])
tree = treeFactorGraph(graph; root = :x1)

nothing # hide
```

Lookup and print helpers can be called on the tree object and delegate to the
underlying [`DiscreteFactorGraph`](@ref). Graph-only topology changes through the
tree object refresh the stored tree orientation.

---

## Updating Factor Nodes

You can update a factor's table without changing graph topology. The factor can
be selected by integer index or by label:

```@example discrete_factor_graph
updateFactor!(graph, "f1"; table = [0.3, 0.7], initialize = true)
updateFactor!(graph, 1; table = [0.4, 0.6])

nothing # hide
```

The connected variables remain fixed. The updated table size must match the
current factor table size. Keywords that are omitted keep their current values,
so changing factor dimensions requires constructing a new graph.

When `initialize = true` is used, unary factor initialization overrides the
variable's own `probability` initial belief for that variable. The graph-only
`updateFactor!(graph, ...)` method updates only the model data. If inference has
already started, create a new inference object before running message passing
again.

---

## Lookup Functions

Most variable lookup functions accept either a variable identifier or a variable
label:

```@example discrete_factor_graph
variableIndex(graph, :x2)
variableIndex(graph, "x2")
variableDimension(graph, :x2)

nothing # hide
```

State lookup helpers map between state references and one-based table indices:

```@example discrete_factor_graph
stateIndex(graph, :x1, :on)
stateValue(graph, :x2, 2)

nothing # hide
```

Factor references accept integer indices or factor labels:

```@example discrete_factor_graph
factorIndex(graph, "f2")
factorIndex(graph, 2)

nothing # hide
```

Edges can be selected by variable, factor, or their pair:

```@example discrete_factor_graph
edgeIndices(graph; variable = :x1)
edgeIndices(graph; factor = "f2")
edgeIndex(graph; variable = :x1, factor = "f2")

nothing # hide
```

---

## Inspecting the Factor Graph

The print helpers are intended for interactive inspection in the REPL.

[`printModel`](@ref) prints the model data, including variable labels,
identifiers, cardinalities, factor labels, identifiers, connected variables,
table sizes, and initialization flags.

[`printGraph`](@ref) prints a compact table view of the graph. It shows the
factors connected to each variable, the variables connected to each factor, and
the edge identifiers used internally by message passing.

[`printEdges`](@ref) prints one line per edge in the form
`edge_id: variable <-> factor`.

---

## Inference

After constructing the graph, create an inference object and run the selected
discrete belief propagation algorithm.

Use
[Iterative Belief Propagation](@ref iterative-discrete-belief-propagation)
for arbitrary discrete factor graphs, including graphs with cycles. Use
[Forward-Backward Belief Propagation](@ref forward-backward-discrete-belief-propagation)
for exact inference on tree-structured discrete factor graphs.