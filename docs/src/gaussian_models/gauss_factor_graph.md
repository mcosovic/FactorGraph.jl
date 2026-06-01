# [Factor Graph](@id gaussian-factor-graph)

```@meta
CurrentModule = FactorGraph
```

A Gaussian factor graph represents a model with continuous variables and local Gaussian
potentials. It is bipartite: variable nodes store continuous state blocks, factor nodes
store local Gaussian potentials, and edges connect each factor to the variables that appear
in it.

The full model factorizes as a product of local potentials:

```math
p(\mathbf{x}) \propto
\prod_i \psi_i(\mathbf{x}_{\mathcal{V}_i}).
```

In a Gaussian factor graph, each local potential has the form:

```math
\psi_i(\mathbf{x}_{\mathcal{V}_i})
=
\mathcal{N}
\left(
    \mathbf{z}_i ;
    \mathbf{H}_i \mathbf{x}_{\mathcal{V}_i},
    \mathbf{\Sigma}_i
\right).
```

Equivalently:

```math
\psi_i(\mathbf{x}_{\mathcal{V}_i})
\propto
\exp\left(
-\frac{1}{2}
\left(
\mathbf{z}_i -
\mathbf{H}_i \mathbf{x}_{\mathcal{V}_i}
\right)^\top
\mathbf{\Sigma}_i^{-1}
\left(
\mathbf{z}_i -
\mathbf{H}_i \mathbf{x}_{\mathcal{V}_i}
\right)
\right).
```

Here, ``\mathbf{z}_i \in \mathbb{R}^{m_i}`` is the observation or target vector stored in
the factor's `mean` field, ``\mathbf{H}_i \in \mathbb{R}^{m_i \times n_i}`` is the factor
`coefficient` matrix, and ``\mathbf{\Sigma}_i \in \mathbb{R}^{m_i \times m_i}`` is the
factor `covariance` matrix. The vector
``\mathbf{x}_{\mathcal{V}_i} \in \mathbb{R}^{n_i}`` is the stacked state vector for the
variables connected to factor ``i``.

---

## Variable Nodes

Variable nodes are created with [`GaussianVariable`](@ref):

```@example general_gaussian_factor_graph
using FactorGraph # hide

x1 = GaussianVariable(:x1, 1; label = "x1", mean = 0.5, covariance = 1e1)
x2 = GaussianVariable(:x2, 2; label = "x2", mean = [1.0, 0.2], covariance = [1e3, 1e2])
x3 = GaussianVariable(:x3, 2; label = "x3", components = [:position, :velocity])

nothing # hide
```

---

##### Identifiers and Labels

The identifier, such as `:x1`, is the compact programmatic reference used by factors. It can
be either a `Symbol` or an `Integer`.

The `label` is a `String` used for printed output and lookup functions. Factor definitions
may also refer to variables by label.

Gaussian variables may optionally name the entries of their continuous state vector with the
`components` keyword. Component references can be integers, symbols, or strings. For `x2`,
the components are `[1, 2]`, while for `x3`, they are `[:position, :velocity]`. Components
are used only for lookup and inspection, while their order determines the order of the
coefficient matrix columns associated with that variable.

---

##### Dimension

The second argument defines the dimension of the variable node. In this example, `x1` is
scalar, while `x2` and `x3` are vector-valued variable nodes with dimension 2.

---

##### Initial Belief

A variable may optionally define its own initial belief using the `mean` and `covariance`
keywords. This belief is used only to initialize belief propagation messages. It is not a
replacement for a unary factor node when the prior should be part of the probabilistic
model.

!!! note "Info"
    Whenever covariance is accepted as input, it can be given as a scalar, vector, or
    matrix. A scalar is used for the one-dimensional case. A vector is internally converted
    into a diagonal covariance matrix, while a matrix is used as a full covariance matrix.

If a variable has no initial belief defined in this way, [`moment`](@ref),
[`canonical`](@ref), and [`minsum`](@ref) use `mean = 0.0` and
`covariance = 1e6` by default. These scalar defaults are expanded to every dimension of
each variable node: the mean becomes a vector filled with `0.0`, and the covariance
becomes `1e6` times the identity matrix of the variable dimension. For [`minsum`](@ref), this
belief is converted to an initial quadratic cost message. The defaults can be changed by
passing the `mean` and `covariance` keyword arguments to [`moment`](@ref),
[`canonical`](@ref), or [`minsum`](@ref).

---

## Factor Nodes

Factor nodes are created with [`GaussianFactor`](@ref):

```@example general_gaussian_factor_graph
z1 = 0.2; H1 = 1.0; Σ1 = 0.1
f1 = GaussianFactor(x1, z1, H1, Σ1; label = "f1", initialize = true)

z2 = [0.4, 0.2]; H2 = [1.0 0.5 0.3; 0.4 1.2 0.8]; Σ2 = [1.0, 1.0]
f2 = GaussianFactor(x1, x2, z2, H2, Σ2; label = "f2")

z3 = [0.4, 0.2]; H3 = [0.1 0.0; 0.0 0.1]; Σ3 = [1.5, 2.5]
f3 = GaussianFactor(x3, z3, H3, Σ3; label = "f3")

nothing # hide
```

---

##### Variable Refs

The first arguments define the variable nodes connected to the factor node. Their order is
important and must match the column ordering of the coefficient matrix. A variable can be
referenced by id, such as `:x1` or `1`, by label, such as `"x1"`, or by passing the
[`GaussianVariable`](@ref) node itself.

---

##### Mean Values

Each Gaussian factor is defined by a mean value that stores ``\mathbf{z}_i``. A
one-dimensional factor can use a scalar mean, while a higher-dimensional factor uses a
vector mean whose length defines the factor dimension.

---

##### Coefficient Matrix

The coefficient matrix stores ``\mathbf{H}_i``. A scalar is used for a one-dimensional
factor. All higher-dimensional coefficient inputs must be explicit matrices.

The number of coefficient rows must match the factor mean dimension. The coefficient columns
are ordered by the connected variable arguments, and their total count must match the sum of
those variable dimensions.

---

##### Covariance Matrix

The covariance matrix stores ``\mathbf{\Sigma}_i``. The covariance input can be given as a
scalar, vector, or matrix. A scalar is used for a one-dimensional factor, a vector is
converted into a diagonal covariance matrix, and a matrix is used as a full covariance
matrix.

The covariance matrix must be square, positive definite, and its size must match the factor
mean dimension.

---

##### Initialization from a Unary Factor

A unary factor can be marked with `initialize = true`. This keeps the factor as an ordinary
factor in the model, but also uses its mean and covariance to initialize messages on all
edges incident to the connected variable.

When this is used, the initializing unary factor overrides the variable's own
`mean`/`covariance` initial belief and the defaults passed to [`moment`](@ref),
[`canonical`](@ref), or [`minsum`](@ref) for that variable.

Only unary factors can use `initialize = true`, and each variable can have at most one
initializing unary factor.

---

## Constructing the Factor Graph

Use [`factorGraph`](@ref) to validate and build the factor graph.

```@example general_gaussian_factor_graph
graph = factorGraph([x1, x2], [f1, f2])

nothing # hide
```

The resulting [`GaussianFactorGraph`](@ref) stores the validated model, edge list, and
adjacency lists used for message passing.

For quick debugging, the graph structure can be rendered as a SVG with
[`saveGraphFigure`](@ref). This is useful for checking that factors connect to the intended
variables before running inference:

```@example general_gaussian_factor_graph
saveGraphFigure("gfg.svg", graph)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <img src="gfg.svg" alt="Gaussian factor graph" style="width: 45%; height: auto;">
</div>
```

---

## Incremental Construction

The graph can also be built incrementally. Start with an empty [`GaussianFactorGraph`](@ref),
add variables with [`addVariable!`](@ref), and add factors with [`addFactor!`](@ref):

```@example general_gaussian_factor_graph
graph = GaussianFactorGraph()

addVariable!(graph, :x1, 1; label = "x1")
addVariable!(graph, :x2, 2; label = "x2")

addFactor!(graph, :x1, 0.2, 1.0, 0.1; label = "f1")
addFactor!(graph, :x1, :x2, [0.4, 0.2], [1.0 0.5 0.3; 0.4 1.2 0.8], [1, 1]; label = "f2")

nothing # hide
```

Both functions also accept pre-built node objects:

```@example general_gaussian_factor_graph
addVariable!(graph, GaussianVariable(:x3, 1; label = "x3"))
addFactor!(graph, GaussianFactor(:x3, 0.1, 1.0, 0.5; label = "f3"))

nothing # hide
```

Graph-only construction calls may be interleaved while the model is being assembled. If
inference has already started, graph-only topology changes make existing inference objects
stale. To keep a running inference state as a warm start, use the inference-aware forms
`addVariable!(graph, inference, ...)` and `addFactor!(graph, inference, ...)`, which extend
the graph and the inference object together.

!!! note "Info"
    This style is useful when a model is assembled from data or loops. The batch style
    remains useful when variables and factors are naturally available as collections.

---

## [Tree Factor Graph](@id tree-gaussian-factor-graph)

A tree-structured Gaussian factor graph is a connected Gaussian factor graph with no cycles.
When a `root` variable is provided to [`factorGraph`](@ref), the result is a
[`TreeFactorGraph`](@ref) object that can run forward-backward message passing and compute
exact marginals on a tree.

```@example tree_gaussian_factor_graph
using FactorGraph # hide

variables = [
    GaussianVariable(:x1, 1; label = "x1"),
    GaussianVariable(:x2, 1; label = "x2"),
    GaussianVariable(:x3, 1; label = "x3"),
]

factors = [
    GaussianFactor(:x1, 1.0, 1.0, 0.1; label = "f1"),
    GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.2; label = "f2"),
    GaussianFactor(:x2, :x3, 0.0, [1.0 -1.0], 0.2; label = "f3"),
]

tree = factorGraph(variables, factors; root = :x1)

nothing # hide
```

The returned [`TreeFactorGraph`](@ref) is a tree-structured representation of a Gaussian
factor graph. It stores the selected root variable, parent-edge orientation, and the
forward/backward edge orders used by tree inference.

The root controls the order of forward and backward messages, but not the final exact
marginals. The `root` keyword accepts the same variable references used elsewhere: ids such
as `:x3` or labels such as `"x3"`.

If a regular [`GaussianFactorGraph`](@ref) already exists, it can still be converted
explicitly:

```@example tree_gaussian_factor_graph
graph = factorGraph(variables, factors)
tree = treeFactorGraph(graph; root = :x1)

nothing # hide
```

---


## Updating Factor Nodes

You can update a factor's mean, coefficient, or covariance without changing graph topology.
The factor can be selected by integer index or by label:

```@example general_gaussian_factor_graph
updateFactor!(graph, "f1"; mean = 0.30, coefficient = 1.2, covariance = 0.20)
updateFactor!(graph, 1; mean = 0.35)

nothing # hide
```

The connected variables stay fixed. The updated mean length, coefficient matrix size, and
covariance matrix size must match the old factor. Keywords that are omitted keep their old
values, so changing factor dimensions still requires constructing a new graph.

When `initialize = true` is used, unary factor initialization overrides the variable's own
`mean`/`covariance` initial belief. The graph-only `updateFactor!(graph, ...)` method only
updates the model data; use `updateFactor!(graph, inference, ...)` when a running inference
object should also reset warm-start messages from the initializing unary factor.

---

## Lookup Functions

Most lookup functions accept either a variable id or a variable label:

```@example general_gaussian_factor_graph
variableIndex(graph, :x2)
variableIndex(graph, "x2")
variableDimension(graph, :x2)

nothing # hide
```

Factor references accept integer indices or factor labels:

```@example general_gaussian_factor_graph
factorIndex(graph, "f2")
factorIndex(graph, 2)

nothing # hide
```

Edges can be selected by variable, factor, or their pair:

```@example general_gaussian_factor_graph
edgeIndices(graph; variable = :x1)
edgeIndex(graph; variable = :x1, factor = "f2")

nothing # hide
```

For matrix inspection, [`coefficientBlock`](@ref) returns the submatrix of one factor that
corresponds to one variable, while [`coefficientBlocks`](@ref) returns all blocks in the
factor's variable order:

```@example general_gaussian_factor_graph
coefficientBlock(graph; factor = "f2", variable = :x2)
coefficientBlocks(graph, 2)
coefficientBlocks(graph, "f2")

nothing # hide
```

---

## Inspecting the Factor Graph

The print helpers are intended for interactive inspection in the REPL.

[`printModel`](@ref) prints the model data: variable labels, ids, dimensions, factor labels,
ids, connected variables, and matrix sizes.

[`printGraph`](@ref) prints a compact table view of the graph. It shows which factors are
connected to each variable, which variables are connected to each factor, and the edge ids
used internally by message passing.

[`printEdges`](@ref) prints one line per edge, using the form
`edge_id: variable <-> factor`.

---

## Inference

After the graph is formed, create an inference object and run one of the
[Gaussian belief propagation algorithms](@ref gaussian-sum-product-belief-propagation).

For a [`TreeFactorGraph`](@ref), the recommended option is to use the
[forward-backward message schedule](@ref forward-backward-gaussian-belief-propagation),
which exploits the tree structure and computes the exact solution with one
forward and one backward pass.
