# [Gaussian Factor Graph](@id gaussian-factor-graph)

```@meta
CurrentModule = FactorGraph
```

A Gaussian factor graph represents a model with continuous variables and local
Gaussian potentials. It is a bipartite graph in which variable nodes store
continuous state blocks, factor nodes store local Gaussian potentials, and each
edge connects a factor to a variable that appears in that factor.

The full model factorizes as a product of local potentials:
```math
p(\mathbf{x}) \propto
\prod_i \psi_i(\mathbf{x}_{\mathcal{V}_i}).
```

Each local potential has the form:
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
factor `covariance` matrix. The vector ``\mathbf{x}_{\mathcal{V}_i} \in \mathbb{R}^{n_i}``
is the stacked state vector for the variables connected to factor ``i``.

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

The identifier, such as `:x1`, is the compact programmatic reference used by
factors. It can be either a `Symbol` or an `Integer`.

The `label` is a `String` used for printed output and lookup functions. Factor
definitions may also refer to variables by label.

Gaussian variables may optionally name the entries of their continuous state
vector with the `components` keyword. Component references can be integers,
symbols, or strings. In the example above, the components of `x2` are `[1, 2]`,
while the components of `x3` are `[:position, :velocity]`.

Components are used for lookup and inspection only. Their order determines the
order of the coefficient matrix columns associated with that variable.

---

##### Dimension

The second argument defines the dimension of the variable node. In this example,
`x1` is scalar, while `x2` and `x3` are vector-valued variables with dimension 2.

---

##### Initial Belief

A variable may optionally define its initial belief using the `mean` and
`covariance` keywords. This belief is used only to initialize belief propagation
messages. It is not a substitute for a unary factor when the prior should be part
of the probabilistic model.

!!! note "Info"
    Whenever covariance is accepted as input, it can be given as a scalar,
    vector, or matrix. For Gaussian factors and inference defaults, a scalar is
    expanded to an isotropic covariance matrix. A vector is internally converted
    into a diagonal covariance matrix, while a matrix is used as a full
    covariance matrix.

If a variable has no initial belief, [`moment`](@ref), [`canonical`](@ref), and
[`minsum`](@ref) use `mean = 0.0` and `covariance = 1e6` by default. These scalar
defaults are expanded to every dimension of each variable node: the mean becomes
a vector filled with `0.0`, and the covariance becomes `1e6` times the identity
matrix of the variable dimension.

For [`minsum`](@ref), this belief is converted to an initial quadratic cost
message. The defaults can be changed by passing the `mean` and `covariance`
keyword arguments to [`moment`](@ref), [`canonical`](@ref), or [`minsum`](@ref).

---

## Factor Nodes

Factor nodes are created with [`GaussianFactor`](@ref):

```@example general_gaussian_factor_graph
z1 = 0.2
H1 = 1.0
Σ1 = 0.1
f1 = GaussianFactor(x1, z1, H1, Σ1; label = "f1", initialize = true)

z2 = [0.4, 0.2]
H2 = [ 1.0 0.5 0.3 0.2 0.1; 0.4 1.2 0.8 0.3 0.2]
Σ2 = [1.0, 1.0]
f2 = GaussianFactor(x1, x2, x3, z2, H2, Σ2; label = "f2")

z3 = [0.4, 0.2]
H3 = [0.1 0.0; 0.0 0.1]
Σ3 = [1.5, 2.5]
f3 = GaussianFactor(x3, z3, H3, Σ3; label = "f3")

nothing # hide
```

---

##### Variable References

The first arguments define the variable nodes connected to the factor node. Their
order is important and must match the column ordering of the coefficient matrix.

A variable can be referenced by identifier, such as `:x1` or `1`, by label, such
as `"x1"`, or by passing the [`GaussianVariable`](@ref) node itself.

---

##### Mean

Each Gaussian factor stores the observation or target vector ``\mathbf{z}_i`` in
its `mean` field. A one-dimensional factor can use a scalar mean, while a
higher-dimensional factor uses a vector mean whose length defines the factor
dimension.

---

##### Coefficient Matrix

Each Gaussian factor stores the coefficient matrix ``\mathbf{H}_i`` in its
`coefficient` field. A scalar can be used for a one-dimensional factor. All
higher-dimensional coefficient inputs must be explicit matrices.

The number of coefficient rows must match the factor mean dimension. The
coefficient columns are ordered according to the connected variable arguments,
and their total count must match the sum of the corresponding variable
dimensions.

---

##### Covariance Matrix

Each Gaussian factor stores the covariance matrix ``\mathbf{\Sigma}_i`` in its
`covariance` field. The covariance input can be given as a scalar, vector, or
matrix. A scalar is expanded to an isotropic covariance matrix matching the
factor mean dimension, a vector is converted into a diagonal covariance matrix,
and a matrix is used as a full covariance matrix.

The covariance matrix must be square and positive definite. Its size must match
the factor mean dimension.

---

##### Initialization from a Unary Factor

A unary factor can be marked with `initialize = true`. The factor remains an
ordinary factor in the model, but its mean and covariance are also used to
initialize messages on all edges incident to the connected variable.

When this option is used, the initializing unary factor overrides the variable's
own `mean` and `covariance` initial belief, as well as the defaults passed to
[`moment`](@ref), [`canonical`](@ref), or [`minsum`](@ref) for that variable.

Only unary factors can be used for initialization, and each variable can have at
most one initializing unary factor.

---

## Constructing the Factor Graph

Use [`factorGraph`](@ref) to validate and build the factor graph:

```@example general_gaussian_factor_graph
graph = factorGraph([x1, x2, x3], [f1, f2, f3])

nothing # hide
```

The resulting [`GaussianFactorGraph`](@ref) stores the validated model, edge
list, and adjacency lists used for message passing.

For quick debugging, the graph structure can be rendered as an SVG with
[`saveGraphFigure`](@ref). This is useful for checking that factors are connected
to the intended variables before running inference:

```@example general_gaussian_factor_graph
saveGraphFigure("../gfg.svg", graph)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../gfg.svg"
    type="image/svg+xml"
    aria-label="Gaussian factor graph"
    style="width: 47%; height: auto;">
    <a href="../../gfg.svg">Gaussian factor graph</a>
  </object>
</div>
```

The rendered SVG is interactive. Hover over variables, factors, and edges to
inspect summary metadata while checking the graph structure. Full Gaussian
parameter tooltips are available through the graph figure API; see
[`graphFigure`](@ref) and [`saveGraphFigure`](@ref).

---

## Incremental Construction

The graph can also be built incrementally. Start with an empty
[`GaussianFactorGraph`](@ref), add variables with [`addVariable!`](@ref), and add
factors with [`addFactor!`](@ref):

```@example general_gaussian_factor_graph
graph = GaussianFactorGraph()

addVariable!(graph, :x1, 1; label = "x1")
addVariable!(graph, :x2, 2; label = "x2")

z1 = 0.2
H1 = 1.0
Σ1 = 0.1
addFactor!(graph, :x1, z1, H1, Σ1; label = "f1")

z2 = [0.4, 0.2]
H2 = [1.0 0.5 0.3; 0.4 1.2 0.8]
Σ2 = [1.0, 1.0]
addFactor!(graph, :x1, :x2, z2, H2, Σ2; label = "f2")

nothing # hide
```

Both functions also accept pre-built node objects:

```@example general_gaussian_factor_graph
addVariable!(graph, GaussianVariable(:x3, 1; label = "x3"))
addFactor!(graph, GaussianFactor(:x3, 0.1, 1.0, 0.5; label = "f3"))

nothing # hide
```

Graph-only construction calls can be interleaved while the model is being
assembled. If inference has already started, graph-only topology changes make
existing inference objects stale. To keep a running inference state as a warm
start, use the inference-aware forms `addVariable!(graph, inference, ...)` and
`addFactor!(graph, inference, ...)`, which extend the graph and the inference
object together.

!!! note "Info"
    This style is useful when a model is assembled from data or loops. The batch
    style remains useful when variables and factors are naturally available as
    collections.

---

## [Tree Factor Graph](@id tree-gaussian-factor-graph)

A tree-structured Gaussian factor graph is a connected Gaussian factor graph with
no cycles. When a root variable is provided to [`factorGraph`](@ref), the result
is a [`TreeFactorGraph`](@ref) object that supports forward-backward message
passing and exact marginal computation on trees.

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

The returned [`TreeFactorGraph`](@ref) stores a tree-structured representation of
the Gaussian factor graph. It contains the selected root variable, parent-edge
orientation, and the forward/backward edge orders used by tree inference.

The root controls the order of forward and backward messages, but not the final
exact marginals. The `root` keyword accepts the same variable references used
elsewhere, such as identifiers like `:x3` or labels like `"x3"`.

If a regular [`GaussianFactorGraph`](@ref) already exists, it can be converted
explicitly:

```@example tree_gaussian_factor_graph
graph = factorGraph(variables, factors)
tree = treeFactorGraph(graph; root = :x1)

nothing # hide
```

---


## Updating Factor Nodes

You can update a factor's mean, coefficient, or covariance without changing the
graph topology. The factor can be selected by integer index or by label:

```@example general_gaussian_factor_graph
updateFactor!(graph, "f1"; mean = 0.30, coefficient = 1.2, covariance = 0.20)
updateFactor!(graph, 1; mean = 0.35)

nothing # hide
```

The connected variables remain fixed. The updated mean length, coefficient matrix
size, and covariance matrix size must match the current factor. Keywords that
are omitted keep their current values, so changing factor dimensions requires
constructing a new graph.

When `initialize = true` is used, unary factor initialization overrides the
variable's own `mean` and `covariance` initial belief. The graph-only
`updateFactor!(graph, ...)` method updates only the model data. Use
`updateFactor!(graph, inference, ...)` when a running inference object should
also reset warm-start messages from the initializing unary factor.

---

## Lookup Functions

Most variable lookup functions accept either a variable identifier or a variable
label:

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
edgeIndices(graph; factor = "f2")
edgeIndex(graph; variable = :x1, factor = "f2")

nothing # hide
```

For matrix inspection, [`coefficientBlock`](@ref) returns the submatrix of one
factor that corresponds to one variable, while [`coefficientBlocks`](@ref)
returns all blocks in the factor's variable order:

```@example general_gaussian_factor_graph
coefficientBlock(graph; factor = "f2", variable = :x2)
coefficientBlocks(graph, 2)
coefficientBlocks(graph, "f2")

nothing # hide
```

---

## Inspecting the Factor Graph

The print helpers are intended for interactive inspection in the REPL.

[`printModel`](@ref) prints the model data, including variable labels,
identifiers, dimensions, factor labels, identifiers, connected variables, and
matrix sizes.

[`printGraph`](@ref) prints a compact table view of the graph. It shows the
factors connected to each variable, the variables connected to each factor, and
the edge identifiers used internally by message passing.

[`printEdges`](@ref) prints one line per edge in the form
`edge_id: variable <-> factor`.

---

## Inference

After constructing the graph, create an inference object and run the selected
Gaussian belief propagation algorithm.

Use
[Iterative Belief Propagation](@ref gaussian-sum-product-belief-propagation)
for arbitrary Gaussian factor graphs, including graphs with cycles. Use
[Forward-Backward Belief Propagation](@ref forward-backward-gaussian-belief-propagation)
for exact inference on tree-structured Gaussian factor graphs.

---
