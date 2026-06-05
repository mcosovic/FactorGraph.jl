# [Forward-Backward Gaussian Belief Propagation](@id forward-backward-gaussian-belief-propagation)

```@meta
CurrentModule = FactorGraph
```

Forward-backward Gaussian belief propagation solves tree-structured Gaussian
factor graphs with one forward pass and one backward pass. On a tree, this
schedule is exact: after both passes, every edge has sent the messages needed to
compute the final variable marginals or MAP estimates.

The inference object determines which Gaussian belief propagation variant is
run:

* [`moment`](@ref) runs sum-product belief propagation in moment form and stores
  Gaussian marginals.
* [`canonical`](@ref) runs sum-product belief propagation in canonical form and
  stores Gaussian marginals.
* [`minsum`](@ref) runs min-sum belief propagation for MAP inference and stores
  MAP estimates.

The forward-backward message order is the same for these inference states. The
difference is the message representation and the stored variable result:
[`moment`](@ref) and [`canonical`](@ref) store marginals, while [`minsum`](@ref)
stores MAP estimates.

Start from a tree-structured Gaussian factor graph:

```@example tree_gaussian_inference
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

graph = factorGraph(variables, factors; root = :x1)

nothing # hide
```

---

## Sum-Product Moment Inference

The moment form runs sum-product Gaussian belief propagation and stores each
message and marginal as a mean vector and covariance matrix. Create a
moment-form inference state with [`moment`](@ref), then use
[`forwardBackward!`](@ref) to run the exact tree sweep:

```@example tree_gaussian_inference
inference = moment(graph; mean = 0.0, covariance = 1e6)
forwardBackward!(graph, inference)

nothing # hide
```

The `mean` and `covariance` keywords define the default initial belief for
variables created without their own initial `mean` or `covariance`. Scalar
defaults are expanded to each variable dimension.

To inspect results, use [`marginalMean`](@ref) and
[`marginalCovariance`](@ref) when only one field is needed. Use
[`marginal`](@ref) when the full marginal object is needed, or
[`marginals`](@ref) to return all stored variable marginals:

```@example tree_gaussian_inference
marginalMean(graph, inference, :x2)
marginalCovariance(graph, inference, :x2)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to
one variable or factor, and [`printMarginal`](@ref) prints sum-product marginals.

---

## Sum-Product Canonical Inference

The canonical form computes the same exact sum-product Gaussian marginals as the
moment form, but represents messages with information vectors and precision
matrices. Products of Gaussian messages are represented internally as sums of
information vectors and precision matrices.

Create a canonical-form inference state with [`canonical`](@ref), then use
[`forwardBackward!`](@ref) to run the exact tree sweep:

```@example tree_gaussian_inference
inference = canonical(graph; mean = 0.0, covariance = 1e6)
forwardBackward!(graph, inference)

nothing # hide
```

Marginals can still be inspected as means and covariances with the same public
lookup functions used for the moment form.

---

## [Min-Sum MAP Inference](@id forward-backward-gaussian-min-sum-belief-propagation)

The min-sum form computes exact MAP estimates on a tree. It is the negative-log
quadratic form of Gaussian max-product belief propagation: factor potentials
become least-squares costs, messages are quadratic costs, and variable beliefs
are MAP estimates.

This form optimizes the most likely state instead of computing marginal
covariances.

Create a min-sum inference state with [`minsum`](@ref), then use
[`forwardBackward!`](@ref) to run the exact tree sweep:

```@example tree_gaussian_inference
inference = minsum(graph; mean = 0.0, covariance = 1e6)
forwardBackward!(graph, inference)

nothing # hide
```

The `mean` and `covariance` keywords define the default initial belief for
variables created without their own initial `mean` or `covariance`. For
[`minsum`](@ref), this belief is converted to an initial quadratic cost message.

Use [`estimate`](@ref) for one variable and [`estimates`](@ref) for all
variables:

```@example tree_gaussian_inference
estimate(graph, inference, :x2)
estimates(graph, inference)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to
one variable or factor, and [`printEstimate`](@ref) prints min-sum MAP estimates.

---

## Manual Message Updates

Use the manual message functions when you want to control the tree sweep
yourself. They run the same message-passing computation as
[`forwardBackward!`](@ref), but split it into explicit forward, backward, and
result recomputation steps.

[`forward!`](@ref) runs the forward pass, and [`backward!`](@ref) runs the
backward pass. They do not recompute variable results by themselves, so call
[`marginals!`](@ref) for [`moment`](@ref) and [`canonical`](@ref), or
[`estimates!`](@ref) for [`minsum`](@ref):

```@example tree_gaussian_inference
inference = moment(graph)

forward!(graph, inference)
backward!(graph, inference)
marginals!(graph, inference)

nothing # hide
```

For [`minsum`](@ref), recompute MAP estimates after the message passes:

```@example tree_gaussian_inference
inference = minsum(graph)

forward!(graph, inference)
backward!(graph, inference)
estimates!(graph, inference)

nothing # hide
```

For edge-level control, use [`forwardStep!`](@ref) and [`backwardStep!`](@ref).
Each call updates one directed tree message and returns the edge that was
updated, or `nothing` when that pass has no remaining scheduled edge. These
functions work with [`moment`](@ref), [`canonical`](@ref), and [`minsum`](@ref):

```@example tree_gaussian_inference
inference = canonical(graph)

while true
    edge = forwardStep!(graph, inference)
    edge === nothing && break
end

while true
    edge = backwardStep!(graph, inference)
    edge === nothing && break
end

marginals!(graph, inference)

nothing # hide
```

When you already know which edge should be updated, pass the edge selection
directly to [`forwardStep!`](@ref) or [`backwardStep!`](@ref):

```@example tree_gaussian_inference
inference = minsum(graph)

forwardStep!(graph, inference; variable = :x3, factor = "f3")
backwardStep!(graph, inference; variable = :x3, factor = "f3")
estimates!(graph, inference)

nothing # hide
```

Selected-edge calls update only the requested tree message. They are useful for
manual or incremental schedules where the caller decides which edges need to be
recomputed.

Call [`reset!`](@ref) before reusing the tree's automatic step cursor:

```@example tree_gaussian_inference
reset!(graph)

nothing # hide
```

---

## Freezing Updates

Freezing keeps selected messages unchanged while the rest of the tree continues
updating. It uses the same inference-state flags as iterative belief propagation.

To freeze all outgoing messages from a factor or variable, use:

```@example tree_gaussian_inference
freezeFactor!(graph, inference, "f3")
freezeVariable!(graph, inference, :x3)

nothing # hide
```

To freeze both message directions on one edge, use:

```@example tree_gaussian_inference
freezeEdge!(graph, inference; variable = :x3, factor = "f3")

nothing # hide
```

Resume updates by unfreezing the same selections:

```@example tree_gaussian_inference
unfreezeFactor!(graph, inference, "f3")
unfreezeVariable!(graph, inference, :x3)
unfreezeEdge!(graph, inference; variable = :x3, factor = "f3")

nothing # hide
```

Inspection helpers return booleans and use the same references:

```julia
isFrozenFactor(graph, inference, "f3")
isFrozenVariable(graph, inference, :x3)
isFrozenEdge(graph, inference; variable = :x3, factor = "f3")
```

A frozen message is not recomputed. While messages remain frozen, the resulting
beliefs are partial warm-start beliefs rather than exact tree marginals or exact
MAP estimates. Unfreeze the relevant messages and run another complete sweep
when exact tree results are needed.

---

## Dynamic Tree Changes

Use [`addVariable!`](@ref) with a tree inference object when a new variable node
must be added after inference has started. A new variable is usually followed by
at least one new factor that connects it to the tree. The updated graph must
still be a tree:

```@example tree_gaussian_inference
addVariable!(graph, inference, :x4, 1; label = "x4")
addFactor!(graph, inference, :x3, :x4, 0.0, [1.0 -1.0], 0.2; label = "f4")

forwardBackward!(graph, inference)

nothing # hide
```

Passing the inference object makes this a warm-start update. Existing messages
are preserved, and the message arrays are extended when new factor edges are
added.

For [`moment`](@ref) and [`canonical`](@ref), existing marginals are also
preserved, and a new variable receives an initial belief from its own `mean` and
`covariance` or from the inference defaults. For [`minsum`](@ref), this belief
is converted to an initial quadratic cost message and MAP estimate.

Warm-start topology updates extend only the inference object passed to
[`addVariable!`](@ref) or [`addFactor!`](@ref). Any other inference object
created from the same tree before the topology change no longer matches the tree
and should be recreated or updated separately.

For custom incremental schedules, use [`forwardStep!`](@ref) and
[`backwardStep!`](@ref) on selected edges:

```@example tree_gaussian_inference
forwardStep!(graph, inference; variable = :x4, factor = "f4")
backwardStep!(graph, inference; variable = :x4, factor = "f4")
estimates!(graph, inference)

nothing # hide
```

Selected-edge calls update only the requested tree messages. They are useful
when the caller decides which edges need to be recomputed.

Use [`updateFactor!`](@ref) when an existing factor's mean, coefficient, or
covariance changes without changing tree topology:

```@example tree_gaussian_inference
updateFactor!(graph, inference; factor = "f4", mean = 0.1)
forwardBackward!(graph, inference)

nothing # hide
```

This reuses the existing inference state as a warm start. The message state still
contains information from previous sweeps, which is useful for dynamic
inference. For an independent solve, update the tree and then create a fresh
[`moment`](@ref), [`canonical`](@ref), or [`minsum`](@ref) object.

Because [`updateFactor!`](@ref) does not change tree topology, other inference
objects created from the same tree still structurally match the tree. However,
their messages have not been refreshed for the new factor parameters.

---

## Root Selection

The root determines the forward and backward message order. On a connected tree,
it does not change the final exact marginals or MAP estimates after a complete
sweep:

```@example tree_gaussian_inference
refresh!(graph; root = "x3")
inference = canonical(graph)
forwardBackward!(graph, inference)

nothing # hide
```