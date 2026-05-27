# [Forward-Backward Discrete Belief Propagation](@id forward-backward-discrete-belief-propagation)

```@meta
CurrentModule = FactorGraph
```

Forward-backward discrete belief propagation solves tree-structured discrete factor graphs
with one forward pass and one backward pass. On a tree, this schedule is exact: after both
passes, every edge has sent the messages needed to compute the final variable marginals.

The inference object determines which discrete belief propagation variant is run:

- [`sumproduct`](@ref) runs sum-product belief propagation and stores probability
  marginals.
- [`minsum`](@ref) runs min-sum belief propagation for MAP inference and stores MAP
  estimates.

Start from a tree discrete factor graph:

```@example tree_discrete_inference
using FactorGraph # hide

x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
x2 = DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])
x3 = DiscreteVariable(:x3, 2; label = "x3")

f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")
f2 = DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "f2")
f3 = DiscreteFactor(:x2, :x3, [0.9 0.3; 0.2 0.8]; label = "f3")

graph = factorGraph([x1, x2, x3], [f1, f2, f3]; root = :x1)

nothing # hide
```

---

## Sum-Product Inference

Discrete tree inference uses the same [`sumproduct`](@ref) inference state as iterative
discrete belief propagation. Messages and marginals are probability vectors. Create a sum-product
inference state, then use [`forwardBackward!`](@ref) to run the exact tree sweep:

```@example tree_discrete_inference
inference = sumproduct(graph)
forwardBackward!(graph, inference)

nothing # hide
```

Variables created without an initial `probability` use a uniform initial message. A unary
[`DiscreteFactor`](@ref) with `initialize = true` overrides the variable's initial
probability for message initialization.

To inspect results use [`marginal`](@ref) for one variable's probability vector,
[`marginalProbability`](@ref) for one state, or [`marginals`](@ref) to return all stored
variable marginals:

```@example tree_discrete_inference
marginal(graph, inference, :x1)
marginalProbability(graph, inference, :x1, :on)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to one
variable or factor, and [`printMarginal`](@ref) prints one or all variable marginals.

---

## Min-Sum MAP Inference

The min-sum form computes exact MAP estimates on a tree. It is the negative-log form of
discrete max-product belief propagation: factor tables become costs, messages are cost
vectors, and variable beliefs are MAP states.

Create a min-sum inference state with [`minsum`](@ref):

```@example tree_discrete_inference
inference = minsum(graph)
forwardBackward!(graph, inference)

nothing # hide
```

For `minsum`, use [`estimate`](@ref) for one variable and [`estimates`](@ref) for all
variables:

```@example tree_discrete_inference
estimate(graph, inference, :x1)
estimates(graph, inference)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints min-sum messages adjacent to
one variable or factor, and [`printEstimate`](@ref) prints MAP estimates.

---

## Manual Message Updates

Use the manual message functions when you want to control the tree sweep yourself. This
is the same message-passing computation as [`forwardBackward!`](@ref), split into
explicit forward, backward, and result steps.

[`forward!`](@ref) runs the forward pass and [`backward!`](@ref) runs the backward pass.
They do not recompute variable results by themselves, so call [`marginals!`](@ref) for
`sumproduct`, or [`estimates!`](@ref) for `minsum`:

```@example tree_discrete_inference
inference = sumproduct(graph)

forward!(graph, inference)
backward!(graph, inference)
marginals!(graph, inference)

nothing # hide
```

```@example tree_discrete_inference
inference = minsum(graph)

forward!(graph, inference)
backward!(graph, inference)
estimates!(graph, inference)

nothing # hide
```

For edge-level control, use [`forwardStep!`](@ref) and [`backwardStep!`](@ref). Each call
updates one directed tree message and returns the edge that was updated, or `nothing` when
that pass has no remaining scheduled edge:

```@example tree_discrete_inference
inference = sumproduct(graph)

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

When you already know which edge should be updated, pass the edge selection directly to
[`forwardStep!`](@ref) or [`backwardStep!`](@ref):

```@example tree_discrete_inference
inference = sumproduct(graph)

forwardStep!(graph, inference; variable = :x3, factor = "f3")
backwardStep!(graph, inference; variable = :x3, factor = "f3")
marginals!(graph, inference)

nothing # hide
```

Selected-edge calls update only that tree message. They are useful for manual or
incremental schedules where the caller decides which edges need to be recomputed.

Call [`reset!`](@ref) before reusing the tree's automatic step cursor:

```@example tree_discrete_inference
reset!(graph)

nothing # hide
```

---

## Freezing Updates

Freezing keeps selected messages unchanged while the rest of the tree continues updating.
It uses the same inference-state flags as iterative belief propagation.

To freeze all outgoing messages from a factor or variable, use:

```@example tree_discrete_inference
freezeFactor!(graph, inference, "f3")
freezeVariable!(graph, inference, :x3)
```

To freeze both message directions on one edge, use:

```@example tree_discrete_inference
freezeEdge!(graph, inference; variable = :x3, factor = "f3")
```

Resume updates with:

```@example tree_discrete_inference
unfreezeFactor!(graph, inference, "f3")
unfreezeVariable!(graph, inference, :x3)
unfreezeEdge!(graph, inference; variable = :x3, factor = "f3")
```

Inspection helpers return booleans:

```julia
isFrozenFactor(graph, inference, "f3")
isFrozenVariable(graph, inference, :x3)
isFrozenEdge(graph, inference; variable = :x3, factor = "f3")
```

A frozen message is not recomputed. While messages remain frozen, the resulting beliefs
are partial warm-start beliefs rather than exact tree marginals or exact MAP estimates.
Unfreeze the relevant messages and run another complete sweep when exact tree results are
needed.

---

## Dynamic Tree Changes

Use [`addVariable!`](@ref) with a tree inference object when a new variable node must be
added after inference has already started. A new variable is usually followed by at least
one new factor that connects it to the tree. The updated graph must still be a tree:

```@example tree_discrete_inference
addVariable!(graph, inference, :x4, 2; label = "x4", states = [:a, :b])
addFactor!(graph, inference, :x3, :x4, [0.9 0.1; 0.1 0.9]; label = "f4")

forwardBackward!(graph, inference)

nothing # hide
```

Passing the inference object makes this a warm-start update: existing messages and
marginals are preserved, and the message arrays are extended when new factor edges are
added. Adding a factor through the tree view refreshes the tree message order.

Warm-start topology updates extend only the inference object passed to
[`addVariable!`](@ref) or [`addFactor!`](@ref). Any other inference object created from the
same tree before the topology change no longer matches the tree and should be recreated or
updated separately.

Use [`updateFactor!`](@ref) when an existing factor table changes without changing tree
topology:

```@example tree_discrete_inference
updateFactor!(graph, inference; factor = "f4", table = [0.8 0.2; 0.2 0.8])
forwardBackward!(graph, inference)

nothing # hide
```

This reuses the existing inference state as a warm start. The message state still contains
information from previous sweeps, which is useful for dynamic inference. For an independent
solve, update the tree and then create a fresh [`sumproduct`](@ref) object. Because
[`updateFactor!`](@ref) does not change tree topology, other inference objects created
from the same tree still match the tree, but their messages have not been refreshed for
the new factor table. For custom incremental schedules, use [`forwardStep!`](@ref) and
[`backwardStep!`](@ref) on selected edges.

---

## Root Selection

The root determines the forward and backward message order. On a connected tree, it does
not change the final exact marginals after a complete sweep:

```@example tree_discrete_inference
refresh!(graph; root = "x3")
inference = sumproduct(graph)
forwardBackward!(graph, inference)

nothing # hide
```
