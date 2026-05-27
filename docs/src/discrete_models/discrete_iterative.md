# [Iterative Discrete Belief Propagation](@id iterative-discrete-belief-propagation)

```@meta
CurrentModule = FactorGraph
```

Iterative discrete belief propagation runs message updates on a
[`DiscreteFactorGraph`](@ref). Each iteration updates messages from factor nodes to
variable nodes and from variable nodes to factor nodes, then recomputes the stored
variable result. This is the general schedule for arbitrary discrete factor graphs,
including graphs with cycles.

The inference object determines which discrete belief propagation variant is run:

- [`sumproduct`](@ref) runs sum-product belief propagation and stores probability
  marginals.
- [`minsum`](@ref) runs min-sum belief propagation for MAP inference and stores MAP
  estimates.

Use `sumproduct` when marginal probabilities are needed. Use `minsum` when only the MAP
assignment is needed.

Start from a discrete factor graph:

```@example iterative_discrete_inference
using FactorGraph # hide

x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
x2 = DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])
x3 = DiscreteVariable(:x3, 2; label = "x3")

f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1", initialize = true)
f2 = DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "f2")
f3 = DiscreteFactor(:x2, :x3, [0.9 0.3; 0.2 0.8]; label = "f3")
f4 = DiscreteFactor(:x3, :x1, [0.7 0.4; 0.3 0.6]; label = "f4")

graph = factorGraph([x1, x2, x3], [f1, f2, f3, f4])

nothing # hide
```

The example graph contains a cycle, so it uses the general iterative schedule rather than
the specialized tree schedule.

---

## Sum-Product Inference

The sum-product form of discrete belief propagation represents messages and
marginals as probability vectors. Create a sum-product inference state with
[`sumproduct`](@ref), then use [`gbp!`](@ref) to run iterative updates:

```@example iterative_discrete_inference
inference = sumproduct(graph)
gbp!(graph, inference; iterations = 100, tolerance = 1e-6)

nothing # hide
```

Each iteration updates factor-to-variable messages, variable-to-factor
messages, and then the stored variable marginals. Variables created without an
initial `probability` use a uniform initial message. A unary [`DiscreteFactor`](@ref)
with `initialize = true` overrides the variable's initial probability for message
initialization.

If `iterations` is omitted, [`gbp!`](@ref) runs at most 10 iterations. The `iterations`
keyword always acts as the maximum number of iterations. The tolerance is applied to
[`maxMarginalChange`](@ref). If `tolerance` is omitted, no early stopping is used.

To inspect results use [`marginal`](@ref) for one variable's probability vector,
[`marginalProbability`](@ref) for one state, or [`marginals`](@ref) to return all
stored variable marginals:

```@example iterative_discrete_inference
marginal(graph, inference, :x1)
marginalProbability(graph, inference, :x1, :on)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to
one variable or factor, and [`printMarginal`](@ref) prints one or all variable
marginals.

---

## Min-Sum MAP Inference

The min-sum form computes MAP estimates. It is the negative-log form of discrete
max-product belief propagation: factor tables become costs, messages are cost vectors,
and variable beliefs are MAP states. This form optimizes the most likely assignment
instead of computing marginal probabilities.

Create a min-sum inference state with [`minsum`](@ref):

```@example iterative_discrete_inference
inference = minsum(graph)
gbp!(graph, inference; iterations = 50, tolerance = 1e-6)

nothing # hide
```

Each iteration updates factor-to-variable messages, variable-to-factor messages, and then
the stored variable result. For `minsum`, the result is a set of MAP estimates. The
tolerance is applied to [`maxEstimateChange`](@ref).

For `minsum`, use [`estimate`](@ref) for one variable and [`estimates`](@ref) for all
variables:

```@example iterative_discrete_inference
estimate(graph, inference, :x1)
estimates(graph, inference)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints min-sum messages adjacent to
one variable or factor, and [`printEstimate`](@ref) prints MAP estimates.

---

## Manual Message Updates

Use the manual message functions when you want to control the iteration loop yourself.
This is the same message-passing computation as [`gbp!`](@ref), split into explicit
message and result steps.

[`messages!`](@ref) runs one full message step: factor-to-variable messages and
variable-to-factor messages. It does not recompute variable results by itself, so call
[`marginals!`](@ref) for `sumproduct`, or [`estimates!`](@ref) for `minsum`:

```@example iterative_discrete_inference
inference = sumproduct(graph)

for _ in 1:20
    messages!(graph, inference)
end

marginals!(graph, inference)

nothing # hide
```

```@example iterative_discrete_inference
inference = minsum(graph)

for _ in 1:20
    messages!(graph, inference)
end

estimates!(graph, inference)

nothing # hide
```

You can also run the two message directions explicitly:

```@example iterative_discrete_inference
inference = sumproduct(graph)

for _ in 1:20
    factorToVariableMessages!(graph, inference)
    variableToFactorMessages!(graph, inference)
end

marginals!(graph, inference)

nothing # hide
```

The directional functions expose the two message passes separately. They are
useful when experimenting with custom schedules. For normal inference,
[`gbp!`](@ref) is the simpler entry point.

---

## Convergence Diagnostics

Use [`maxMessageChange`](@ref) to compare message states. For stored variable results,
use [`maxMarginalChange`](@ref) with `sumproduct`, and
[`maxEstimateChange`](@ref) with `minsum`.


```@example iterative_discrete_inference
inference = sumproduct(graph)

for _ in 1:100
    oldVarMessages = deepcopy(inference.variableToFactor)
    oldFacMessages = deepcopy(inference.factorToVariable)
    oldMarginals = deepcopy(inference.marginal)

    messages!(graph, inference)
    marginals!(graph, inference)

    if maxMessageChange(graph, inference, oldVarMessages, oldFacMessages) < 1e-6 &&
            maxMarginalChange(graph, inference, oldMarginals) < 1e-6
        break
    end
end

nothing # hide
```

The message change includes both factor-to-variable and variable-to-factor
messages. Use [`maxFactorMessageChange`](@ref) or
[`maxVariableMessageChange`](@ref) to inspect one direction.
The marginal change compares marginal probability vectors.

---

## Sequential Schedule

Sequential scheduling is the default update order for [`gbp!`](@ref) and
[`messages!`](@ref). Factor-to-variable messages are updated in place, then
variable-to-factor messages are updated in place, so later updates in the same pass can
use values computed earlier in that pass.

You can request it explicitly with `schedule = :sequential`:

```@example iterative_discrete_inference
inference = sumproduct(graph)
gbp!(graph, inference; iterations = 20, schedule = :sequential)

nothing # hide
```

For a manual loop, call [`messages!`](@ref) directly; sequential scheduling is the
default:

```@example iterative_discrete_inference
inference = sumproduct(graph)

for _ in 1:20
    messages!(graph, inference)
    marginals!(graph, inference)
end

nothing # hide
```

For directional control, call the two sequential message passes directly:

```@example iterative_discrete_inference
inference = sumproduct(graph)

for _ in 1:20
    factorToVariableMessages!(graph, inference)
    variableToFactorMessages!(graph, inference)
    marginals!(graph, inference)
end

nothing # hide
```

---

## Flooding Schedule

Use [`gbp!`](@ref) with `schedule = :flooding` when new messages should be written into a
next-message buffer before being committed:

```@example iterative_discrete_inference
inference = sumproduct(graph)
gbp!(graph, inference; iterations = 20, schedule = :flooding)

nothing # hide
```

For a manual flooding loop, create a flooding schedule and pass it to [`messages!`](@ref):

```@example iterative_discrete_inference
inference = sumproduct(graph)
schedule = floodingSchedule(graph, inference)

for _ in 1:20
    messages!(graph, inference, schedule)
    marginals!(graph, inference)
end

nothing # hide
```

The returned schedule is mutable step state for this graph and inference object; recreate it
after changing graph topology.

Flooding scheduling works with `sumproduct` and `minsum` inference states.
Use [`messages!`](@ref) with a [`FloodingSchedule`](@ref) as the manual flooding step.
There are no separate public flooding calls for the two directions, because committing
one direction before the other would change the schedule semantics.

With `schedule = :flooding`, updates within a pass read from the old message arrays and
commit together.

---

## Residual Scheduling

Residual scheduling updates only the directed messages whose recomputed probability
vectors differ most from their current values. This can reduce work on loopy graphs where
a small subset of messages dominates the current error.

Use [`gbp!`](@ref) with `schedule = :residual` to run a complete residual-scheduled loop:

```@example iterative_discrete_inference
inference = sumproduct(graph)
gbp!(graph, inference; schedule = :residual, updateFraction = 0.2)

nothing # hide
```

The `updateFraction` keyword controls the fraction of directed messages updated in each
residual step. Use `updateCount` instead when you want a fixed number of largest-residual
messages:

```@example iterative_discrete_inference
inference = sumproduct(graph)
gbp!(graph, inference; schedule = :residual, updateCount = 5)

nothing # hide
```

For a manual residual loop, create a [`ResidualSchedule`](@ref) and pass it to
[`messages!`](@ref):

```@example iterative_discrete_inference
inference = sumproduct(graph)
schedule = residualSchedule(graph, inference; updateFraction = 0.1)

for _ in 1:20
    messages!(graph, inference, schedule)
    marginals!(graph, inference)
end

nothing # hide
```

The returned schedule is mutable step state for this graph and inference object; recreate it
after changing graph topology.

Residual scheduling works with `sumproduct` and `minsum` inference states.
Use [`messages!`](@ref) with a [`ResidualSchedule`](@ref) as the manual residual step.
There are no separate public residual calls for complete factor-to-variable or
variable-to-factor passes, because each residual step selects individual directed
messages.

---

## Damping

Discrete damping mixes a new factor-to-variable message with the previous
message and renormalizes the result. Enable global damping through [`gbp!`](@ref):

```@example iterative_discrete_inference
inference = sumproduct(graph)
gbp!(graph, inference; iterations = 50, damping = true, prob = 1.0, alpha = 0.35)

nothing # hide
```

The `alpha` keyword is the weight of the previous message. The `prob` keyword is the
probability that damping is applied on an eligible update. Damping is applied to
factor-to-variable message updates.

Use [`dampEdges!`](@ref) to damp selected edges:

```@example iterative_discrete_inference
dampEdges!(graph, inference; variable = :x2, factor = "f2", prob = 1.0, alpha = 0.35)
isDampedEdge(graph, inference; variable = :x2, factor = "f2")

nothing # hide
```

Use [`undampEdges!`](@ref) to remove damping from selected edges.

---

## Freezing Updates

Freezing keeps selected messages unchanged while the rest of the graph continues updating.

To freeze all outgoing messages from a factor or variable, use:

```@example iterative_discrete_inference
freezeFactor!(graph, inference, "f2")
freezeVariable!(graph, inference, :x2)
```

To freeze both message directions on one edge, use:

```@example iterative_discrete_inference
freezeEdge!(graph, inference; variable = :x2, factor = "f2")
```

Resume updates with:

```@example iterative_discrete_inference
unfreezeFactor!(graph, inference, "f2")
unfreezeVariable!(graph, inference, :x2)
unfreezeEdge!(graph, inference; variable = :x2, factor = "f2")
```

Inspection helpers return booleans:

```julia
isFrozenFactor(graph, inference, "f2")
isFrozenVariable(graph, inference, :x2)
isFrozenEdge(graph, inference; variable = :x2, factor = "f2")
```

---

## Dynamic Graph Changes

Use [`addVariable!`](@ref) with an inference object when a new variable node must be added
after inference has already started. A new variable is usually followed by at least one new
factor that connects it to the graph:

```@example iterative_discrete_inference
addVariable!(graph, inference, :x4, 2; label = "x4", states = [:a, :b])
addFactor!(graph, inference, :x3, :x4, [0.9 0.1; 0.1 0.9]; label = "f5")

gbp!(graph, inference; iterations = 50)

nothing # hide
```

Messages on new edges are initialized from the connected variables' initial messages,
matching the initialization used by [`sumproduct`](@ref). Existing messages and marginals
are kept.

Warm-start topology updates extend only the inference object passed to
[`addVariable!`](@ref) or [`addFactor!`](@ref). Any other inference object created from the
same graph before the topology change no longer matches the graph and should be recreated
or updated separately.

Use [`updateFactor!`](@ref) when an existing factor table changes without changing graph
topology:

```@example iterative_discrete_inference
updateFactor!(graph, inference; factor = "f1", table = [0.7, 0.3], initialize = true)
gbp!(graph, inference; iterations = 50)

nothing # hide
```

This reuses the existing inference state as a warm start. The message state still contains
information from previous iterations, which is useful for dynamic inference. For an
independent solve, update the graph and then create a fresh [`sumproduct`](@ref) object.
Because [`updateFactor!`](@ref) does not change graph topology, other inference objects
created from the same graph still match the graph, but their messages have not been
refreshed for the new factor table.
