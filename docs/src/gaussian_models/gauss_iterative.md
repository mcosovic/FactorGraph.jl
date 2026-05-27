# [Iterative Belief Propagation](@id gaussian-sum-product-belief-propagation)

```@meta
CurrentModule = FactorGraph
```

Iterative Gaussian belief propagation solves arbitrary Gaussian factor graphs by repeated
message updates on a [`GaussianFactorGraph`](@ref). Each iteration updates messages from
factor nodes to variable nodes and from variable nodes back to factor nodes. Scheduling,
damping, freezing, and warm-start graph updates are available for moment, canonical, and
min-sum inference states.

The inference object determines which Gaussian belief propagation variant is run:

- [`moment`](@ref) runs sum-product belief propagation in moment form and stores
  Gaussian marginals.
- [`canonical`](@ref) runs sum-product belief propagation in canonical form and stores
  Gaussian marginals.
- [`minsum`](@ref) runs min-sum belief propagation for MAP inference and stores MAP
  estimates.

Use `moment` or `canonical` when posterior means and covariances are needed. Use `minsum`
when only the MAP estimate is needed.

Start from a Gaussian factor graph:

```@example iterative_gaussian_inference
using FactorGraph # hide

variables = [
    GaussianVariable(:x1, 1),
    GaussianVariable(:x2, 2),
    GaussianVariable(:x3, 1),
]

factors = [
    GaussianFactor(:x1, 1.0, 1.0, 0.5; label = "f1"),
    GaussianFactor(:x2, [2.0, -0.5], [1.0 0.0; 0.0 1.0], [1, 1]; label = "f2"),
    GaussianFactor(:x1, :x2, [1.0, 0.2], [1.0 1.0 0.3; 0.5 0.2 1.0], [2, 3]; label = "f3"),
    GaussianFactor(:x2, :x3, [0.5, 0.1], [1.0 0.5 1.0; 0.3 1.0 0.4], [2, 3]; label = "f4"),
    GaussianFactor(:x3, :x1, 0.5, [1.0 -1.0], 0.6; label = "f5"),
]

graph = factorGraph(variables, factors)

nothing # hide
```

The example graph contains a cycle, so it uses the general iterative schedule rather than
the specialized tree schedule.

---

## Sum-Product Moment Inference

The moment form runs sum-product Gaussian belief propagation and stores each message and marginal as a
mean vector and covariance matrix. Create a moment-form inference state with
[`moment`](@ref), then use [`gbp!`](@ref) to run iterative updates. Each iteration updates
factor-to-variable messages, variable-to-factor messages, and then the stored variable
result. The result is a set of Gaussian marginals.

```@example iterative_gaussian_inference
inference = moment(graph; mean = 0.0, covariance = 1e6)
gbp!(graph, inference; iterations = 50, tolerance = 1e-6)

nothing # hide
```

The `mean` and `covariance` keywords define the default initial belief for variables
created without their own initial `mean` or `covariance`. Scalar defaults are expanded to
each variable dimension.

If `iterations` is omitted, [`gbp!`](@ref) runs at most 10 iterations. The `iterations`
keyword always acts as the maximum number of iterations. The tolerance is applied to
[`maxMarginalChange`](@ref). If `tolerance` is omitted, no early stopping is used.

To inspect results use [`marginalMean`](@ref) and [`marginalCovariance`](@ref) when only
one field is needed. Use [`marginal`](@ref) when the full marginal object is needed, or
[`marginals`](@ref) to return all stored variable marginals:

```@example iterative_gaussian_inference
marginalMean(graph, inference, :x1)
marginalCovariance(graph, inference, :x1)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to one
variable or factor, and [`printMarginal`](@ref) prints sum-product marginals.

---

## Sum-Product Canonical Inference

The canonical form computes the same sum-product Gaussian marginals as moment form, but
represents messages with information vectors and precision matrices.
Products of Gaussian messages are represented internally as sums of information vectors
and precision matrices.
Create a canonical-form inference state with [`canonical`](@ref):

```@example iterative_gaussian_inference
inference = canonical(graph; mean = 0.0, covariance = 1e6)
gbp!(graph, inference; iterations = 50, tolerance = 1e-6)
nothing # hide
```

Marginals can still be inspected as means and covariances with the same public lookup
functions used for moment form.

---

## [Min-Sum MAP Inference](@id gaussian-min-sum-belief-propagation)

The min-sum form computes MAP estimates. It is the negative-log quadratic form of Gaussian
max-product belief propagation: factor potentials become least-squares costs, messages are
quadratic costs, and variable beliefs are MAP estimates. This form optimizes the most
likely state instead of computing marginal covariances. The min-sum inference state stores
[`QuadraticMessage`](@ref) objects, which represent quadratic costs rather than Gaussian
marginals. Use moment or canonical when Gaussian marginals are needed.

Create a min-sum inference state with [`minsum`](@ref):

```@example iterative_gaussian_inference
inference = minsum(graph; mean = 0.0, covariance = 1e6)
gbp!(graph, inference; iterations = 50, tolerance = 1e-6)
nothing # hide
```

The `mean` and `covariance` keywords define the default initial belief for variables
created without their own initial `mean` or `covariance`. For `minsum`, this belief is
converted to an initial quadratic cost message.

Each iteration updates factor-to-variable messages, variable-to-factor messages, and then
the stored variable result. For `minsum`, the result is a set of MAP estimates.

The tolerance is applied to [`maxEstimateChange`](@ref).

For `minsum`, use [`estimate`](@ref) for one variable and [`estimates`](@ref) for all
variables:

```@example iterative_gaussian_inference
estimate(graph, inference, :x1)
estimates(graph, inference)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to one
variable or factor, and [`printEstimate`](@ref) prints min-sum MAP estimates.

---

## Manual Message Updates

Use the manual message functions when you want to control the iteration loop yourself.
This is the same message-passing computation as [`gbp!`](@ref), split into explicit
message and result steps.

[`messages!`](@ref) runs one full message step: factor-to-variable messages and
variable-to-factor messages. It does not recompute variable results by itself, so call
[`marginals!`](@ref) for `moment` and `canonical`, or [`estimates!`](@ref) for `minsum`:

```@example iterative_gaussian_inference
inference = canonical(graph)

for _ in 1:20
    messages!(graph, inference)
end

marginals!(graph, inference)

nothing # hide
```

```@example iterative_gaussian_inference
inference = minsum(graph)

for _ in 1:20
    messages!(graph, inference)
end

estimates!(graph, inference)

nothing # hide
```

For directional control, call the factor-to-variable and variable-to-factor passes
explicitly:

```@example iterative_gaussian_inference
inference = canonical(graph)

for _ in 1:20
    factorToVariableMessages!(graph, inference)
    variableToFactorMessages!(graph, inference)
end

marginals!(graph, inference)

nothing # hide
```

The directional functions expose the two message passes separately. They are useful when
experimenting with custom schedules. For normal inference, [`gbp!`](@ref) is the simpler
entry point.

---

## Convergence Diagnostics

Use [`maxMessageChange`](@ref) to compare message states. For stored variable results,
use [`maxMarginalChange`](@ref) with `moment` or `canonical`, and
[`maxEstimateChange`](@ref) with `minsum`.

```@example iterative_gaussian_inference
inference = canonical(graph)

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

The message change includes both factor-to-variable and variable-to-factor messages.
Use [`maxFactorMessageChange`](@ref) or [`maxVariableMessageChange`](@ref) to
inspect one direction.

---

## Sequential Schedule

Sequential scheduling is the default update order for [`gbp!`](@ref) and
[`messages!`](@ref). Factor-to-variable messages are updated in place, then
variable-to-factor messages are updated in place, so later updates in the same pass can
use values computed earlier in that pass.

You can request it explicitly with `schedule = :sequential`:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 20, schedule = :sequential)

nothing # hide
```

For a manual loop, call [`messages!`](@ref) directly; sequential scheduling is the
default:

```@example iterative_gaussian_inference
inference = moment(graph)

for _ in 1:20
    messages!(graph, inference)
    marginals!(graph, inference)
end

nothing # hide
```

For directional control, call the two sequential message passes directly:

```@example iterative_gaussian_inference
inference = moment(graph)

for _ in 1:20
    factorToVariableMessages!(graph, inference)
    variableToFactorMessages!(graph, inference)
    marginals!(graph, inference)
end

nothing # hide
```

---

## Flooding Schedule

Use [`gbp!`](@ref) with `schedule = :flooding` when new messages should be written into
a next-message buffer before being committed:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 20, schedule = :flooding)

nothing # hide
```

For a manual flooding loop, create a flooding schedule and pass it to [`messages!`](@ref):

```@example iterative_gaussian_inference
inference = moment(graph)
schedule = floodingSchedule(graph, inference)

for _ in 1:20
    messages!(graph, inference, schedule)
    marginals!(graph, inference)
end

nothing # hide
```

The returned schedule is mutable step state for this graph and inference object; recreate it
after changing graph topology.

Flooding scheduling works with `moment`, `canonical`, and `minsum` inference states.
Use [`messages!`](@ref) with a [`FloodingSchedule`](@ref) as the manual flooding step.
There are no separate public flooding calls for the two directions, because committing
one direction before the other would change the schedule semantics.

With `schedule = :flooding`, updates within a pass read from the old message arrays and
commit together.

---

## Residual Scheduling

Residual scheduling updates only the directed messages whose recomputed values differ most
from their current values. This can reduce work on loopy graphs where a small subset of
messages dominates the current error.

Use [`gbp!`](@ref) with `schedule = :residual` to run a complete residual-scheduled loop:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 50, schedule = :residual, updateFraction = 0.2)

nothing # hide
```

The `updateFraction` keyword controls the fraction of directed messages updated in each
residual step. Use `updateCount` instead when you want a fixed number of largest-residual
messages:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 50, schedule = :residual, updateCount = 5)

nothing # hide
```

For a manual residual loop, create a [`ResidualSchedule`](@ref) and pass it to
[`messages!`](@ref):

```@example iterative_gaussian_inference
inference = moment(graph)
schedule = residualSchedule(graph, inference; updateFraction = 0.1)

for _ in 1:20
    messages!(graph, inference, schedule)
    marginals!(graph, inference)
end

nothing # hide
```

The returned schedule is mutable step state for this graph and inference object; recreate it
after changing graph topology.

Residual scheduling works with `moment`, `canonical`, and `minsum` inference states.
Use [`messages!`](@ref) with a [`ResidualSchedule`](@ref) as the manual residual step.
There are no separate public residual calls for complete factor-to-variable or
variable-to-factor passes, because each residual step selects individual directed
messages.

---

## Broadcast Message Computation

The `broadcast` keyword selects grouped message computation for `moment`, `canonical`,
and `minsum`.
Pass it to [`gbp!`](@ref) for a full inference run:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 20, broadcast = true)

nothing # hide
```

For a manual loop, pass the same keyword to [`messages!`](@ref):

```@example iterative_gaussian_inference
inference = moment(graph)

for _ in 1:20
    messages!(graph, inference; broadcast = true)
    marginals!(graph, inference)
end

nothing # hide
```

For manual directional loops, pass `broadcast = true` to each message direction:

```@example iterative_gaussian_inference
inference = moment(graph)

for _ in 1:20
    factorToVariableMessages!(graph, inference; broadcast = true)
    variableToFactorMessages!(graph, inference; broadcast = true)
    marginals!(graph, inference)
end

nothing # hide
```

Grouped computation first accumulates all incoming contributions for a factor or variable,
then forms each outgoing message by removing the contribution from the target edge. It can
be combined with sequential or flooding schedules. It is not part of residual scheduling,
because residual scheduling selects individual directed messages by residual rather than
running complete grouped message passes.

---

## Damping

Damping blends a new factor-to-variable message with the previous message. It is often
useful on loopy graphs where undamped message updates may oscillate.

Enable global damping through [`gbp!`](@ref):

```@example iterative_gaussian_inference
inference = minsum(graph)
gbp!(graph, inference; iterations = 50, damping = true, prob = 1.0, alpha = 0.35)

nothing # hide
```

The `alpha` keyword is the weight of the previous message. The `prob` keyword is the
probability that damping is applied on an eligible update.

Damping uses the same `damping`, `prob`, and `alpha` keywords with sequential, flooding,
and residual schedules. It is applied to factor-to-variable message updates.

The same damping keywords can be used with the manual [`messages!`](@ref) API:

```@example iterative_gaussian_inference
inference = canonical(graph)

for _ in 1:20
    messages!(graph, inference; damping = true, prob = 1.0, alpha = 0.35)
end

marginals!(graph, inference)

nothing # hide
```

When running the two message directions manually, pass damping keywords to
[`factorToVariableMessages!`](@ref):

```@example iterative_gaussian_inference
inference = canonical(graph)

for _ in 1:20
    factorToVariableMessages!(graph, inference; damping = true, prob = 1.0, alpha = 0.35)
    variableToFactorMessages!(graph, inference)
end

marginals!(graph, inference)

nothing # hide
```

Use [`dampEdges!`](@ref) to damp selected edges:

```@example iterative_gaussian_inference
dampEdges!(graph, inference; variable = :x1, factor = "f3", alpha = 0.35)

nothing # hide
```

Use [`undampEdges!`](@ref) to remove damping from the selected edges:

```@example iterative_gaussian_inference
undampEdges!(graph, inference; variable = :x1, factor = "f3")

nothing # hide
```

---

## Freezing Updates

Freezing keeps selected messages unchanged while the rest of the graph continues updating.

To freeze all outgoing messages from a factor or variable, use:

```@example iterative_gaussian_inference
freezeFactor!(graph, inference, "f3")
freezeVariable!(graph, inference, :x2)
```

To freeze both message directions on one edge, use:

```@example iterative_gaussian_inference
freezeEdge!(graph, inference; variable = :x1, factor = "f3")
```

Resume updates with:

```@example iterative_gaussian_inference
unfreezeFactor!(graph, inference, "f3")
unfreezeVariable!(graph, inference, :x2)
unfreezeEdge!(graph, inference; variable = :x1, factor = "f3")
```

Inspection helpers return booleans:

```julia
isFrozenFactor(graph, inference, "f3")
isFrozenVariable(graph, inference, :x2)
isFrozenEdge(graph, inference; variable = :x1, factor = "f3")
```

---

## Dynamic Graph Changes

Use [`addVariable!`](@ref) with an inference object when a new variable node must be added
after inference has already started. A new variable is usually followed by at least one new
factor that connects it to the graph:

```@example iterative_gaussian_inference
addVariable!(graph, inference, :x4, 1; label = "x4")
addFactor!(graph, inference, :x3, :x4, 0.4, [1.0 -1.0], 0.2; label = "f6")

gbp!(graph, inference; iterations = 50)
```

Passing the inference object makes this a warm-start update: existing messages are
preserved, and the message arrays are extended when new factor edges are added. For
`moment` and `canonical`, existing marginals are also preserved and a new variable receives
an initial belief from its own `mean`/`covariance` or from the inference defaults. For
`minsum`, this belief is converted to an initial quadratic cost message and MAP estimate.

Warm-start topology updates extend only the inference object passed to
[`addVariable!`](@ref) or [`addFactor!`](@ref). Any other inference object created from the
same graph before the topology change no longer matches the graph and should be recreated
or updated separately.

Use [`updateFactor!`](@ref) when an existing factor's mean, coefficient, or covariance
changes without changing graph topology:

```@example iterative_gaussian_inference
updateFactor!(graph, inference; factor = "f2", mean = [1.8, -0.2], covariance = [0.1, 0.1])
gbp!(graph, inference; iterations = 50)
```

This reuses the existing inference state as a warm start. The message state still contains
information from previous iterations, which is useful for dynamic inference. For an
independent solve, update the graph and then create a fresh `moment`, `canonical`, or
`minsum` object.
Because [`updateFactor!`](@ref) does not change graph topology, other inference objects
created from the same graph still match the graph, but their messages have not been
refreshed for the new factor parameters.
