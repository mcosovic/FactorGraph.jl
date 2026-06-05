# [Iterative Gaussian Belief Propagation](@id gaussian-sum-product-belief-propagation)

```@meta
CurrentModule = FactorGraph
```

Iterative Gaussian belief propagation solves arbitrary Gaussian factor graphs by
repeated message updates on a [`GaussianFactorGraph`](@ref). Each iteration
updates messages from factor nodes to variable nodes and from variable nodes
back to factor nodes. Scheduling, damping, freezing, and warm-start graph updates
are available for moment, canonical, and min-sum inference states.

The inference object determines which Gaussian belief propagation variant is
run:

* [`moment`](@ref) runs sum-product belief propagation in moment form and stores
  Gaussian marginals.
* [`canonical`](@ref) runs sum-product belief propagation in canonical form and
  stores Gaussian marginals.
* [`minsum`](@ref) runs min-sum belief propagation for MAP inference and stores
  MAP estimates.

Use [`moment`](@ref) or [`canonical`](@ref) when posterior means and covariances
are needed. Use [`minsum`](@ref) when only the MAP estimate is needed.

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
    GaussianFactor(:x2, [2.0, -0.5], [1.0 0.0; 0.0 1.0], [1.0, 1.0]; label = "f2"),
    GaussianFactor(:x1, :x2, [1.0, 0.2], [1.0 1.0 0.3; 0.5 0.2 1.0], 2.0; label = "f3"),
    GaussianFactor(:x2, :x3, [0.5, 0.1], [1.0 0.5 1.0; 0.3 1.0 0.4], 2.0; label = "f4"),
    GaussianFactor(:x3, :x1, 0.5, [1.0 -1.0], 0.6; label = "f5"),
]

graph = factorGraph(variables, factors)

nothing # hide
```

The example graph contains a cycle, so it uses the general iterative schedule
rather than the specialized tree schedule.

---

## Sum-Product Moment Inference

The moment form runs sum-product Gaussian belief propagation and stores each
message and marginal as a mean vector and covariance matrix. Create a
moment-form inference state with [`moment`](@ref), then use [`gbp!`](@ref) to run
iterative updates.

Each iteration updates factor-to-variable messages, variable-to-factor messages,
and the stored variable marginals. The result is a set of Gaussian marginals.

```@example iterative_gaussian_inference
inference = moment(graph; mean = 0.0, covariance = 1e6)
gbp!(graph, inference; iterations = 50, tolerance = 1e-6)

nothing # hide
```

The `mean` and `covariance` keywords define the default initial belief for
variables created without their own initial `mean` or `covariance`. Scalar
defaults are expanded to each variable dimension.

If `iterations` is omitted, [`gbp!`](@ref) runs at most 10 iterations. The
`iterations` keyword always acts as the maximum number of iterations. The
`tolerance` keyword is applied to [`maxMarginalChange`](@ref). If `tolerance` is
omitted, no early stopping is used.

Use [`marginalMean`](@ref) and [`marginalCovariance`](@ref) when only one field is
needed. Use [`marginal`](@ref) when the full marginal object is needed, or
[`marginals`](@ref) to return all stored variable marginals:

```@example iterative_gaussian_inference
marginalMean(graph, inference, :x1)
marginalCovariance(graph, inference, :x1)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to
one variable or factor, and [`printMarginal`](@ref) prints sum-product marginals.

---

## Sum-Product Canonical Inference

The canonical form computes the same sum-product Gaussian marginals as the
moment form, but represents messages with information vectors and precision
matrices. Products of Gaussian messages are represented internally as sums of
information vectors and precision matrices.

Create a canonical-form inference state with [`canonical`](@ref), then use
[`gbp!`](@ref) to run iterative updates:

```@example iterative_gaussian_inference
inference = canonical(graph; mean = 0.0, covariance = 1e6)
gbp!(graph, inference; iterations = 50, tolerance = 1e-6)

nothing # hide
```

Marginals can still be inspected as means and covariances with the same public
lookup functions used for the moment form.

---

## [Min-Sum MAP Inference](@id gaussian-min-sum-belief-propagation)

The min-sum form computes MAP estimates. It is the negative-log quadratic form
of Gaussian max-product belief propagation: factor potentials become
least-squares costs, messages are quadratic costs, and variable beliefs are MAP
estimates.

This form optimizes the most likely state instead of computing marginal
covariances. The min-sum inference state stores [`QuadraticMessage`](@ref)
objects, which represent quadratic costs rather than Gaussian marginals. Use
[`moment`](@ref) or [`canonical`](@ref) when Gaussian marginals are needed.

Create a min-sum inference state with [`minsum`](@ref), then use [`gbp!`](@ref)
to run iterative updates:

```@example iterative_gaussian_inference
inference = minsum(graph; mean = 0.0, covariance = 1e6)
gbp!(graph, inference; iterations = 50, tolerance = 1e-6)

nothing # hide
```

The `mean` and `covariance` keywords define the default initial belief for
variables created without their own initial `mean` or `covariance`. For
[`minsum`](@ref), this belief is converted to an initial quadratic cost message.

Each iteration updates factor-to-variable messages, variable-to-factor messages,
and the stored variable MAP estimates. For [`minsum`](@ref), the result is a set
of MAP estimates.

The `tolerance` keyword is applied to [`maxEstimateChange`](@ref).

Use [`estimate`](@ref) for one variable and [`estimates`](@ref) for all
variables:

```@example iterative_gaussian_inference
estimate(graph, inference, :x1)
estimates(graph, inference)

nothing # hide
```

For interactive inspection, [`printMessages`](@ref) prints messages adjacent to
one variable or factor, and [`printEstimate`](@ref) prints min-sum MAP estimates.

---

## Manual Message Updates

Use the manual message functions when you want to control the iteration loop
yourself. They run the same message-passing computation as [`gbp!`](@ref), but
split it into explicit message and result steps.

[`messages!`](@ref) runs one full message step: factor-to-variable messages and
variable-to-factor messages. It does not recompute variable results by itself, so
call [`marginals!`](@ref) for [`moment`](@ref) and [`canonical`](@ref), or
[`estimates!`](@ref) for [`minsum`](@ref):

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

For directional control, call the factor-to-variable and variable-to-factor
passes explicitly:

```@example iterative_gaussian_inference
inference = canonical(graph)

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

Use [`maxMessageChange`](@ref) to compare message states between two iterations.
For stored variable results, use [`maxMarginalChange`](@ref) with [`moment`](@ref)
or [`canonical`](@ref), and [`maxEstimateChange`](@ref) with [`minsum`](@ref).

```@example iterative_gaussian_inference
inference = canonical(graph)

for _ in 1:100
    oldVariableMessages = deepcopy(inference.variableToFactor)
    oldFactorMessages = deepcopy(inference.factorToVariable)
    oldMarginals = deepcopy(inference.marginal)

    messages!(graph, inference)
    marginals!(graph, inference)

    messageDiff = maxMessageChange(graph, inference, oldVariableMessages, oldFactorMessages)
    marginalDiff = maxMarginalChange(graph, inference, oldMarginals)

    if messageDiff < 1e-6 && marginalDiff < 1e-6
        break
    end
end

nothing # hide
```

The message change includes both factor-to-variable and variable-to-factor
messages. Use [`maxFactorMessageChange`](@ref) or
[`maxVariableMessageChange`](@ref) to inspect one direction.

---

## Sequential Schedule

Sequential scheduling is the default update order for [`gbp!`](@ref) and
[`messages!`](@ref). Factor-to-variable messages are updated first, followed by
variable-to-factor messages. Messages are updated in place, so later updates in
the same pass can use values computed earlier in that pass.

You can request sequential scheduling explicitly with `schedule = :sequential`:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 20, schedule = :sequential)

nothing # hide
```

For a manual loop, call [`messages!`](@ref) directly. Sequential scheduling is
the default:

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

Flooding scheduling updates messages in a next-message buffer before committing
them to the inference state. With `schedule = :flooding`, updates within a pass
read from the old message arrays and are committed together.

Use [`gbp!`](@ref) with `schedule = :flooding` to run iterative inference with
this update order:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 20, schedule = :flooding)

nothing # hide
```

For a manual flooding loop, pass the same schedule keyword to [`messages!`](@ref):

```@example iterative_gaussian_inference
inference = moment(graph)

for _ in 1:20
    messages!(graph, inference; schedule = :flooding)
    marginals!(graph, inference)
end

nothing # hide
```

Flooding scheduling works with [`moment`](@ref), [`canonical`](@ref), and
[`minsum`](@ref) inference states. The keyword form is the normal public entry
point. Create a [`FloodingSchedule`](@ref) explicitly only for low-level
[`messages!`](@ref) calls that need a schedule object:

```@example iterative_gaussian_inference
inference = moment(graph)
schedule = floodingSchedule(graph, inference)

for _ in 1:20
    messages!(graph, inference, schedule)
    marginals!(graph, inference)
end

nothing # hide
```

The returned schedule is mutable step state for the selected graph and inference
object. Recreate it after changing graph topology.

There are no separate public flooding calls for the two message directions,
because committing one direction before the other would change the flooding
semantics.

---

## Residual Scheduling

Residual scheduling recomputes directed messages and updates only those whose
new values differ most from their current values. This can reduce work on loopy
graphs where a small subset of messages dominates the current error.

Use [`gbp!`](@ref) with `schedule = :residual` to run a complete
residual-scheduled loop:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 50, schedule = :residual, updateFraction = 0.2)

nothing # hide
```

The `updateFraction` keyword controls the fraction of directed messages updated
in each residual step. Use `updateCount` when you want a fixed number of
largest-residual messages:

```@example iterative_gaussian_inference
inference = moment(graph)
gbp!(graph, inference; iterations = 50, schedule = :residual, updateCount = 5)

nothing # hide
```

For a manual residual loop, pass the same residual schedule keywords to
[`messages!`](@ref):

```@example iterative_gaussian_inference
inference = moment(graph)

for _ in 1:20
    messages!(graph, inference; schedule = :residual, updateFraction = 0.1)
    marginals!(graph, inference)
end

nothing # hide
```

Residual scheduling works with [`moment`](@ref), [`canonical`](@ref), and
[`minsum`](@ref) inference states. The keyword form is the normal public entry
point. Create a [`ResidualSchedule`](@ref) explicitly when you want to inspect or
reuse the selected low-level schedule object:

```@example iterative_gaussian_inference
inference = moment(graph)
schedule = residualSchedule(graph, inference; updateFraction = 0.1)

for _ in 1:20
    messages!(graph, inference, schedule)
    marginals!(graph, inference)
end

nothing # hide
```

The residual schedule object stores mutable step state for the selected graph
and inference object. Recreate it after changing graph topology or replacing the
inference object.

There are no separate public residual calls for complete factor-to-variable or
variable-to-factor passes, because each residual step selects individual
directed messages.

---

## Broadcast Message Computation

The `broadcast` keyword selects grouped message computation for [`moment`](@ref),
[`canonical`](@ref), and [`minsum`](@ref). Pass it to [`gbp!`](@ref) for a full
inference run:

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

Grouped computation first accumulates all incoming contributions for a factor or
variable. Each outgoing message is then formed by removing the contribution from
the target edge.

This computation mode can be combined with sequential or flooding schedules. It
is not used with residual scheduling, because residual scheduling selects
individual directed messages by residual rather than running complete grouped
message passes.

---

## Damping

Damping blends a new factor-to-variable message with the previous message. It is
often useful on loopy graphs where undamped message updates may oscillate.

Enable global damping through [`gbp!`](@ref):

```@example iterative_gaussian_inference
inference = minsum(graph)
gbp!(graph, inference; iterations = 50, damping = true, prob = 1.0, alpha = 0.35)

nothing # hide
```

The `alpha` keyword is the weight of the previous message. The `prob` keyword is
the probability that damping is applied to an eligible update.

Damping uses the same `damping`, `prob`, and `alpha` keywords with sequential,
flooding, and residual schedules. It is applied to factor-to-variable message
updates.

The same damping keywords can be used with the manual [`messages!`](@ref) API:

```@example iterative_gaussian_inference
inference = canonical(graph)

for _ in 1:20
    messages!(graph, inference; damping = true, prob = 1.0, alpha = 0.35)
end

marginals!(graph, inference)

nothing # hide
```

When running the two message directions manually, pass the damping keywords to
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

Use [`undampEdges!`](@ref) to remove damping from selected edges:

```@example iterative_gaussian_inference
undampEdges!(graph, inference; variable = :x1, factor = "f3")

nothing # hide
```

---

## Freezing Updates

Freezing keeps selected messages unchanged while the rest of the graph continues
updating.

To freeze all outgoing messages from a factor or variable, use:

```@example iterative_gaussian_inference
freezeFactor!(graph, inference, "f3")
freezeVariable!(graph, inference, :x2)

nothing # hide
```

To freeze both message directions on one edge, use:

```@example iterative_gaussian_inference
freezeEdge!(graph, inference; variable = :x1, factor = "f3")

nothing # hide
```

Resume updates by unfreezing the same selections:

```@example iterative_gaussian_inference
unfreezeFactor!(graph, inference, "f3")
unfreezeVariable!(graph, inference, :x2)
unfreezeEdge!(graph, inference; variable = :x1, factor = "f3")

nothing # hide
```

Inspection helpers return booleans and use the same references:

```julia
isFrozenFactor(graph, inference, "f3")
isFrozenVariable(graph, inference, :x2)
isFrozenEdge(graph, inference; variable = :x1, factor = "f3")
```

---

## Dynamic Graph Changes

Use [`addVariable!`](@ref) with an inference object when a new variable node must
be added after inference has started. A new variable is usually followed by at
least one new factor that connects it to the graph:

```@example iterative_gaussian_inference
addVariable!(graph, inference, :x4, 1; label = "x4")
addFactor!(graph, inference, :x3, :x4, 0.4, [1.0 -1.0], 0.2; label = "f6")

gbp!(graph, inference; iterations = 50)

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
created from the same graph before the topology change no longer matches the
graph and should be recreated or updated separately.

Use [`updateFactor!`](@ref) when an existing factor's mean, coefficient, or
covariance changes without changing graph topology:

```@example iterative_gaussian_inference
updateFactor!(graph, inference; factor = "f2", mean = [1.8, -0.2], covariance = [0.1, 0.1])

gbp!(graph, inference; iterations = 50)

nothing # hide
```

This reuses the existing inference state as a warm start. The message state still
contains information from previous iterations, which is useful for dynamic
inference. For an independent solve, update the graph and then create a fresh
[`moment`](@ref), [`canonical`](@ref), or [`minsum`](@ref) object.

Because [`updateFactor!`](@ref) does not change graph topology, other inference
objects created from the same graph still structurally match the graph. However,
their messages have not been refreshed for the new factor parameters.