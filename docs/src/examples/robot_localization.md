# [Robot Localization on a Chain](@id robot-localization-example)

```@meta
CurrentModule = FactorGraph
```

This example builds a one-dimensional robot localization problem. The robot
moves along a line, and the state variables are robot positions at consecutive
time steps. Odometry measurements connect neighboring positions, so the factor
graph is a chain. This tree structure is a natural fit for forward-backward
Gaussian belief propagation.

---

## Variable Nodes

The state variables are robot positions at four time steps. Each position is
modeled as a scalar variable:

```@example robot_localization
using FactorGraph

variables = [
    GaussianVariable(:x1, 1),
    GaussianVariable(:x2, 1),
    GaussianVariable(:x3, 1),
    GaussianVariable(:x4, 1)
]

nothing # hide
```

---

## Prior Factor Node

The initial position is measured or known approximately. This is represented
as a unary prior factor:

```@example robot_localization
prior = GaussianFactor(:x1, 0.0, 1.0, 0.01)

nothing # hide
```

---

## Odometry Factor Nodes

Odometry measurements connect consecutive robot positions. For example, a
measured displacement between time steps 1 and 2 is represented as a factor
connected to the first two position variables:

```@example robot_localization
odom12 = GaussianFactor(:x1, :x2, 1.05, [-1.0 1.0], 0.04)
odom23 = GaussianFactor(:x2, :x3, 0.95, [-1.0 1.0], 0.04)
odom34 = GaussianFactor(:x3, :x4, 1.10, [-1.0 1.0], 0.04)

nothing # hide
```

These odometry factors form a chain, so the graph is tree-structured.

---

## Landmark Factor Node

Suppose the robot also receives a position measurement at the final time step,
for example from a landmark or GPS-like sensor:

```@example robot_localization
landmark = GaussianFactor(:x4, 3.05, 1.0, 0.02)

nothing # hide
```

This factor does not create a cycle because it is unary.

---

## Factor Graph Construction

Collect the factors and build the graph with a selected root for
forward-backward Gaussian belief propagation:

```@example robot_localization
factors = [prior, odom12, odom23, odom34, landmark]

graph = factorGraph(variables, factors; root = :x1)

nothing # hide
```

The graph can be rendered as an SVG factor graph figure:

```@example robot_localization
saveGraphFigure("../rl.svg", graph; layout = (columnSpacing = 120,))

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../rl.svg"
    type="image/svg+xml"
    aria-label="Robot localization factor graph">
    <a href="../../rl.svg">Robot localization factor graph</a>
  </object>
</div>
```

---

## Running Belief Propagation

Run Gaussian belief propagation on the graph:

```@example robot_localization
inference = canonical(graph)

forwardBackward!(graph, inference)

nothing # hide
```

And print results:

```@example robot_localization
printMarginal(graph, inference)
```

---

## Adding a New Time Step

In an online localization problem, the robot later moves again and a new
position state is needed. Add the new variable node, connect it with the latest
odometry measurement, and continue from the current inference state:

```@example robot_localization
addVariable!(graph, inference, :x5, 1)
addFactor!(graph, inference, :x4, :x5, 0.90, [-1.0 1.0], 0.04; label = "odom45")

forwardStep!(graph, inference; variable = :x5, factor = "odom45")
backward!(graph, inference)
marginals!(graph, inference)
```

The previous messages and marginals are kept as a warm start. The new factor
extends the chain, so the updated graph is still a tree. Use the inference-aware
`addVariable!` and `addFactor!` forms for this warm start; graph-only topology
changes require a fresh inference object. The old forward messages are still
valid because the old part of the chain did not change. The only new forward
message comes from the new terminal state through the new odometry factor. After
that, a backward pass sends the updated information from the root side of the
chain toward all later positions, including the new terminal state. A full
[`forwardBackward!`](@ref) sweep can be used instead when all messages in the
updated tree should be recomputed.

```@example robot_localization
printMarginal(graph, inference)
```

---

## Validation

Compare the current Gaussian belief propagation result with the centralized weighted least-squares
solution:

```@example robot_localization
reference = solveWLS(graph)
maxMeanError(graph, inference, reference)
```
