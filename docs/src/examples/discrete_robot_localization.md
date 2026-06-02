# [Discrete Robot Localization](@id discrete-robot-localization-example)

```@meta
CurrentModule = FactorGraph
```

This example builds a grid-cell robot localization problem. The robot moves
through three corridor cells, and the state variables are the possible cell at
consecutive time steps. Motion and beacon observations form a chain, so the
factor graph is tree-structured and forward-backward discrete belief
propagation is exact.

---

## Variable Nodes

Each position variable has the same three possible cell states:

```@example discrete_robot_localization
using FactorGraph

cells = [:A, :B, :C]

variables = [
    DiscreteVariable(:x1, length(cells); label = "x1", states = cells),
    DiscreteVariable(:x2, length(cells); label = "x2", states = cells),
    DiscreteVariable(:x3, length(cells); label = "x3", states = cells),
    DiscreteVariable(:x4, length(cells); label = "x4", states = cells),
]

nothing # hide
```

---

## Prior Factor Node

The robot is believed to start near cell `A`:

```@example discrete_robot_localization
prior = DiscreteFactor(:x1, [0.85, 0.10, 0.05]; label = "prior_x1", initialize = true)

nothing # hide
```

---

## Motion Factor Nodes

The motion model favors moving one cell to the right, while still allowing the
robot to stay in place or slip:

```@example discrete_robot_localization
motion = [
    0.20 0.75 0.05
    0.05 0.20 0.75
    0.02 0.18 0.80
]

motion12 = DiscreteFactor(:x1, :x2, motion; label = "motion12")
motion23 = DiscreteFactor(:x2, :x3, motion; label = "motion23")
motion34 = DiscreteFactor(:x3, :x4, motion; label = "motion34")

nothing # hide
```

Rows follow the previous cell and columns follow the next cell.

---

## Beacon Factor Nodes

Beacon observations are unary likelihood factors. A beacon near cell `C` is
detected at time step 3, and another observation again favors cell `C` at time
step 4:

```@example discrete_robot_localization
beacon3 = DiscreteFactor(:x3, [0.05, 0.15, 0.80]; label = "beacon3")
beacon4 = DiscreteFactor(:x4, [0.03, 0.17, 0.80]; label = "beacon4")

nothing # hide
```

---

## Factor Graph Construction

Collect the factors and build a tree graph with a selected root:

```@example discrete_robot_localization
factors = [prior, motion12, motion23, motion34, beacon3, beacon4]

graph = factorGraph(variables, factors; root = :x1)

nothing # hide
```

The graph can be rendered as an SVG factor graph figure:

```@example discrete_robot_localization
saveGraphFigure("../drl.svg", graph; layout = (columnSpacing = 120,))

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../drl.svg"
    type="image/svg+xml"
    aria-label="Discrete robot localization factor graph">
    <a href="../drl.svg">Discrete robot localization factor graph</a>
  </object>
</div>
```

---

## Running Belief Propagation

Run discrete belief propagation on the graph:

```@example discrete_robot_localization
inference = sumproduct(graph)

forwardBackward!(graph, inference)

nothing # hide
```

And inspect the posterior cell probabilities:

```@example discrete_robot_localization
printMarginal(graph, inference)
```

---

## Adding a New Time Step

In an online localization problem, the robot later moves again and a new
position state is needed. Add the new variable, connect it with the latest
motion model, add a new beacon observation, and continue from the current
inference state:

```@example discrete_robot_localization
addVariable!(graph, inference, :x5, length(cells); label = "x5", states = cells)
addFactor!(graph, inference, :x4, :x5, motion; label = "motion45")
addFactor!(graph, inference, :x5, [0.02, 0.18, 0.80]; label = "beacon5")

forwardStep!(graph, inference; variable = :x5, factor = "beacon5")
forwardStep!(graph, inference; variable = :x5, factor = "motion45")
backward!(graph, inference)

marginals!(graph, inference)
```

The previous messages and marginals are kept as a warm start. Use the
inference-aware `addVariable!` and `addFactor!` forms for this warm start;
graph-only topology changes require a fresh inference object. A full
[`forwardBackward!`](@ref) sweep can be used instead when all messages in the
updated tree should be recomputed.

```@example discrete_robot_localization
printMarginal(graph, inference; variable = :x5)
```
