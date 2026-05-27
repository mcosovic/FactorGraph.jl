# [DC State Estimation](@id dc-state-estimation-example)

```@meta
CurrentModule = FactorGraph
```

This example builds a three-bus DC power system state estimation problem, where all susceptances
of the branches are equal 10 per-unit. The buses are connected in a triangle, so the factor graph
contains a cycle and is a natural fit for iterative Gaussian belief propagation.

---

## Variable Nodes

The state variables are the bus voltage angles at all three buses. Each angle is modeled as a
scalar variable:

```@example dc_state_estimation
using FactorGraph

variables = [
    GaussianVariable(:θ₁, 1; label = "Bus 1"),
    GaussianVariable(:θ₂, 1; label = "Bus 2"),
    GaussianVariable(:θ₃, 1; label = "Bus 3")
]

nothing # hide
```

---

## Slack Factor Node

Bus 1 is the slack bus. The slack angle measurement is represented as a unary
factor:

```@example dc_state_estimation
slack = GaussianFactor(:θ₁, 0.0, 1.0, 1e-6; label = "Slack")

nothing # hide
```

---

## Flow Factor Nodes

Branch-flow measurements connect the two endpoint angle variables. Here we use
two flow active measurements between buses 1 and 2 and buses 3 and 1.

```@example dc_state_estimation
P12 = GaussianFactor(:θ₁, :θ₂, 0.70, [10.0 -10.0], 0.0016; label = "P12")
P31 = GaussianFactor(:θ₃, :θ₁, 0.30, [10.0 -10.0], 0.0015; label = "P31")

nothing # hide
```

---

## Injection Factor Node

An active-power injection measurement at bus 2 depends on all neighboring bus
angles:

```@example dc_state_estimation
P2 = GaussianFactor(:θ₁, :θ₂, :θ₃, -1.70, [-10.0 20.0 -10.0], 0.0025; label = "P2")

nothing # hide
```

---

## Running Belief Propagation

Collect the factors, build the graph, and run Gaussian belief propagation:

```@example dc_state_estimation
factors = [slack, P12, P31, P2]

graph = factorGraph(variables, factors)
inference = canonical(graph)

gbp!(graph, inference; iterations = 60)
```

And print results:

```@example dc_state_estimation
printMarginal(graph, inference)
```

---

## Dynamic Measurement Update

If a measurement changes, update the corresponding factor and continue Gaussian belief propagation
from the current message state. For example, suppose the `P12` flow
measurement changes:

```@example dc_state_estimation
updateFactor!(graph, "P12"; mean = 0.74)
gbp!(graph, inference; iterations = 20)

printMarginal(graph, inference)
```

The same `inference` object is reused, so the messages from the previous run
act as a warm start.

---

## Validation
Compare the loopy Gaussian belief propagation result with the centralized weighted least-squares
solution:

```@example dc_state_estimation
reference = solveWLS(graph)
maxMeanError(graph, inference, reference)
```

Residuals compare each measurement factor with the estimate implied by the
current marginal means. Normalized residuals scale them by the measurement
standard deviation:

```@example dc_state_estimation
residuals(graph, inference)
normalizedResiduals(graph, inference)
```
