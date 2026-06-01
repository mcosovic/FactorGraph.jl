# [PMU-Based State Estimation](@id pmu-state-estimation-example)

```@meta
CurrentModule = FactorGraph
```

This example builds a PMU-based state estimation problem in rectangular
coordinates. Each bus voltage is represented by its real and imaginary
components, and each branch has susceptance 10 per-unit. The measured voltage
and current phasors are modeled as Gaussian factors.

---

## Variable Nodes

The state variables are bus voltage phasors in rectangular coordinates. Each
bus voltage stores its real and imaginary components, so each bus is modeled as
a two-dimensional variable:

```@example pmu_state_estimation
using FactorGraph

variables = [
    GaussianVariable(:V₁, 2; label = "Bus 1"),
    GaussianVariable(:V₂, 2; label = "Bus 2"),
    GaussianVariable(:V₃, 2; label = "Bus 3")
]

nothing # hide
```

---

## Voltage Factor Nodes

A PMU voltage measurement is a unary vector-valued factor. For example, a
voltage phasor measurement at bus 1 is represented by a factor connected only
to bus 1. The factor mean stores the measured real and imaginary voltage
components, and the identity coefficient matrix maps the variable directly to
that measurement:

```@example pmu_state_estimation
V1 = GaussianFactor(:V₁, [1.02, 0.00], [1.0 0.0; 0.0 1.0], [1e-4, 1e-4])

nothing # hide
```

We also include a voltage measurement at bus 2:

```@example pmu_state_estimation
V2 = GaussianFactor(:V₂, [1.01, -0.04], [1.0 0.0; 0.0 1.0], [1e-4, 1e-4])

nothing # hide
```

---

## Current Factor Nodes

For a current phasor measurement, the factor connects the voltage variables at
the two branch endpoints. The factor mean stores the measured real and
imaginary current components. The coefficient matrix stores the linear relation
between endpoint voltages and branch current. Here each branch has susceptance
10 per-unit.

The current measurement on branch 1-2:

```@example pmu_state_estimation
I12 = GaussianFactor(:V₁, :V₂, [0.4, -0.1], [0 10.0 0 -10.0; -10.0 0 10.0 0], [4e-4, 4e-4])

nothing # hide
```

The current measurement on branch 3-1 closes the cycle:

```@example pmu_state_estimation
I31 = GaussianFactor(:V₃, :V₁, [0.2, 0.1], [0 10.0 0 -10.0; -10.0 0 10.0 0], [4e-4, 4e-4])

nothing # hide
```

---

## Factor Graph Construction

Collect the factors and build the factor graph:

```@example pmu_state_estimation
factors = [V1, V2, I12, I31]

graph = factorGraph(variables, factors)

nothing # hide
```

The graph can be rendered as an SVG factor graph figure:

```@example pmu_state_estimation
saveGraphFigure("pmuse.svg", graph)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <img src="../pmuse.svg" alt="PMU state estimation factor graph" style="width: 45%; height: auto;">
</div>
```

---

## Running Belief Propagation

Run Gaussian belief propagation on the graph:

```@example pmu_state_estimation
inference = canonical(graph)

gbp!(graph, inference; iterations = 30)
```

And print results:

```@example pmu_state_estimation
printMarginal(graph, inference)
```

---

## Validation

Compare the loopy Gaussian belief propagation result with the centralized weighted least-squares
solution:

```@example pmu_state_estimation
reference = solveWLS(graph)
maxMeanError(graph, inference, reference)
```
