# [Binary Pump Diagnosis](@id binary-pump-diagnosis-example)

```@meta
CurrentModule = FactorGraph
```

This example builds a small binary fault-diagnosis model for an industrial pump
skid. Each variable has two states. Latent variables represent component health,
and observed variables represent binary alarms from the control system. The
graph contains a cycle, so iterative discrete belief propagation is used.

---

## Variable Nodes

Each variable is binary:

```@example binary_pump_diagnosis
using FactorGraph

states = [:ok, :bad]
alarmStates = [:no, :yes]

variables = [
    DiscreteVariable(:power, 2; label = "power", states = states),
    DiscreteVariable(:fuse, 2; label = "fuse", states = states),
    DiscreteVariable(:motor, 2; label = "motor", states = states),
    DiscreteVariable(:vibration, 2; label = "vibration", states = alarmStates),
    DiscreteVariable(:temperature, 2; label = "temperature", states = alarmStates),
    DiscreteVariable(:trip, 2; label = "trip", states = alarmStates)
]

nothing # hide
```

---

## Component Priors

Unary factors encode prior failure rates. The power feed, fuse, and motor are
expected to be healthy, but not guaranteed:

```@example binary_pump_diagnosis
f1 = DiscreteFactor(:power, [0.97, 0.03]; label = "power_prior", initialize = true)
f2 = DiscreteFactor(:fuse, [0.96, 0.04]; label = "fuse_prior", initialize = true)
f3 = DiscreteFactor(:motor, [0.94, 0.06]; label = "motor_prior", initialize = true)

nothing # hide
```

---

## Fault Propagation Factors

Pairwise factors encode simple engineering relationships. A bad power feed makes
a blown fuse more likely, a blown fuse makes motor trouble more likely, and a bad
motor makes vibration and temperature alarms more likely:

```@example binary_pump_diagnosis
f4 = DiscreteFactor(:power, :fuse, [0.98 0.02; 0.30 0.70]; label = "power_fuse")
f5 = DiscreteFactor(:fuse, :motor, [0.96 0.04; 0.25 0.75]; label = "fuse_motor")
f6 = DiscreteFactor(:motor, :vibration, [0.90 0.10; 0.15 0.85]; label = "motor_vibration")
f7 = DiscreteFactor(:motor, :temperature, [0.92 0.08; 0.2 0.8]; label = "motor_temperature")

nothing # hide
```

---

## Protection Logic

The trip indication depends on both vibration and temperature. This ternary
factor keeps all variables binary while still using a multidimensional table.
The first slice is `trip = :no`, and the second slice is `trip = :yes`:

```@example binary_pump_diagnosis
f8 = DiscreteFactor(
    :vibration, :temperature, :trip,
    cat(
        [0.98 0.30; 0.25 0.05],
        [0.02 0.70; 0.75 0.95];
        dims = 3
    );
    label = "trip_logic"
)

nothing # hide
```

---

## Observed Evidence

The operator sees vibration and temperature alarms, and the protection relay has
tripped:

```@example binary_pump_diagnosis
f9 = DiscreteFactor(:vibration, [0.05, 0.95]; label = "obs_vibration")
f10 = DiscreteFactor(:temperature, [0.10, 0.90]; label = "obs_temperature")
f11 = DiscreteFactor(:trip, [0.02, 0.98]; label = "obs_trip")

nothing # hide
```

---

## Factor Graph Construction

Collect the factors and build the factor graph:

```@example binary_pump_diagnosis
factors = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11]
graph = factorGraph(variables, factors)

nothing # hide
```

The graph can be rendered as an SVG factor graph figure:

```@example binary_pump_diagnosis
saveGraphFigure("../bpd.svg", graph; label = (showEdgeIds = true, tooltipDetail = :full))

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../bpd.svg"
    type="image/svg+xml"
    aria-label="Binary pump diagnosis factor graph"
    style="width: 70%; height: auto;">
    <a href="../../bpd.svg">Binary pump diagnosis factor graph</a>
  </object>
</div>
```

---

## Sum-Product Diagnosis

Run damped loopy sum-product belief propagation to estimate posterior fault
probabilities:

```@example binary_pump_diagnosis
inference = sumproduct(graph)

gbp!(graph, inference; iterations = 80, tolerance = 1e-8, damping = true)

nothing # hide
```

Inspect the posterior health estimates:

```@example binary_pump_diagnosis
printMarginal(graph, inference; variable = :power)
printMarginal(graph, inference; variable = :fuse)
printMarginal(graph, inference; variable = :motor)
```

---

## MAP Diagnosis

Use min-sum when only the most likely binary assignment is needed:

```@example binary_pump_diagnosis
mapInference = minsum(graph)
gbp!(graph, mapInference; iterations = 80, tolerance = 0.0, schedule = :residual)

estimates(graph, mapInference)
```

The MAP estimate is a single most likely explanation, while the sum-product
marginals show how much uncertainty remains around each component.