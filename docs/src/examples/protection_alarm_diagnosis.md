# [Protection Alarm Diagnosis](@id protection-alarm-diagnosis-example)

```@meta
CurrentModule = FactorGraph
```

This example builds a small protection-alarm diagnosis model. A control center
receives relay, breaker, and voltage indications from a feeder section. The
variables are finite-state health or indication states, and the factors encode
domain likelihoods. The graph contains cycles, so iterative discrete
belief propagation is the natural schedule.

---

## Variable Nodes

The latent fault location, breaker state, relay indications, and voltage alarm
are modeled as discrete variables:

```@example protection_alarm_diagnosis
using FactorGraph

variables = [
    DiscreteVariable(:fault, 3; label = "fault", states = [:none, :line12, :line23]),
    DiscreteVariable(:breaker, 2; label = "breaker", states = [:closed, :open]),
    DiscreteVariable(:relay12, 2; label = "relay12", states = [:quiet, :trip]),
    DiscreteVariable(:relay23, 2; label = "relay23", states = [:quiet, :trip]),
    DiscreteVariable(:voltage, 2; label = "voltage", states = [:normal, :low]),
]

nothing # hide
```

---

## Prior and Sensor Model

The fault prior is a unary factor. Pairwise factors encode how likely each
relay or voltage alarm is under each fault state:

```@example protection_alarm_diagnosis
f1 = DiscreteFactor(:fault, [0.92, 0.05, 0.03]; label = "prior_fault", initialize = true)
f2 = DiscreteFactor(:fault, :relay12, [0.9 0.1; 0.1 0.9; 0.7 0.2]; label = "fault_relay12")
f3 = DiscreteFactor(:fault, :relay23, [0.9 0.1; 0.7 0.3; 0.1 0.9]; label = "fault_relay23")
f4 = DiscreteFactor(:fault, :voltage, [0.9 0.1; 0.2 0.7; 0.2 0.8]; label = "fault_voltage")

nothing # hide
```

---

## Operational Coupling

The breaker and measured voltage are not independent: an open breaker makes a
low-voltage indication more likely. Relay 12 can also be affected by breaker
operation, which creates a cycle in the factor graph:

```@example protection_alarm_diagnosis
f5 = DiscreteFactor(:breaker, :voltage, [0.9 0.1; 0.2 0.8]; label = "breaker_voltage")
f6 = DiscreteFactor(:breaker, :relay12, [0.8 0.1; 0.3 0.7]; label = "breaker_relay12")

nothing # hide
```

---

## Observed Evidence

Observed SCADA indications are represented as unary likelihood factors. Here
relay 12 trips, relay 23 stays quiet, the breaker is open, and the voltage is
low:

```@example protection_alarm_diagnosis
f7 = DiscreteFactor(:relay12, [0.05, 0.95]; label = "obs_relay12")
f8 = DiscreteFactor(:relay23, [0.95, 0.05]; label = "obs_relay23")
f9 = DiscreteFactor(:breaker, [0.05, 0.95]; label = "obs_breaker")
f10 = DiscreteFactor(:voltage, [0.10, 0.90]; label = "obs_voltage")

nothing # hide
```

---

## Factor Graph Construction

Collect the factors and build the factor graph:

```@example protection_alarm_diagnosis
factors = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]
graph = factorGraph(variables, factors)

nothing # hide
```

The graph can be rendered as an SVG factor graph figure:

```@example protection_alarm_diagnosis
saveGraphFigure("../pad.svg", graph)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../pad.svg"
    type="image/svg+xml"
    aria-label="Protection alarm diagnosis factor graph"
    style="width: 60%; height: auto;">
    <a href="../pad.svg">Protection alarm diagnosis factor graph</a>
  </object>
</div>
```

---

## Running Belief Propagation

Run damped sum-product belief propagation on the graph:

```@example protection_alarm_diagnosis
inference = sumproduct(graph)

gbp!(
    graph, inference;
    iterations = 60, tolerance = 1e-8, schedule = :flooding, damping = true
)

nothing # hide
```

And inspect the posterior probabilities:

```@example protection_alarm_diagnosis
printMarginal(graph, inference; variable = :fault)
printMarginal(graph, inference; variable = :breaker)
```

---

## Alarm Update

If a later relay 23 indication changes to trip, update the evidence factor and
continue from the current messages:

```@example protection_alarm_diagnosis
updateFactor!(graph, inference; factor = "obs_relay23", table = [0.05, 0.95])
gbp!(graph, inference; iterations = 40, tolerance = 1e-8, damping = true)

printMarginal(graph, inference; variable = :fault)
```

The inference object is reused, so the previous messages act as a warm start.
