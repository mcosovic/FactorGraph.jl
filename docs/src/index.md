# FactorGraph

```@meta
CurrentModule = FactorGraph
```

FactorGraph is a Julia package for constructing factor graphs and running
belief propagation algorithms.

The package supports Gaussian factor graphs with linear Gaussian factors and
discrete factor graphs with finite-state variables and factor tables. It
provides sum-product and min-sum inference, Gaussian belief propagation in
moment and canonical forms, several message schedules, forward-backward
inference for tree-structured factor graphs, and dynamic graph updates.

---

### Documentation

The documentation is organized around the two model families supported by the
package. The Gaussian Models section covers linear Gaussian factor graphs,
iterative Gaussian belief propagation, min-sum inference, and forward-backward
inference on trees. The Discrete Models section covers finite-state factor
graphs, sum-product inference, min-sum inference, and tree inference.

The Examples section shows complete workflows built from these APIs, including
how to construct, update, visualize, and run inference on factor graphs. The API
section collects the public types, graph construction functions, inference
objects, validation utilities, visualization utilities, and printing helpers.

---

### Installation

FactorGraph is available through the Julia package manager. Install it with:

```julia-repl
pkg> add FactorGraph
```

Then load it with:

```julia
using FactorGraph
```

---

### Gaussian Quick Start

```@example gaussian_quick_start
using FactorGraph

variables = [
    GaussianVariable(:x1, 1),
    GaussianVariable(:x2, 1)
]

factors = [
    GaussianFactor(:x1, 1.0, 1.0, 0.1),
    GaussianFactor(:x2, 2.0, 1.0, 0.1),
    GaussianFactor(:x1, :x2, -1.0, [1.0 -1.0], 0.2)
]

graph = factorGraph(variables, factors)
inference = moment(graph)

gbp!(graph, inference; iterations = 30)

nothing # hide
```

---

### Discrete Quick Start

```@example discrete_quick_start
using FactorGraph

variables = [
    DiscreteVariable(:x1, 2; states = [:off, :on]),
    DiscreteVariable(:x2, 2; states = [:low, :high])
]

factors = [
    DiscreteFactor(:x1, [0.8, 0.2]),
    DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9])
]

graph = factorGraph(variables, factors)
inference = sumproduct(graph)

gbp!(graph, inference; iterations = 10)

nothing # hide
```