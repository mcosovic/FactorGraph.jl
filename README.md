# FactorGraph

[![Build][build-badge]][build]
[![Documentation][documentation-badge]][documentation]

<a href="https://mcosovic.github.io/FactorGraph.jl/stable/">
  <img align="left" width="180" src="/docs/src/assets/logo2.png" />
</a>

FactorGraph is a Julia package for constructing factor graphs and running
belief propagation algorithms.

Gaussian models can use scalar or vector variables, linear Gaussian factors,
belief propagation in moment, canonical, and min-sum form, and weighted
least-squares validation utilities.

Discrete finite-state models are represented with variables over finite state
spaces and factor tables, with iterative sum-product inference for marginals
and min-sum inference for MAP estimates.

The Gaussian and discrete APIs share support for sequential, flooding, and
residual schedules, damping, exact forward-backward inference on
tree-structured graphs, dynamic graph updates, and inspection helpers for
graphs, messages, marginals, and estimates.

<br clear="left"/>

## Installation

FactorGraph is registered as `FactorGraph.jl`. Install it from the Julia package
manager:

```julia-repl
pkg> add FactorGraph
```

Then load it in Julia:

```julia
using FactorGraph
```

## Quick Start

```julia
using FactorGraph

variables = [
    GaussianVariable(:x1, 1),
    GaussianVariable(:x2, 1),
]

factors = [
    GaussianFactor(:x1, 1.0, 1.0, 0.1; label = "prior_x1"),
    GaussianFactor(:x2, 2.0, 1.0, 0.1; label = "prior_x2"),
    GaussianFactor(:x1, :x2, -1.0, [1.0 -1.0], 0.2; label = "link"),
]

graph = factorGraph(variables, factors)
inference = moment(graph)

gbp!(graph, inference; iterations = 30)

marginalMean(graph, inference, :x1)
marginalCovariance(graph, inference, :x1)
```

## Tree Inference

When the graph is tree-structured, pass a root variable to construct a
`TreeFactorGraph` and run exact forward-backward belief propagation:

```julia
tree = treeFactorGraph(graph; root = :x1)
inference = moment(graph)

forwardBackward!(tree, inference)
```

## Discrete Models

```julia
variables = [
    DiscreteVariable(:x1, 2; states = [:off, :on]),
    DiscreteVariable(:x2, 2; states = [:low, :high]),
]

factors = [
    DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1"),
    DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "link"),
]

graph = factorGraph(variables, factors)
inference = sumproduct(graph)

gbp!(graph, inference; iterations = 10)

marginalProbability(graph, inference, :x1, :on)
```

## Documentation

The full documentation is available at
[mcosovic.github.io/FactorGraph.jl/stable][documentation].

It includes:

- Gaussian and discrete factor graph construction
- Iterative sum-product and min-sum belief propagation
- Forward-backward inference on tree-structured graphs
- Dynamic graph updates and stale-inference checks
- API references for types, graphs, inference, validation, and printing helpers

## License

FactorGraph is released under the MIT License.

[build-badge]: https://github.com/mcosovic/FactorGraph.jl/actions/workflows/Build.yml/badge.svg
[build]: https://github.com/mcosovic/FactorGraph.jl/actions/workflows/Build.yml
[documentation-badge]: https://github.com/mcosovic/FactorGraph.jl/actions/workflows/Documentation.yml/badge.svg
[documentation]: https://mcosovic.github.io/FactorGraph.jl/stable/
