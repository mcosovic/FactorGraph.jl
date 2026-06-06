# FactorGraph

[![Build][build-badge]][build]
[![Documentation][documentation-badge]][documentation]
[![Codecov][codecov-badge]][codecov]

<a href="https://mcosovic.github.io/FactorGraph.jl/stable/">
  <img align="left" width="200" src="./docs/src/assets/logo2.png" />
</a>

FactorGraph is a Julia package for constructing factor graphs and running
belief propagation algorithms.

Gaussian models can use scalar or vector variables, linear Gaussian factors,
belief propagation in moment, canonical, and min-sum form.

Discrete finite-state models are represented with variables over finite state
spaces and factor tables, with iterative sum-product inference for marginals
and min-sum inference for MAP estimates.

The Gaussian and discrete APIs share support for sequential, flooding, and
residual schedules, damping, exact forward-backward inference on
tree-structured graphs, dynamic graph updates, and inspection helpers for
graphs, messages, marginals, and estimates.

## Installation

Install FactorGraph from the Julia package manager:

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

## Features

- Gaussian factor graphs with moment, canonical, and min-sum inference
- Discrete finite-state factor graphs with sum-product and min-sum inference
- Sequential, flooding, and residual message schedules
- Exact forward-backward inference on tree-structured graphs
- Dynamic graph updates with stale-inference checks
- Freezing, damping, diagnostics, WLS validation, and graph visualization

## Documentation

The full documentation is available at
[mcosovic.github.io/FactorGraph.jl/stable][documentation].

It includes:

- Gaussian and discrete model guides
- Iterative and tree-structured inference
- Domain examples and graph visualization
- Dynamic graph updates and validation helpers
- API references and release notes

## License

FactorGraph is released under the MIT License.

[build-badge]: https://github.com/mcosovic/FactorGraph.jl/actions/workflows/Build.yml/badge.svg
[build]: https://github.com/mcosovic/FactorGraph.jl/actions/workflows/Build.yml
[documentation-badge]: https://github.com/mcosovic/FactorGraph.jl/actions/workflows/Documentation.yml/badge.svg
[documentation]: https://mcosovic.github.io/FactorGraph.jl/stable/
[codecov-badge]: https://codecov.io/github/mcosovic/FactorGraph.jl/graph/badge.svg?token=CWT8PIUXWB
[codecov]: https://codecov.io/github/mcosovic/FactorGraph.jl
