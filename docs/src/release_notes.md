# Release Notes

---

## Version 0.2.5 - Unreleased

#### Fixed

- Fixed Gaussian min-sum frozen variable and edge message updates so frozen
  outgoing quadratic messages remain unchanged until explicitly unfrozen.
- Fixed inference-aware discrete factor updates so invalid initializing unary
  updates are rejected before mutating the graph.

#### Tests

- Added numerical correctness coverage for Gaussian schedules, damping,
  freeze/unfreeze message preservation, tree exactness and selected
  step-by-step sweeps, root refreshes, warm-start updates, dynamic stale-state
  rejection, residual partial updates, WLS agreement for correlated vector
  models, and discrete brute-force references.

---

## Version 0.2.4 - 2026-06-05

#### Changed

- Removed the redundant factor node `ID` line from `graphFigure` tooltips.
- Replaced graph figure view `depth` with `hops`, using true bipartite graph
  hop expansion from selected variables and factors.
- Replaced separate unary and multi-factor graph spacing layout options with
  scalar or per-gap `rowSpacing` and `columnSpacing` values.
- Expanded scalar Gaussian factor covariances to isotropic covariance matrices
  matching the factor mean dimension.
- Added `messages!(...; schedule = :residual, updateFraction/updateCount)` as
  the one-step counterpart to residual-scheduled `gbp!`.

#### Fixed

- Fixed residual-scheduled `messages!` for discrete min-sum inference.

#### Documentation

- Added a graph visualization example page covering figure layout, labels,
  focused views, styles, highlights, and tree views.
- Clarified when Gaussian and discrete flooding and residual schedule objects
  are needed instead of the `schedule` keyword.
- Tightened Gaussian graph docstrings for labels, component lookup, factor
  updates, and stale inference object guidance.
- Tightened package docstrings for discrete graph construction, residual
  schedules, inference updates, WLS helpers, and graph visualization.

#### Tests

- Added residual scheduling coverage for stale schedules, discrete min-sum tree
  overloads, and zero/negative iteration handling.

---

## Version 0.2.3 - 2026-06-02

#### Added

- Added optional SVG hover tooltips to `graphFigure` output for variables,
  factors, and edges, with configurable summary and full detail levels.
- Added depth-limited graph figure views for drawing selected variables and
  factors with nearby factors and variables.
- Documented graph figure tooltip controls in the Gaussian and discrete factor
  graph guides.

#### Changed

- Embedded generated graph SVGs with `object` so SVG hover tooltips remain
  available in generated documentation pages.

---

## Version 0.2.2 - 2026-06-01

#### Fixed

- Fixed generated SVG graph figures so Documenter includes them in deployed HTML
  documentation.

---

## Version 0.2.1 - 2026-06-01

#### Added

- Added `graphFigure` and `saveGraphFigure` for SVG factor
  graph visualization.

#### Changed

- Changed automatically assigned graph labels to compact forms such as `x1` and
  `f1` instead of `x_1` and `f_1`.
- Renamed `isDampedEdge` to `areDampedEdges` to reflect that the selector can
  match one or more edges.

### Documentation

- Added rendered factor graph figures to the examples.
- Added quick graph drawing sections to the Gaussian and discrete factor graph
  guides.
- Documented damping selector rules for `areDampedEdges`, `dampEdges!`, and
  `undampEdges!`.

---

## Version 0.2.0 - 2026-05-27

#### Added

- Added discrete finite-state factor graph construction and inference utilities.
- Added iterative sum-product and min-sum belief propagation for discrete models.
- Added tree-oriented factor graph views and exact forward-backward inference on
  tree-structured graphs.
- Added flooding, sequential, and residual message scheduling utilities.
- Added dynamic graph update support with stale-inference checks.
- Added weighted least-squares validation, diagnostics, and printing helpers.

#### Changed

- Reorganized the package source into graph, inference, schedule, and utility
  components.
- Expanded Gaussian factor graph support for scalar and vector variables, linear
  Gaussian factors, and Gaussian belief propagation in moment, canonical, and
  min-sum form.
- Reworked and expanded the documentation, including Gaussian models, discrete
  models, examples, schedules, validation, and API references.
- Expanded the test suite across graph construction, inference, schedules, tree
  inference, validation, diagnostics, and printing.
