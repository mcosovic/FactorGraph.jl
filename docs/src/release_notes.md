# Release Notes

## Version 0.2.1 - Unreleased

#### Changed

- Renamed `isDampedEdge` to `areDampedEdges` to reflect that the selector can
  match one or more edges.

### Documentation

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