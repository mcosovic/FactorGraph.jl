# [Factor Graphs API](@id factor-graphs-api)

```@meta
CurrentModule = FactorGraph
```

This page describes the public API for building, updating, and inspecting factor
graph structures.

---

## Construction and Updates

```@docs
factorGraph
treeFactorGraph
```

---

### Gaussian Models

```@docs
addVariable!(::GaussianFactorGraph, ::GaussianVariable)
addFactor!(::GaussianFactorGraph, ::GaussianFactor)
updateFactor!(::GaussianFactorGraph, ::FactorRef)
```

---

### Discrete Models

```@docs
addVariable!(::DiscreteFactorGraph, ::DiscreteVariable)
addFactor!(::DiscreteFactorGraph, ::DiscreteFactor)
updateFactor!(::DiscreteFactorGraph, ::FactorRef)
```

---

## Graph Lookup

```@docs
variableIndex
variableDimension
factorIndex
edgeIndex
edgeIndices
```

---

### Gaussian Models

```@docs
componentIndex
componentValue
coefficientBlock
coefficientBlocks
```

---

### Discrete Models

```@docs
stateIndex
stateValue
factorAxis
factorAxes
```