# [Inference API](@id inference-api)

```@meta
CurrentModule = FactorGraph
```

This page describes the public API for constructing inference states, updating
matching graph/inference pairs, running message passing, and reading results.

---

## Construction

### Gaussian Models

```@docs
moment
canonical
minsum(::GaussianFactorGraph)
```

---

### Discrete Models

```@docs
sumproduct
minsum(::DiscreteFactorGraph)
```

---

## Graph Updates

### Gaussian Models

```@docs
addVariable!(::GaussianFactorGraph, ::GaussianInference, ::GaussianVariable)
addFactor!(::GaussianFactorGraph, ::GaussianInference, ::GaussianFactor)
updateFactor!(::GaussianFactorGraph, ::GaussianInference, ::FactorRef)
```

---

### Discrete Models

```@docs
addVariable!(::DiscreteFactorGraph, ::DiscreteSumProductInference, ::DiscreteVariable)
addVariable!(::DiscreteFactorGraph, ::DiscreteMinSumInference, ::DiscreteVariable)
addFactor!(::DiscreteFactorGraph, ::DiscreteSumProductInference, ::DiscreteFactor)
addFactor!(::DiscreteFactorGraph, ::DiscreteMinSumInference, ::DiscreteFactor)
updateFactor!(::DiscreteFactorGraph, ::DiscreteSumProductInference, ::FactorRef)
updateFactor!(::DiscreteFactorGraph, ::DiscreteMinSumInference, ::FactorRef)
```

---

## Freezing and Damping

```@docs
freezeFactor!
unfreezeFactor!
isFrozenFactor
freezeVariable!
unfreezeVariable!
isFrozenVariable
freezeEdge!
unfreezeEdge!
isFrozenEdge
dampEdges!
undampEdges!
isDampedEdge
```

---

## Message Passing

### Gaussian Models

```@docs
factorToVariableMessages!(::GaussianFactorGraph, ::GaussianInference)
variableToFactorMessages!(::GaussianFactorGraph, ::GaussianInference)
messages!(::GaussianFactorGraph, ::GaussianInference)
marginals!
estimates!(::GaussianFactorGraph, ::GaussianMinSumInference)
gbp!(::GaussianFactorGraph, ::GaussianInference)
```

---

### Discrete Models

```@docs
factorToVariableMessages!(::DiscreteFactorGraph, ::DiscreteSumProductInference)
factorToVariableMessages!(::DiscreteFactorGraph, ::DiscreteMinSumInference)
variableToFactorMessages!(::DiscreteFactorGraph, ::DiscreteSumProductInference)
variableToFactorMessages!(::DiscreteFactorGraph, ::DiscreteMinSumInference)
messages!(::DiscreteFactorGraph, ::DiscreteSumProductInference)
messages!(::DiscreteFactorGraph, ::DiscreteMinSumInference)
marginals!(::DiscreteFactorGraph, ::DiscreteSumProductInference)
estimates!(::DiscreteFactorGraph, ::DiscreteMinSumInference)
gbp!(::DiscreteFactorGraph, ::DiscreteSumProductInference)
gbp!(::DiscreteFactorGraph, ::DiscreteMinSumInference)
```

---

## Iterative Scheduling

```@docs
sequentialSchedule
floodingSchedule
residualSchedule
residualStep!
messages!(::GaussianFactorGraph, ::GaussianInference, ::ResidualSchedule)
```

---

## Tree Inference

```@docs
forwardBackwardSchedule(::TreeFactorGraph)
reset!(::ForwardBackwardSchedule)
reset!(::TreeFactorGraph)
refresh!(::TreeFactorGraph)
forwardStep!
backwardStep!
forward!
backward!
forwardBackward!
```

---

## Results

```@docs
marginals
marginal
estimates
estimate
```

---

### Gaussian Models

```@docs
marginalMean(::GaussianFactorGraph, ::GaussianSumProductInference, ::VariableRef)
marginalCovariance(::GaussianFactorGraph, ::GaussianSumProductInference, ::VariableRef)
```

---

### Discrete Models

```@docs
marginalProbability
```
