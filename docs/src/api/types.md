# [Types](@id types-api)

```@meta
CurrentModule = FactorGraph
```

This page lists the public types used by the factor graph, inference, and scheduling APIs.

The following hierarchy gives a quick overview of the main type relationships.

Factor graph types separate graph models from tree views over those models:

```text
AbstractFactorGraph
├─ GaussianFactorGraph
├─ DiscreteFactorGraph
└─ TreeFactorGraph{<:AbstractFactorGraph}
```

Inference types are grouped first by model family and then by message-passing objective:

```text
AbstractInference
├─ GaussianInference
│  ├─ GaussianSumProductInference
│  │  ├─ GaussianMomentInference
│  │  └─ GaussianCanonicalInference
│  └─ GaussianMinSumInference
└─ DiscreteInference
   ├─ DiscreteSumProductInference
   └─ DiscreteMinSumInference
```

Cross-model aliases collect equivalent inference objectives across model families:

```text
Cross-model inference aliases
├─ AbstractSumProductInference
│  ├─ GaussianSumProductInference
│  └─ DiscreteSumProductInference
└─ AbstractMinSumInference
   ├─ GaussianMinSumInference
   └─ DiscreteMinSumInference
```

---

## Factor Graphs

```@docs
AbstractFactorGraph
TreeFactorGraph
Edge
```

---

### Gaussian Models

```@docs
GaussianVariable
GaussianFactor
GaussianFactorGraph
```

---

### Discrete Models

```@docs
DiscreteVariable
DiscreteFactor
DiscreteFactorGraph
```

---

## Inference

```@docs
AbstractInference
AbstractSumProductInference
AbstractMinSumInference
```

---

### Gaussian Models

```@docs
GaussianInference
GaussianSumProductInference
GaussianMomentInference
GaussianCanonicalInference
GaussianMinSumInference
QuadraticMessage
```

---

### Discrete Models

```@docs
DiscreteInference
DiscreteSumProductInference
DiscreteMinSumInference
```

---

## Schedule

```@docs
SequentialSchedule
FloodingSchedule
ResidualSchedule
ForwardBackwardSchedule
```

---

## References

```@docs
VariableId
VariableRef
FactorRef
StateRef
ComponentRef
```

---

## Solvers

```@docs
WeightedLeastSquaresResult
```
