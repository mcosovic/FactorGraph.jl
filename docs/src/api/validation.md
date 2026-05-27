# [Validation API](@id validation-api)

```@meta
CurrentModule = FactorGraph
```

The groups below follow the validation workflow: message-change checks first,
then stored-result changes, weighted least-squares reference solves, residual
inspection, and WLS comparison utilities.

---



## Message Changes

```@docs
maxVariableMessageChange
maxFactorMessageChange
maxMessageChange
```

---

## Result Changes

```@docs
maxMarginalChange
maxEstimateChange
```

---

## Gaussian Residuals

```@docs
residuals
normalizedResiduals
```

---

## Weighted Least Squares Comparisons

```@docs
solveWLS
compareMeanWithWLS
maxMeanError
compareCovarianceWithWLS
maxCovarianceError
```
