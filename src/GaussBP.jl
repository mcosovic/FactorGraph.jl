module GaussBP

using SparseArrays, LinearAlgebra
using HDF5, XLSX
using Random
using PrettyTables
using Printf

### Form a factor graph and initialize messages and marginals
include("graphicalmodel.jl")
export graphicalModel, damping!

### Factor graph manipulation
include("graphmanipulation.jl")
export freezeFactor!, defreezeFactor!, freezeVariable!, defreezeVariable!,
       freezeVariableFactor!, defreezeVariableFactor!, freezeFactorVariable!, defreezeFactorVariable!,
       hideFactor!, addFactor!

### Inference
include("inference.jl")
export marginal, dynamicFactor!, ageingVariance!

### Vanilla GBP algorithm
include("vanillaGBP.jl")
export messageFactorVariableVanilla, meanFactorVariableVanilla, varianceFactorVariableVanilla,
       messageDampFactorVariableVanilla, meanDampFactorVariableVanilla,
       messageVariableFactorVanilla, meanVariableFactorVanilla, varianceVariableFactorVanilla

### Computation-efficient GBP algorithm
include("efficientGBP.jl")
export messageFactorVariableEfficient, meanFactorVariableEfficient, varianceFactorVariableEfficient,
       messageDampFactorVariableEfficient, meanDampFactorVariableEfficient,
       messageVariableFactorEfficient, meanVariableFactorEfficient, varianceVariableFactorEfficient

### Computation-efficient GBP algorithm with Kahan-Babuska algorithm
include("kahanGBP.jl")
export messageFactorVariableKahan, meanFactorVariableKahan, varianceFactorVariableKahan,
       messageDampFactorVariableKahan, meanDampFactorVariableKahan,
       messageVariableFactorKahan, meanVariableFactorKahan, varianceVariableFactorKahan

# Compute and show results
include("utility.jl")
export wls, errorMetric, displayData

end # GaussBP