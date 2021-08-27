module GaussBP

using SparseArrays, LinearAlgebra
using HDF5, XLSX
using Random
using PrettyTables
using Printf

### Form a factor graph and initialize messages and marginals
include("gbp.jl")
export graphicalModel, dynamicInference!, ageingInference!, damping!,
       freezeFactor!, defreezeFactor!, freezeVariable!, defreezeVariable!

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
include("results.jl")
export marginal, wls, errorMetric, displayData

end # GaussBP