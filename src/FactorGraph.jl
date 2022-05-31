module FactorGraph

using SparseArrays, LinearAlgebra
using Random

### Form a continuous factor graph and initialize messages and marginals
include("continuousGraph.jl")
export continuousModel, continuousTreeModel

### Form a discrete factor graph and initialize messages and marginals
include("discreteGraph.jl")
export discreteTreeModel

### Factor graph manipulation
include("graphManipulation.jl")
export freezeFactor!, defreezeFactor!, freezeVariable!, defreezeVariable!,
       freezeVariableFactor!, defreezeVariableFactor!, freezeFactorVariable!, defreezeFactorVariable!,
       hideFactor!, addFactors!

### Inference
include("inference.jl")
export marginal, marginalUnnormalized, dynamicFactor!, ageingVariance!, damping!

### Vanilla GBP algorithm
include("vanillaGBP.jl")
export messageFactorVariable, meanFactorVariable, varianceFactorVariable,
       messageDampFactorVariable, meanDampFactorVariable,
       messageVariableFactor, meanVariableFactor, varianceVariableFactor

### Broadcast GBP algorithm
include("broadcastGBP.jl")
export messageFactorVariableBroadcast, meanFactorVariableBroadcast, varianceFactorVariableBroadcast,
       messageDampFactorVariableBroadcast, meanDampFactorVariableBroadcast,
       messageVariableFactorBroadcast, meanVariableFactorBroadcast, varianceVariableFactorBroadcast

### Broadcast GBP with Kahan-Babuska algorithm
include("kahanGBP.jl")
export messageFactorVariableKahan, meanFactorVariableKahan, varianceFactorVariableKahan,
       messageDampFactorVariableKahan, meanDampFactorVariableKahan,
       messageVariableFactorKahan, meanVariableFactorKahan, varianceVariableFactorKahan

# Tree factor graph
include("treeGraph.jl")
export forwardVariableFactor, forwardFactorVariable, backwardVariableFactor, backwardFactorVariable, isTree

# Compute and show results
include("utility.jl")
export wls, errorMetric

end # FactorGraph
