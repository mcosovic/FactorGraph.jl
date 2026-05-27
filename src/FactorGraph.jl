module FactorGraph

using LinearAlgebra, Random

# Factor Graph
include("graph/core.jl")
include("graph/gaussian.jl")
include("graph/discrete.jl")

# Inference
include("inference/core.jl")
include("inference/gaussianMoment.jl")
include("inference/gaussianCanonical.jl")
include("inference/discreteSumProduct.jl")
include("inference/discreteMinSum.jl")
include("inference/gaussianMinSum.jl")
include("schedule/forwardBackward.jl")
include("inference/tree.jl")

# Schedule
include("schedule/core.jl")
include("schedule/sequential.jl")
include("schedule/flooding.jl")
include("schedule/residual.jl")

# Utility
include("utility/diagnostics.jl")
include("utility/wls.jl")
include("utility/print.jl")

export AbstractFactorGraph,
    GaussianFactorGraph,
    DiscreteFactorGraph,
    TreeFactorGraph

export VariableId,
    VariableRef,
    FactorRef,
    StateRef,
    ComponentRef

export GaussianVariable,
    GaussianFactor,
    DiscreteVariable,
    DiscreteFactor,
    Edge

export factorGraph,
    treeFactorGraph,
    addVariable!,
    addFactor!,
    updateFactor!

export variableIndex,
    variableDimension,
    factorIndex,
    edgeIndex,
    edgeIndices

export componentIndex,
    componentValue,
    coefficientBlock,
    coefficientBlocks

export stateIndex,
    stateValue,
    factorAxis,
    factorAxes

export AbstractInference,
    GaussianInference,
    GaussianSumProductInference,
    DiscreteInference,
    AbstractSumProductInference,
    AbstractMinSumInference,
    GaussianMomentInference,
    GaussianCanonicalInference,
    GaussianMinSumInference,
    DiscreteSumProductInference,
    DiscreteMinSumInference,
    QuadraticMessage

export moment,
    canonical,
    minsum,
    sumproduct

export factorToVariableMessages!,
    variableToFactorMessages!,
    messages!,
    marginals!,
    gbp!

export marginals,
    marginal,
    marginalMean,
    marginalCovariance,
    marginalProbability,
    estimates,
    estimate,
    estimates!,
    maxEstimateChange

export ForwardBackwardSchedule,
    SequentialSchedule,
    FloodingSchedule,
    ResidualSchedule

export forwardBackwardSchedule,
    sequentialSchedule,
    floodingSchedule,
    residualSchedule

export reset!,
    refresh!,
    forwardStep!,
    backwardStep!,
    forward!,
    backward!,
    forwardBackward!,
    residualStep!

export isFrozenFactor,
    isFrozenVariable,
    freezeFactor!,
    freezeVariable!,
    unfreezeFactor!,
    unfreezeVariable!,
    isFrozenEdge,
    freezeEdge!,
    unfreezeEdge!

export areDampedEdges,
    dampEdges!,
    undampEdges!

export WeightedLeastSquaresResult,
    solveWLS

export residuals,
    normalizedResiduals,
    compareMeanWithWLS,
    maxMeanError,
    compareCovarianceWithWLS,
    maxCovarianceError

export maxVariableMessageChange,
    maxFactorMessageChange,
    maxMessageChange,
    maxMarginalChange

export printModel,
    printGraph,
    printEdges,
    printMessages,
    printMarginal,
    printEstimate,
    printWLS
end
