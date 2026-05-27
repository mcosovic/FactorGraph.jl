function maximumVectorChange(current::AbstractVector, reference::AbstractVector)
    return maximumElementChange(current, reference, "Compared vectors")
end

function maximumMatrixChange(current::AbstractVector, reference::AbstractVector)
    return maximumElementChange(current, reference, "Compared matrix vectors")
end

function maximumElementChange(
    current::AbstractVector,
    reference::AbstractVector,
    label::String
)
    if length(current) != length(reference)
        error("$label must have the same length.")
    end

    maximumChange = 0.0

    for index in eachindex(current)
        maximumChange = max(maximumChange, norm(current[index] - reference[index]))
    end

    return maximumChange
end

function assertSameGaussianMomentInferenceShape(
    current::GaussianMomentInference,
    reference::GaussianMomentInference
)
    return assertSameInferenceShape(
        current,
        reference,
        (:variableToFactor, :factorToVariable, :marginal),
        "Moment inference objects"
    )
end

function assertSameGaussianCanonicalInferenceShape(
    current::GaussianCanonicalInference,
    reference::GaussianCanonicalInference
)
    return assertSameInferenceShape(
        current,
        reference,
        (:variableToFactor, :factorToVariable, :marginal),
        "Canonical inference objects"
    )
end

function assertSameDiscreteSumProductInferenceShape(
    current::DiscreteSumProductInference,
    reference::DiscreteSumProductInference
)
    return assertSameInferenceShape(
        current,
        reference,
        (:variableToFactor, :factorToVariable, :marginal),
        "DiscreteSumProductInference objects"
    )
end

function assertSameDiscreteMinSumInferenceShape(
    current::DiscreteMinSumInference,
    reference::DiscreteMinSumInference
)
    return assertSameInferenceShape(
        current,
        reference,
        (:variableToFactor, :factorToVariable, :estimate),
        "DiscreteMinSumInference objects"
    )
end

function assertSameInferenceShape(current, reference, fields::Tuple, label::String)
    for field in fields
        if length(getfield(current, field)) != length(getfield(reference, field))
            error("$label do not have the same shape.")
        end
    end

    return nothing
end

function assertSameLength(current::AbstractVector, reference::AbstractVector, label::String)
    if length(current) != length(reference)
        error("$label do not have the same length.")
    end

    return nothing
end

function maximumGaussianMessageChange(
    current::Vector{GaussianMessage},
    reference::Vector{GaussianMessage}
)
    assertSameLength(current, reference, "Gaussian message vectors")

    meanChange = maximumVectorChange(
        [message.mean for message in current],
        [message.mean for message in reference]
    )
    covarianceChange = maximumMatrixChange(
        [message.covariance for message in current],
        [message.covariance for message in reference]
    )

    return max(meanChange, covarianceChange)
end

function maximumCanonicalMessageChange(
    current::Vector{CanonicalMessage},
    reference::Vector{CanonicalMessage}
)
    assertSameLength(current, reference, "Canonical message vectors")

    informationChange = maximumVectorChange(
        [message.information for message in current],
        [message.information for message in reference]
    )
    precisionChange = maximumMatrixChange(
        [message.precision for message in current],
        [message.precision for message in reference]
    )

    return max(informationChange, precisionChange)
end

function maximumCanonicalMarginalChange(
    current::Vector{CanonicalMarginal},
    reference::Vector{CanonicalMarginal}
)
    assertSameLength(current, reference, "Canonical marginal vectors")

    meanChange = maximumVectorChange(
        [marginal.mean for marginal in current],
        [marginal.mean for marginal in reference]
    )
    covarianceChange = maximumMatrixChange(
        [marginal.covariance for marginal in current],
        [marginal.covariance for marginal in reference]
    )

    return max(meanChange, covarianceChange)
end

function maximumQuadraticMessageChange(
    current::Vector{QuadraticMessage},
    reference::Vector{QuadraticMessage}
)
    assertSameLength(current, reference, "Quadratic message vectors")

    matrixChange = maximumMatrixChange(
        [message.J for message in current],
        [message.J for message in reference]
    )
    vectorChange = maximumVectorChange(
        [message.h for message in current],
        [message.h for message in reference]
    )
    scalarChange = 0.0

    for index in eachindex(current)
        scalarChange = max(scalarChange, abs(current[index].c - reference[index].c))
    end

    return max(matrixChange, vectorChange, scalarChange)
end

"""
    maxVariableMessageChange(
        graph::AbstractFactorGraph, current::AbstractInference, reference::Vector
    )

Return the largest variable-to-factor message change from a previous message snapshot.

# Arguments

- `graph`: Factor graph or tree view.
- `current`: Current inference object.
- `reference`: Previous `current.variableToFactor` snapshot.

# Returns

The maximum message change.

# Notes

The snapshot type must match `current`: `Vector{GaussianMessage}` for moment
inference, `Vector{CanonicalMessage}` for canonical inference,
`Vector{QuadraticMessage}` for Gaussian min-sum inference, and
`Vector{Vector{Float64}}` for discrete inference.

# Example

```julia
oldMessages = deepcopy(inference.variableToFactor)
messages!(graph, inference)

maxVariableMessageChange(graph, inference, oldMessages)
```
"""
function maxVariableMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMomentInference,
    reference::Vector{GaussianMessage}
)
    assertGaussianMomentInferenceMatchesGraph(graph, current)

    return maximumGaussianMessageChange(current.variableToFactor, reference)
end

function maxVariableMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianCanonicalInference,
    reference::Vector{CanonicalMessage}
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, current)

    return maximumCanonicalMessageChange(current.variableToFactor, reference)
end

function maxVariableMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMinSumInference,
    reference::Vector{QuadraticMessage}
)
    assertGaussianMinSumInferenceMatchesGraph(graph, current)

    return maximumQuadraticMessageChange(current.variableToFactor, reference)
end

function maxVariableMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteSumProductInference,
    reference::Vector{Vector{Float64}}
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, current)

    return maximumVectorChange(current.variableToFactor, reference)
end

function maxVariableMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteMinSumInference,
    reference::Vector{Vector{Float64}}
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, current)

    return maximumVectorChange(current.variableToFactor, reference)
end

function maxVariableMessageChange(tree::TreeFactorGraph, current, reference)
    return maxVariableMessageChange(tree.graph, current, reference)
end

"""
    maxFactorMessageChange(
        graph::AbstractFactorGraph, current::AbstractInference, reference::Vector
    )

Return the largest factor-to-variable message change from a previous message snapshot.

# Arguments

- `graph`: Factor graph or tree view.
- `current`: Current inference object.
- `reference`: Previous `current.factorToVariable` snapshot.

# Returns

The maximum message change.

# Notes

The snapshot type must match `current`: `Vector{GaussianMessage}` for moment
inference, `Vector{CanonicalMessage}` for canonical inference,
`Vector{QuadraticMessage}` for Gaussian min-sum inference, and
`Vector{Vector{Float64}}` for discrete inference.

# Example

```julia
oldMessages = deepcopy(inference.factorToVariable)
messages!(graph, inference)

maxFactorMessageChange(graph, inference, oldMessages)
```
"""
function maxFactorMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMomentInference,
    reference::Vector{GaussianMessage}
)
    assertGaussianMomentInferenceMatchesGraph(graph, current)

    return maximumGaussianMessageChange(current.factorToVariable, reference)
end

function maxFactorMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianCanonicalInference,
    reference::Vector{CanonicalMessage}
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, current)

    return maximumCanonicalMessageChange(current.factorToVariable, reference)
end

function maxFactorMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMinSumInference,
    reference::Vector{QuadraticMessage}
)
    assertGaussianMinSumInferenceMatchesGraph(graph, current)

    return maximumQuadraticMessageChange(current.factorToVariable, reference)
end

function maxFactorMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteSumProductInference,
    reference::Vector{Vector{Float64}}
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, current)

    return maximumVectorChange(current.factorToVariable, reference)
end

function maxFactorMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteMinSumInference,
    reference::Vector{Vector{Float64}}
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, current)

    return maximumVectorChange(current.factorToVariable, reference)
end

function maxFactorMessageChange(tree::TreeFactorGraph, current, reference)
    return maxFactorMessageChange(tree.graph, current, reference)
end

"""
    maxMessageChange(
        graph::AbstractFactorGraph, current::AbstractInference,
        referenceVariableToFactor::Vector, referenceFactorToVariable::Vector
    )
    maxMessageChange(
        graph::AbstractFactorGraph, current::AbstractInference, reference::AbstractInference
    )

Return the largest message change from previous message snapshots or an inference state.

# Arguments

- `graph`: Factor graph or tree view.
- `current`: Current inference object.
- `referenceVariableToFactor`: Previous variable-to-factor message snapshot.
- `referenceFactorToVariable`: Previous factor-to-variable message snapshot.
- `reference`: Previous inference object.

# Returns

The maximum message change across both message directions.

# Notes

The current inference object must match `graph`. For moment-form inference this compares
message means and covariances. For canonical-form inference this compares information
vectors and precision matrices. Both variable-to-factor and factor-to-variable messages
are included.

Snapshot arguments must match the message storage type used by `current`:
`Vector{GaussianMessage}` for moment inference, `Vector{CanonicalMessage}` for
canonical inference, `Vector{QuadraticMessage}` for Gaussian min-sum inference,
and `Vector{Vector{Float64}}` for discrete inference.

This is useful for custom convergence loops. Prefer passing previous
`variableToFactor` and `factorToVariable` message snapshots to avoid copying an
entire inference object.

# Example

```julia
oldVarMessages = deepcopy(inference.variableToFactor)
oldFacMessages = deepcopy(inference.factorToVariable)
messages!(graph, inference)

maxMessageChange(graph, inference, oldVarMessages, oldFacMessages)
```
"""
function maxMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMomentInference,
    referenceVariableToFactor::Vector{GaussianMessage},
    referenceFactorToVariable::Vector{GaussianMessage}
)
    return max(
        maxVariableMessageChange(graph, current, referenceVariableToFactor),
        maxFactorMessageChange(graph, current, referenceFactorToVariable)
    )
end

function maxMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianCanonicalInference,
    referenceVariableToFactor::Vector{CanonicalMessage},
    referenceFactorToVariable::Vector{CanonicalMessage}
)
    return max(
        maxVariableMessageChange(graph, current, referenceVariableToFactor),
        maxFactorMessageChange(graph, current, referenceFactorToVariable)
    )
end

function maxMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMinSumInference,
    referenceVariableToFactor::Vector{QuadraticMessage},
    referenceFactorToVariable::Vector{QuadraticMessage}
)
    return max(
        maxVariableMessageChange(graph, current, referenceVariableToFactor),
        maxFactorMessageChange(graph, current, referenceFactorToVariable)
    )
end

function maxMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteSumProductInference,
    referenceVariableToFactor::Vector{Vector{Float64}},
    referenceFactorToVariable::Vector{Vector{Float64}}
)
    return max(
        maxVariableMessageChange(graph, current, referenceVariableToFactor),
        maxFactorMessageChange(graph, current, referenceFactorToVariable)
    )
end

function maxMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteMinSumInference,
    referenceVariableToFactor::Vector{Vector{Float64}},
    referenceFactorToVariable::Vector{Vector{Float64}}
)
    return max(
        maxVariableMessageChange(graph, current, referenceVariableToFactor),
        maxFactorMessageChange(graph, current, referenceFactorToVariable)
    )
end

function maxMessageChange(
    tree::TreeFactorGraph,
    current,
    referenceVariableToFactor,
    referenceFactorToVariable
)
    return maxMessageChange(
        tree.graph,
        current,
        referenceVariableToFactor,
        referenceFactorToVariable
    )
end

function maxMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMomentInference,
    reference::GaussianMomentInference
)
    assertGaussianMomentInferenceMatchesGraph(graph, current)
    assertGaussianMomentInferenceMatchesGraph(graph, reference)
    assertSameGaussianMomentInferenceShape(current, reference)

    return maxMessageChange(
        graph,
        current,
        reference.variableToFactor,
        reference.factorToVariable
    )
end

function maxMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianCanonicalInference,
    reference::GaussianCanonicalInference
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, current)
    assertGaussianCanonicalInferenceMatchesGraph(graph, reference)
    assertSameGaussianCanonicalInferenceShape(current, reference)

    return maxMessageChange(
        graph,
        current,
        reference.variableToFactor,
        reference.factorToVariable
    )
end

function maxMessageChange(
    graph::GaussianFactorGraph,
    current::GaussianMinSumInference,
    reference::GaussianMinSumInference
)
    assertGaussianMinSumInferenceMatchesGraph(graph, current)
    assertGaussianMinSumInferenceMatchesGraph(graph, reference)

    if length(current.variableToFactor) != length(reference.variableToFactor) ||
            length(current.factorToVariable) != length(reference.factorToVariable)
        error("GaussianMinSumInference objects do not have the same message shape.")
    end

    return maxMessageChange(
        graph,
        current,
        reference.variableToFactor,
        reference.factorToVariable
    )
end

function maxMessageChange(
    tree::TreeFactorGraph{GaussianFactorGraph},
    current,
    reference
)
    return maxMessageChange(tree.graph, current, reference)
end

function maxMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteSumProductInference,
    reference::DiscreteSumProductInference
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, current)
    assertDiscreteSumProductInferenceMatchesGraph(graph, reference)
    assertSameDiscreteSumProductInferenceShape(current, reference)

    return maxMessageChange(
        graph,
        current,
        reference.variableToFactor,
        reference.factorToVariable
    )
end

function maxMessageChange(
    graph::DiscreteFactorGraph,
    current::DiscreteMinSumInference,
    reference::DiscreteMinSumInference
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, current)
    assertDiscreteMinSumInferenceMatchesGraph(graph, reference)
    assertSameDiscreteMinSumInferenceShape(current, reference)

    return maxMessageChange(
        graph,
        current,
        reference.variableToFactor,
        reference.factorToVariable
    )
end

function maxMessageChange(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    current::DiscreteSumProductInference,
    reference::DiscreteSumProductInference
)
    return maxMessageChange(tree.graph, current, reference)
end

function maxMessageChange(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    current::DiscreteMinSumInference,
    reference::DiscreteMinSumInference
)
    return maxMessageChange(tree.graph, current, reference)
end

"""
    maxMarginalChange(
        graph::AbstractFactorGraph, current::AbstractInference, referenceMarginal::Vector
    )
    maxMarginalChange(
        graph::AbstractFactorGraph, current::AbstractInference, reference::AbstractInference
    )

Return the largest marginal change from previous marginals or an inference state.

# Arguments

- `graph`: Factor graph or tree view.
- `current`: Current inference object.
- `referenceMarginal`: Previous marginal snapshot.
- `reference`: Previous inference object.

# Returns

The maximum marginal change.

# Notes

The current inference object must match `graph`. The comparison uses marginal means and
covariance matrices for moment and canonical inference states, and probability vectors for
sum-product inference.

# Example

```julia
oldMarginals = deepcopy(inference.marginal)
gbp!(graph, inference; iterations = 1)

maxMarginalChange(graph, inference, oldMarginals)
```
"""
function maxMarginalChange(
    graph::GaussianFactorGraph,
    current::GaussianMomentInference,
    reference::Vector{GaussianMessage}
)
    assertGaussianMomentInferenceMatchesGraph(graph, current)

    return maximumGaussianMessageChange(current.marginal, reference)
end

function maxMarginalChange(
    graph::GaussianFactorGraph,
    current::GaussianCanonicalInference,
    reference::Vector{CanonicalMarginal}
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, current)

    return maximumCanonicalMarginalChange(current.marginal, reference)
end

function maxMarginalChange(
    graph::DiscreteFactorGraph,
    current::DiscreteSumProductInference,
    reference::Vector{Vector{Float64}}
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, current)

    return maximumVectorChange(current.marginal, reference)
end

function maxMarginalChange(
    tree::TreeFactorGraph{GaussianFactorGraph},
    current::GaussianMomentInference,
    reference::Vector{GaussianMessage}
)
    return maxMarginalChange(tree.graph, current, reference)
end

function maxMarginalChange(
    tree::TreeFactorGraph{GaussianFactorGraph},
    current::GaussianCanonicalInference,
    reference::Vector{CanonicalMarginal}
)
    return maxMarginalChange(tree.graph, current, reference)
end

function maxMarginalChange(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    current::DiscreteSumProductInference,
    reference::Vector{Vector{Float64}}
)
    return maxMarginalChange(tree.graph, current, reference)
end

function maxMarginalChange(
    graph::GaussianFactorGraph,
    current::GaussianMomentInference,
    reference::GaussianMomentInference
)
    assertGaussianMomentInferenceMatchesGraph(graph, current)
    assertGaussianMomentInferenceMatchesGraph(graph, reference)
    assertSameGaussianMomentInferenceShape(current, reference)

    return maxMarginalChange(graph, current, reference.marginal)
end

function maxMarginalChange(
    graph::GaussianFactorGraph,
    current::GaussianCanonicalInference,
    reference::GaussianCanonicalInference
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, current)
    assertGaussianCanonicalInferenceMatchesGraph(graph, reference)
    assertSameGaussianCanonicalInferenceShape(current, reference)

    return maxMarginalChange(graph, current, reference.marginal)
end

function maxMarginalChange(
    tree::TreeFactorGraph{GaussianFactorGraph},
    current,
    reference
)
    return maxMarginalChange(tree.graph, current, reference)
end

function maxMarginalChange(
    graph::DiscreteFactorGraph,
    current::DiscreteSumProductInference,
    reference::DiscreteSumProductInference
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, current)
    assertDiscreteSumProductInferenceMatchesGraph(graph, reference)
    assertSameDiscreteSumProductInferenceShape(current, reference)

    return maxMarginalChange(graph, current, reference.marginal)
end

function maxMarginalChange(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    current::DiscreteSumProductInference,
    reference::DiscreteSumProductInference
)
    return maxMarginalChange(tree.graph, current, reference)
end

function validateTolerance(tolerance)
    if tolerance === nothing
        return nothing
    end

    if !(tolerance isa Real)
        error("Tolerance must be a real number or nothing.")
    end

    if tolerance < 0
        error("Tolerance must be nonnegative.")
    end

    return nothing
end
