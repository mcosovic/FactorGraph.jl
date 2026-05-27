"""
    AbstractInference

Abstract supertype for belief-propagation inference states.

# Notes

[`GaussianInference`](@ref) and [`DiscreteInference`](@ref) are subtypes of
`AbstractInference`. Concrete inference states store algorithm-specific messages,
results, update-control flags, and the graph topology version they were built for.
"""
abstract type AbstractInference end

"""
    GaussianInference

Abstract supertype for Gaussian belief-propagation inference states.

# Notes

[`GaussianMomentInference`](@ref), [`GaussianCanonicalInference`](@ref), and
[`GaussianMinSumInference`](@ref) are subtypes of `GaussianInference`, so APIs
that work with any Gaussian inference state can dispatch on this type.
"""
abstract type GaussianInference <: AbstractInference end

"""
    DiscreteInference

Abstract supertype for discrete belief-propagation inference states.

# Notes

[`DiscreteSumProductInference`](@ref) and [`DiscreteMinSumInference`](@ref) are
subtypes of `DiscreteInference`, so APIs that work with any discrete inference
state can dispatch on this type.
"""
abstract type DiscreteInference <: AbstractInference end

"""
    GaussianSumProductInference

Abstract supertype for Gaussian sum-product inference states.

# Notes

[`GaussianMomentInference`](@ref) and [`GaussianCanonicalInference`](@ref) are
subtypes of `GaussianSumProductInference`. Gaussian min-sum inference is not a
sum-product inference state because it stores MAP estimates instead of marginals.
"""
abstract type GaussianSumProductInference <: GaussianInference end

"""
    DiscreteSumProductInference

Inference state for sum-product discrete belief propagation.

# Fields

- `initial`: Initial variable messages.
- `variableToFactor`: Variable-to-factor probability messages.
- `factorToVariable`: Factor-to-variable probability messages.
- `marginal`: Stored variable marginals.
- `nextVariableToFactor`: Flooding buffer for variable-to-factor messages.
- `nextFactorToVariable`: Flooding buffer for factor-to-variable messages.
- `frozenEdges`: Edge freeze flags.
- `frozenVariables`: Variable freeze flags.
- `frozenFactors`: Factor freeze flags.
- `dampedEdges`: Edge damping flags.
- `edgeDampingProb`: Per-edge damping probabilities.
- `edgeDampingAlpha`: Per-edge damping mixing weights.
- `graphVersion`: Topology version used to detect stale inference states.

# Notes

Stores DiscreteVariable-to-DiscreteFactor messages, DiscreteFactor-to-DiscreteVariable
messages, marginals, and update-control flags. Create instances with [`sumproduct`](@ref).

# Example

```julia
x1 = DiscreteVariable(:x1, 2; states = [:off, :on])
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])
inference = sumproduct(graph)
```
"""
mutable struct DiscreteSumProductInference <: DiscreteInference
    initial::Vector{Vector{Float64}}
    variableToFactor::Vector{Vector{Float64}}
    factorToVariable::Vector{Vector{Float64}}
    marginal::Vector{Vector{Float64}}
    nextVariableToFactor::Vector{Vector{Float64}}
    nextFactorToVariable::Vector{Vector{Float64}}
    frozenEdges::Vector{Bool}
    frozenVariables::Vector{Bool}
    frozenFactors::Vector{Bool}
    dampedEdges::Vector{Bool}
    edgeDampingProb::Vector{Float64}
    edgeDampingAlpha::Vector{Float64}
    graphVersion::Int
end

"""
    AbstractSumProductInference

Union alias for sum-product belief-propagation inference states.

# Notes

`AbstractSumProductInference` includes [`GaussianSumProductInference`](@ref) and
[`DiscreteSumProductInference`](@ref). It is a union alias rather than an
abstract supertype because Julia types have a single immediate supertype, and
[`DiscreteSumProductInference`](@ref) already belongs to the
[`DiscreteInference`](@ref) hierarchy.
"""
const AbstractSumProductInference = Union{GaussianSumProductInference, DiscreteSumProductInference}

"""
    QuadraticMessage

Quadratic message used by Gaussian min-sum belief propagation.

# Fields

- `J`: Quadratic precision matrix.
- `h`: Linear information vector.
- `c`: Constant offset.

# Notes

The message stores a quadratic cost over one variable. `J` controls curvature,
`h` controls the linear term, and `c` stores a constant offset that shifts the
cost without changing the minimizing value.

# Example

```julia
message = QuadraticMessage([1.0;;], [0.0], 0.0)
```
"""
mutable struct QuadraticMessage
    J::Matrix{Float64}
    h::Vector{Float64}
    c::Float64
end

"""
    GaussianMinSumInference

Inference state for Gaussian min-sum belief propagation.

# Fields

- `variableToFactor`: Variable-to-factor quadratic messages.
- `factorToVariable`: Factor-to-variable quadratic messages.
- `estimate`: Stored MAP estimates.
- `nextVariableToFactor`: Flooding buffer for variable-to-factor messages.
- `nextFactorToVariable`: Flooding buffer for factor-to-variable messages.
- `frozenEdges`: Edge freeze flags.
- `frozenVariables`: Variable freeze flags.
- `frozenFactors`: Factor freeze flags.
- `dampedEdges`: Edge damping flags.
- `edgeDampingProb`: Per-edge damping probabilities.
- `edgeDampingAlpha`: Per-edge damping mixing weights.
- `defaultMean`: Default mean used for warm-start variable additions.
- `defaultCovariance`: Default covariance used for warm-start variable additions.
- `graphVersion`: Topology version used to detect stale inference states.

# Notes

Gaussian min-sum is the negative-log MAP form of Gaussian max-product belief
propagation. It stores quadratic cost messages and MAP estimates, not
marginal variances. Create instances with [`minsum`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = minsum(graph)
```
"""
mutable struct GaussianMinSumInference <: GaussianInference
    variableToFactor::Vector{QuadraticMessage}
    factorToVariable::Vector{QuadraticMessage}
    estimate::Vector{Vector{Float64}}
    nextVariableToFactor::Vector{QuadraticMessage}
    nextFactorToVariable::Vector{QuadraticMessage}
    frozenEdges::Vector{Bool}
    frozenVariables::Vector{Bool}
    frozenFactors::Vector{Bool}
    dampedEdges::Vector{Bool}
    edgeDampingProb::Vector{Float64}
    edgeDampingAlpha::Vector{Float64}
    defaultMean
    defaultCovariance
    graphVersion::Int
end

"""
    DiscreteMinSumInference

Inference state for discrete min-sum MAP belief propagation.

# Fields

- `initial`: Initial variable cost messages.
- `variableToFactor`: Variable-to-factor cost messages.
- `factorToVariable`: Factor-to-variable cost messages.
- `estimate`: Stored MAP state estimates.
- `nextVariableToFactor`: Flooding buffer for variable-to-factor messages.
- `nextFactorToVariable`: Flooding buffer for factor-to-variable messages.
- `frozenEdges`: Edge freeze flags.
- `frozenVariables`: Variable freeze flags.
- `frozenFactors`: Factor freeze flags.
- `dampedEdges`: Edge damping flags.
- `edgeDampingProb`: Per-edge damping probabilities.
- `edgeDampingAlpha`: Per-edge damping mixing weights.
- `graphVersion`: Topology version used to detect stale inference states.

# Notes

Stores DiscreteVariable-to-DiscreteFactor cost messages, DiscreteFactor-to-DiscreteVariable
cost messages, MAP estimates, and update-control flags. Create instances with
[`minsum`](@ref).

# Example

```julia
x1 = DiscreteVariable(:x1, 2; states = [:off, :on])
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])
inference = minsum(graph)
```
"""
mutable struct DiscreteMinSumInference <: DiscreteInference
    initial::Vector{Vector{Float64}}
    variableToFactor::Vector{Vector{Float64}}
    factorToVariable::Vector{Vector{Float64}}
    estimate::Vector{StateRef}
    nextVariableToFactor::Vector{Vector{Float64}}
    nextFactorToVariable::Vector{Vector{Float64}}
    frozenEdges::Vector{Bool}
    frozenVariables::Vector{Bool}
    frozenFactors::Vector{Bool}
    dampedEdges::Vector{Bool}
    edgeDampingProb::Vector{Float64}
    edgeDampingAlpha::Vector{Float64}
    graphVersion::Int
end

function defaultGaussianMean(variable::GaussianVariable, mean)
    if mean isa Number
        return fill(Float64(mean), variable.dimension)
    elseif mean isa AbstractVector
        meanVector = Vector{Float64}(mean)

        if length(meanVector) != variable.dimension
            error(
                "Default mean dimension is not valid for GaussianVariable $(variable.label). " *
                "Expected $(variable.dimension), but got $(length(meanVector))."
            )
        end

        return meanVector
    else
        error("Default mean must be a scalar or a vector.")
    end
end

function defaultGaussianCovariance(variable::GaussianVariable, covariance)
    if covariance isa Number
        return Float64(covariance) *
            Matrix{Float64}(I, variable.dimension, variable.dimension)
    end

    covarianceMatrix = symmetricPart(asCovarianceMatrix(covariance))

    if size(covarianceMatrix, 1) != variable.dimension
        error(
            "Default covariance row dimension is not valid for GaussianVariable $(variable.label). " *
            "Expected $(variable.dimension), but got $(size(covarianceMatrix, 1))."
        )
    end

    if size(covarianceMatrix, 2) != variable.dimension
        error(
            "Default covariance column dimension is not valid for GaussianVariable $(variable.label). " *
            "Expected $(variable.dimension), but got $(size(covarianceMatrix, 2))."
        )
    end

    if !isposdef(Symmetric(covarianceMatrix))
        error("Default covariance must be positive definite for GaussianVariable $(variable.label).")
    end

    return covarianceMatrix
end

"""
    addVariable!(
        graph::GaussianFactorGraph, inference::GaussianInference, variable::GaussianVariable;
        mean = inference default, covariance = inference default
    )
    addVariable!(
        graph::GaussianFactorGraph, inference::GaussianInference,
        id::VariableId, dimension::Int;
        label = string(id), components = 1:dimension, mean = nothing, covariance = nothing
    )

Add a GaussianVariable node and extend an existing Gaussian inference state.

Existing messages and marginals are kept as a warm start. The new GaussianVariable is
initialized using its own `mean`/`covariance` if provided, otherwise the
defaults stored in `inference` are used. The returned value is the added
[`GaussianVariable`](@ref).

In the `id`/`dimension` convenience form, `mean = nothing` and
`covariance = nothing` create a GaussianVariable without an explicit prior; the
added variable is then initialized from the inference defaults.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

addVariable!(graph, inference, :x2, 1; label = "x2")
```
"""
function addVariable!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    variable::GaussianVariable;
    mean = inference.defaultMean,
    covariance = inference.defaultCovariance
)
    assertInferenceMatchesGraph(graph, inference)

    added = addVariableNode!(graph, variable)
    extendInferenceForAddedVariable!(
        inference,
        added;
        mean = mean,
        covariance = covariance
    )
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addVariable!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    id::VariableId,
    dimension::Int;
    label::String = string(id),
    components = defaultComponentRefs(dimension),
    mean = nothing,
    covariance = nothing
)
    variable = GaussianVariable(
        id,
        dimension;
        label = label,
        components = components,
        mean = mean,
        covariance = covariance
    )

    return addVariable!(
        graph,
        inference,
        variable
    )
end

"""
    AbstractMinSumInference

Union alias for min-sum belief-propagation inference states.

# Notes

`AbstractMinSumInference` includes [`GaussianMinSumInference`](@ref) and
[`DiscreteMinSumInference`](@ref). Gaussian min-sum stores numeric MAP estimate
vectors, while discrete min-sum stores MAP state references.
"""
const AbstractMinSumInference = Union{GaussianMinSumInference, DiscreteMinSumInference}

"""
    marginals(inference::AbstractSumProductInference)

Return the current marginals stored in a sum-product inference object.
"""
function marginals(inference::AbstractSumProductInference)
    return inference.marginal
end

"""
    marginal(
        graph::AbstractFactorGraph, inference::AbstractSumProductInference,
        variable::VariableRef
    )

Validate and return the current marginal for one variable.
"""
function marginal(
    graph::AbstractFactorGraph,
    inference::AbstractSumProductInference,
    variableRef::VariableRef
)
    if graph isa TreeFactorGraph
        return marginal(graph.graph, inference, variableRef)
    end

    assertInferenceMatchesGraph(graph, inference)

    return inference.marginal[variableIndex(graph, variableRef)]
end

"""
    marginals!(graph::GaussianFactorGraph, inference::GaussianSumProductInference)

Recompute GaussianVariable marginals from the current sum-product messages.
"""
function marginals!(graph::AbstractFactorGraph, inference::AbstractSumProductInference)
    if graph isa TreeFactorGraph
        return marginals!(graph.graph, inference)
    end

    error("marginals! is implemented by concrete sum-product inference methods")
end

"""
    estimates(graph::AbstractFactorGraph, inference::AbstractMinSumInference)

Validate and return all stored MAP estimates from a min-sum inference object.
"""
function estimates(graph::AbstractFactorGraph, inference::AbstractMinSumInference)
    if graph isa TreeFactorGraph
        return estimates(graph.graph, inference)
    end

    assertInferenceMatchesGraph(graph, inference)

    return inference.estimate
end

function estimates!(tree::TreeFactorGraph, inference::AbstractMinSumInference)
    return estimates!(tree.graph, inference)
end

"""
    estimate(
        graph::AbstractFactorGraph, inference::AbstractMinSumInference,
        variable::VariableRef
    )

Validate and return the stored MAP estimate for one variable.
"""
function estimate(
    graph::AbstractFactorGraph,
    inference::AbstractMinSumInference,
    variableRef::VariableRef
)
    if graph isa TreeFactorGraph
        return estimate(graph.graph, inference, variableRef)
    end

    assertInferenceMatchesGraph(graph, inference)

    return inference.estimate[variableIndex(graph.referenceIndex, variableRef)]
end

function maxEstimateChange(tree::TreeFactorGraph, inference::AbstractMinSumInference, previous)
    return maxEstimateChange(tree.graph, inference, previous)
end

function setInferenceGraphVersion!(
    graph::AbstractFactorGraph,
    inference::AbstractInference
)
    inference.graphVersion = graph.topologyVersion

    return inference
end

function assertInferenceCurrent(graph::AbstractFactorGraph, inference::AbstractInference)
    if inference.graphVersion != graph.topologyVersion
        error(
            "$(typeof(inference)) is stale because graph topology changed. " *
            "Create a fresh inference object or use graph/inference update APIs."
        )
    end

    return nothing
end

function assertCommonInferenceStorageMatchesGraph(
    graph::AbstractFactorGraph,
    inference::AbstractInference,
    message::String
)
    if length(inference.variableToFactor) != length(graph.edges) ||
            length(inference.factorToVariable) != length(graph.edges) ||
            length(inference.frozenEdges) != length(graph.edges) ||
            length(inference.frozenVariables) != length(graph.variables) ||
            length(inference.frozenFactors) != length(graph.factors) ||
            length(inference.dampedEdges) != length(graph.edges) ||
            length(inference.edgeDampingProb) != length(graph.edges) ||
            length(inference.edgeDampingAlpha) != length(graph.edges)
        error(message)
    end

    return nothing
end

function assertVariableStorageMatchesGraph(
    graph::AbstractFactorGraph,
    message::String,
    storages...
)
    variableCount = length(graph.variables)

    for storage in storages
        if length(storage) != variableCount
            error(message)
        end
    end

    return nothing
end

function copyMessages(messages::AbstractVector, copyMessageFunction)
    return [copyMessageFunction(message) for message in messages]
end

# Discrete inference helpers

function normalizeDiscreteVector!(data::Vector{Float64}, label::String)
    total = sum(data)

    if !isfinite(total) || total <= 0.0
        error("$label must contain at least one positive finite entry.")
    end

    data ./= total

    return data
end

function normalizedDiscreteVector(data, expectedLength::Int, label::String)
    vector = Vector{Float64}(data)

    if length(vector) != expectedLength
        error("$label length must be $expectedLength, but got $(length(vector)).")
    end

    if any(value -> !isfinite(value) || value < 0.0, vector)
        error("$label entries must be finite and nonnegative.")
    end

    return normalizeDiscreteVector!(vector, label)
end

function uniformDiscreteMessage(variable::DiscreteVariable)
    return fill(1.0 / variable.cardinality, variable.cardinality)
end

function initialMessage(variable::DiscreteVariable)
    if variable.probability === nothing
        return uniformDiscreteMessage(variable)
    end

    return copy(variable.probability)
end

function initialMessageFromUnaryFactor(factorData::DiscreteFactor)
    return normalizedDiscreteVector(
        factorData.table,
        length(factorData.table),
        "Initializing DiscreteFactor table"
    )
end

function normalizeDiscreteCost!(data::Vector{Float64}, label::String)
    offset = minimum(data)

    if !isfinite(offset)
        error("$label must contain at least one finite entry.")
    end

    data .-= offset

    return data
end

function probabilityToCost(probability::Vector{Float64})
    cost = Vector{Float64}(undef, length(probability))

    for index in eachindex(probability)
        cost[index] = probability[index] > 0.0 ? -log(probability[index]) : Inf
    end

    return normalizeDiscreteCost!(cost, "Discrete min-sum initial cost")
end

function tableToCost(table)
    cost = Array{Float64}(undef, size(table))

    for index in eachindex(table)
        value = table[index]
        cost[index] = value > 0.0 ? -log(value) : Inf
    end

    return cost
end

function initialMessage(graph::DiscreteFactorGraph, variableIndexValue::Int)
    factorData = initializingUnaryFactor(graph, variableIndexValue)

    if factorData !== nothing
        return initialMessageFromUnaryFactor(factorData)
    end

    return initialMessage(graph.variables[variableIndexValue])
end

function initialCost(graph::DiscreteFactorGraph, variableIndexValue::Int)
    return probabilityToCost(initialMessage(graph, variableIndexValue))
end

function copyDiscreteMessages(messages::Vector{Vector{Float64}})
    return copyMessages(messages, copy)
end

function currentEstimate(
    graph::DiscreteFactorGraph,
    costs::Vector{Float64},
    variableIndexValue::Int
)
    stateIdx = argmin(costs)

    return graph.variables[variableIndexValue].states[stateIdx]
end

function validateDamping(prob::Float64, alpha::Float64)
    if prob < 0.0 || prob > 1.0
        error("Damping probability prob must be between 0 and 1.")
    end

    if alpha < 0.0 || alpha > 1.0
        error("Damping weight alpha must be between 0 and 1.")
    end

    return nothing
end

function validateDiscreteDamping(prob::Float64, alpha::Float64)
    return validateDamping(prob, alpha)
end

function edgeDimension(graph::DiscreteFactorGraph, edge::Edge)
    return graph.variables[edge.variableIndex].cardinality
end

function ensureMessageDimensions(
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference}
)
    for edge in graph.edges
        dimension = edgeDimension(graph, edge)

        if length(inference.variableToFactor[edge.id]) != dimension ||
                length(inference.factorToVariable[edge.id]) != dimension
            error("$(typeof(inference)) message dimensions do not match the graph.")
        end
    end

    return nothing
end

function factorVariableIndices(graph::DiscreteFactorGraph, factorData::DiscreteFactor)
    indices = Vector{Int}(undef, length(factorData.variables))

    for index in eachindex(factorData.variables)
        indices[index] = variableIndex(graph.referenceIndex, factorData.variables[index])
    end

    return indices
end

function factorEdgeIds(
    graph::DiscreteFactorGraph,
    factorIndexValue::Int,
    variableIndices::Vector{Int}
)
    edgeIds = Vector{Int}(undef, length(variableIndices))

    for position in eachindex(variableIndices)
        variableIdx = variableIndices[position]

        for edgeId in graph.factorEdges[factorIndexValue]
            edge = graph.edges[edgeId]

            if edge.variableIndex == variableIdx
                edgeIds[position] = edgeId
                break
            end
        end
    end

    return edgeIds
end

function variableToFactorMessage!(
    output::Vector{Float64},
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    edge::Edge
)
    hasIncomingMessage = false
    fill!(output, 1.0)

    for otherEdgeId in graph.variableEdges[edge.variableIndex]
        if otherEdgeId == edge.id
            continue
        end

        hasIncomingMessage = true
        output .*= inference.factorToVariable[otherEdgeId]
    end

    if !hasIncomingMessage
        copyto!(output, inference.initial[edge.variableIndex])
    end

    return normalizeDiscreteVector!(output, "DiscreteVariable-to-DiscreteFactor message")
end

function variableToFactorMessage!(
    output::Vector{Float64},
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    edge::Edge
)
    hasIncomingMessage = false
    fill!(output, 0.0)

    for otherEdgeId in graph.variableEdges[edge.variableIndex]
        if otherEdgeId == edge.id
            continue
        end

        hasIncomingMessage = true
        output .+= inference.factorToVariable[otherEdgeId]
    end

    if !hasIncomingMessage
        copyto!(output, inference.initial[edge.variableIndex])
    end

    return normalizeDiscreteCost!(output, "DiscreteVariable-to-DiscreteFactor cost message")
end

function variableToFactorMessage!(
    output::Vector{Float64},
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    edgeId::Int
)
    return variableToFactorMessage!(output, graph, inference, graph.edges[edgeId])
end

function factorToVariableMessage!(
    output::Vector{Float64},
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    edge::Edge
)
    factorData = graph.factors[edge.factorIndex]
    variableIndices = factorVariableIndices(graph, factorData)
    factorEdgeIdsValue = factorEdgeIds(graph, edge.factorIndex, variableIndices)
    targetPosition = findfirst(==(edge.variableIndex), variableIndices)

    fill!(output, 0.0)

    for index in CartesianIndices(factorData.table)
        weight = factorData.table[index]

        for position in eachindex(variableIndices)
            if position == targetPosition
                continue
            end

            incomingEdgeId = factorEdgeIdsValue[position]
            weight *= inference.variableToFactor[incomingEdgeId][index[position]]
        end

        output[index[targetPosition]] += weight
    end

    return normalizeDiscreteVector!(output, "DiscreteFactor-to-DiscreteVariable message")
end

function factorToVariableMessage!(
    output::Vector{Float64},
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    edge::Edge
)
    factorData = graph.factors[edge.factorIndex]
    variableIndices = factorVariableIndices(graph, factorData)
    factorEdgeIdsValue = factorEdgeIds(graph, edge.factorIndex, variableIndices)
    targetPosition = findfirst(==(edge.variableIndex), variableIndices)
    factorCost = tableToCost(factorData.table)

    fill!(output, Inf)

    for index in CartesianIndices(factorData.table)
        cost = factorCost[index]

        for position in eachindex(variableIndices)
            if position == targetPosition
                continue
            end

            incomingEdgeId = factorEdgeIdsValue[position]
            cost += inference.variableToFactor[incomingEdgeId][index[position]]
        end

        targetStateIndex = index[targetPosition]
        output[targetStateIndex] = min(output[targetStateIndex], cost)
    end

    return normalizeDiscreteCost!(output, "DiscreteFactor-to-DiscreteVariable cost message")
end

function factorToVariableMessage!(
    output::Vector{Float64},
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    edgeId::Int
)
    return factorToVariableMessage!(output, graph, inference, graph.edges[edgeId])
end

function updateVariableToFactorMessages!(
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    output::Vector{Vector{Float64}}
)
    for edge in graph.edges
        if inference.frozenVariables[edge.variableIndex] || inference.frozenEdges[edge.id]
            copyto!(output[edge.id], inference.variableToFactor[edge.id])
        else
            variableToFactorMessage!(output[edge.id], graph, inference, edge)
        end
    end

    return output
end

function updateFactorToVariableMessages!(
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    output::Vector{Vector{Float64}}
)
    for edge in graph.edges
        if inference.frozenFactors[edge.factorIndex] || inference.frozenEdges[edge.id]
            copyto!(output[edge.id], inference.factorToVariable[edge.id])
        else
            factorToVariableMessage!(output[edge.id], graph, inference, edge)
        end
    end

    return output
end

function commitFactorToVariableMessages!(
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference}
)
    inference.factorToVariable, inference.nextFactorToVariable =
        inference.nextFactorToVariable, inference.factorToVariable

    return nothing
end

function commitVariableToFactorMessages!(
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference}
)
    inference.variableToFactor, inference.nextVariableToFactor =
        inference.nextVariableToFactor, inference.variableToFactor

    return nothing
end

function dampMessage!(
    message::Vector{Float64},
    previous::Vector{Float64},
    alpha::Float64
)
    message .= alpha .* previous .+ (1.0 - alpha) .* message

    return normalizeDiscreteVector!(message, "Damped discrete message")
end

function dampCostMessage!(
    message::Vector{Float64},
    previous::Vector{Float64},
    alpha::Float64
)
    message .= alpha .* previous .+ (1.0 - alpha) .* message

    return normalizeDiscreteCost!(message, "Damped discrete cost message")
end

function applyFactorToVariableDamping!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    output::Vector{Vector{Float64}},
    previous::Vector{Vector{Float64}};
    damping::Bool,
    prob::Float64,
    alpha::Float64,
    rng
)
    for edge in graph.edges
        shouldDamp = damping
        edgeProb = prob
        edgeAlpha = alpha

        if inference.dampedEdges[edge.id]
            shouldDamp = true
            edgeProb = inference.edgeDampingProb[edge.id]
            edgeAlpha = inference.edgeDampingAlpha[edge.id]
        end

        if shouldDamp && rand(rng) <= edgeProb
            dampMessage!(output[edge.id], previous[edge.id], edgeAlpha)
        end
    end

    return nothing
end

function applyFactorToVariableDamping!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    output::Vector{Vector{Float64}},
    previous::Vector{Vector{Float64}};
    damping::Bool,
    prob::Float64,
    alpha::Float64,
    rng
)
    for edge in graph.edges
        shouldDamp = damping
        edgeProb = prob
        edgeAlpha = alpha

        if inference.dampedEdges[edge.id]
            shouldDamp = true
            edgeProb = inference.edgeDampingProb[edge.id]
            edgeAlpha = inference.edgeDampingAlpha[edge.id]
        end

        if shouldDamp && rand(rng) <= edgeProb
            dampCostMessage!(output[edge.id], previous[edge.id], edgeAlpha)
        end
    end

    return nothing
end

function unsupportedDiscreteControlError(inference)
    error("Freezing and damping controls do not support $(typeof(inference)).")
end

function addVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    variable::DiscreteVariable
)
    error("Unsupported discrete inference type $(typeof(inference)).")
end

function addVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    id::VariableId,
    cardinality::Int;
    label::String = string(id),
    states = defaultStateRefs(cardinality),
    probability = nothing
)
    return addVariable!(
        graph,
        inference,
        DiscreteVariable(
            id,
            cardinality;
            label = label,
            states = states,
            probability = probability
        )
    )
end

function addFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    factorData::DiscreteFactor
)
    error("Unsupported discrete inference type $(typeof(inference)).")
end

function addFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    args...;
    label::String = "",
    initialize::Bool = false
)
    return addFactor!(
        graph,
        inference,
        DiscreteFactor(args...; label = label, initialize = initialize)
    )
end

function updateFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    factorRef::FactorRef;
    table = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    error("Unsupported discrete inference type $(typeof(inference)).")
end

function updateFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference;
    factor::FactorRef,
    table = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    return updateFactor!(
        graph,
        inference,
        factor;
        table = table,
        initialize = initialize
    )
end

function isFrozenFactor(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    factorRef::FactorRef
)
    unsupportedDiscreteControlError(inference)
end

function isFrozenVariable(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    variableRef::VariableRef
)
    unsupportedDiscreteControlError(inference)
end

function isFrozenEdge(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference;
    variable::VariableRef,
    factor::FactorRef
)
    unsupportedDiscreteControlError(inference)
end

function freezeFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    factorRef::FactorRef
)
    unsupportedDiscreteControlError(inference)
end

function freezeVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    variableRef::VariableRef
)
    unsupportedDiscreteControlError(inference)
end

function freezeEdge!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference;
    variable::VariableRef,
    factor::FactorRef
)
    unsupportedDiscreteControlError(inference)
end

function unfreezeFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    factorRef::FactorRef
)
    unsupportedDiscreteControlError(inference)
end

function unfreezeVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference,
    variableRef::VariableRef
)
    unsupportedDiscreteControlError(inference)
end

function unfreezeEdge!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference;
    variable::VariableRef,
    factor::FactorRef
)
    unsupportedDiscreteControlError(inference)
end

function areDampedEdges(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    unsupportedDiscreteControlError(inference)
end

function dampEdges!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4
)
    unsupportedDiscreteControlError(inference)
end

function undampEdges!(
    graph::DiscreteFactorGraph,
    inference::DiscreteInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    unsupportedDiscreteControlError(inference)
end

const GaussianControlledInference = Union{GaussianSumProductInference, GaussianMinSumInference}
const DiscreteControlledInference = Union{DiscreteSumProductInference, DiscreteMinSumInference}
const ControlledInference = Union{GaussianControlledInference, DiscreteControlledInference}

function copyStoredMessage!(output, source)
    return copyMessage!(output, source)
end

function copyStoredMessage!(output::Vector{Float64}, source::Vector{Float64})
    copyto!(output, source)

    return output
end

function copyFactorToVariableMessages!(
    graph::AbstractFactorGraph,
    output::AbstractVector,
    source::AbstractVector,
    factorIdx::Int
)
    for edgeId in graph.factorEdges[factorIdx]
        copyStoredMessage!(output[edgeId], source[edgeId])
    end

    return nothing
end

function copyVariableToFactorMessages!(
    graph::AbstractFactorGraph,
    output::AbstractVector,
    source::AbstractVector,
    variableIdx::Int
)
    for edgeId in graph.variableEdges[variableIdx]
        copyStoredMessage!(output[edgeId], source[edgeId])
    end

    return nothing
end

function isFrozenFactorImpl(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    factorRef::FactorRef
)
    assertInferenceMatchesGraph(graph, inference)

    return inference.frozenFactors[factorIndex(graph, factorRef)]
end

function freezeFactorImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    factorRef::FactorRef
)
    assertInferenceMatchesGraph(graph, inference)

    factorIdx = factorIndex(graph, factorRef)
    inference.frozenFactors[factorIdx] = true

    if length(inference.nextFactorToVariable) == length(inference.factorToVariable)
        copyFactorToVariableMessages!(
            graph,
            inference.nextFactorToVariable,
            inference.factorToVariable,
            factorIdx
        )
    end

    return nothing
end

function freezeVariableImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variableRef::VariableRef
)
    assertInferenceMatchesGraph(graph, inference)

    variableIdx = variableIndex(graph.referenceIndex, variableRef)
    inference.frozenVariables[variableIdx] = true

    if length(inference.nextVariableToFactor) == length(inference.variableToFactor)
        copyVariableToFactorMessages!(
            graph,
            inference.nextVariableToFactor,
            inference.variableToFactor,
            variableIdx
        )
    end

    return nothing
end

function freezeEdgeImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variable::VariableRef,
    factor::FactorRef
)
    assertInferenceMatchesGraph(graph, inference)

    edgeId = edgeIndex(graph; variable = variable, factor = factor)
    inference.frozenEdges[edgeId] = true

    if length(inference.nextFactorToVariable) == length(inference.factorToVariable)
        copyStoredMessage!(
            inference.nextFactorToVariable[edgeId],
            inference.factorToVariable[edgeId]
        )
    end

    if length(inference.nextVariableToFactor) == length(inference.variableToFactor)
        copyStoredMessage!(
            inference.nextVariableToFactor[edgeId],
            inference.variableToFactor[edgeId]
        )
    end

    return nothing
end

function isFrozenVariableImpl(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variableRef::VariableRef
)
    assertInferenceMatchesGraph(graph, inference)

    return inference.frozenVariables[variableIndex(graph.referenceIndex, variableRef)]
end

function isFrozenEdgeImpl(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variable::VariableRef,
    factor::FactorRef
)
    assertInferenceMatchesGraph(graph, inference)

    return inference.frozenEdges[edgeIndex(graph; variable = variable, factor = factor)]
end

function unfreezeFactorImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    factorRef::FactorRef
)
    assertInferenceMatchesGraph(graph, inference)
    inference.frozenFactors[factorIndex(graph, factorRef)] = false

    return nothing
end

function unfreezeVariableImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variableRef::VariableRef
)
    assertInferenceMatchesGraph(graph, inference)
    inference.frozenVariables[variableIndex(graph.referenceIndex, variableRef)] = false

    return nothing
end

function unfreezeEdgeImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variable::VariableRef,
    factor::FactorRef
)
    assertInferenceMatchesGraph(graph, inference)
    inference.frozenEdges[edgeIndex(graph; variable = variable, factor = factor)] = false

    return nothing
end

function areDampedEdgesImpl(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variable::Union{Nothing, VariableRef},
    factor::Union{Nothing, FactorRef}
)
    assertInferenceMatchesGraph(graph, inference)

    edgeIds = edgeIndices(graph; variable = variable, factor = factor)

    return all(edgeId -> inference.dampedEdges[edgeId], edgeIds)
end

function selectedOrAllEdgeIds(
    graph::AbstractFactorGraph,
    variable::Union{Nothing, VariableRef},
    factor::Union{Nothing, FactorRef}
)
    if variable === nothing && factor === nothing
        return eachindex(graph.edges)
    end

    return edgeIndices(graph; variable = variable, factor = factor)
end

function dampEdgesImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variable::Union{Nothing, VariableRef},
    factor::Union{Nothing, FactorRef},
    prob::Float64,
    alpha::Float64
)
    assertInferenceMatchesGraph(graph, inference)
    validateDamping(prob, alpha)

    for edgeId in selectedOrAllEdgeIds(graph, variable, factor)
        inference.dampedEdges[edgeId] = true
        inference.edgeDampingProb[edgeId] = prob
        inference.edgeDampingAlpha[edgeId] = alpha
    end

    return nothing
end

function undampEdgesImpl!(
    graph::AbstractFactorGraph,
    inference::ControlledInference,
    variable::Union{Nothing, VariableRef},
    factor::Union{Nothing, FactorRef}
)
    assertInferenceMatchesGraph(graph, inference)

    for edgeId in selectedOrAllEdgeIds(graph, variable, factor)
        inference.dampedEdges[edgeId] = false
    end

    return nothing
end

function isFrozenFactor(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference,
    factorRef::FactorRef
)
    return isFrozenFactorImpl(graph, inference, factorRef)
end

function isFrozenFactor(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference,
    factorRef::FactorRef
)
    return isFrozenFactorImpl(graph, inference, factorRef)
end

function isFrozenVariable(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference,
    variableRef::VariableRef
)
    return isFrozenVariableImpl(graph, inference, variableRef)
end

function isFrozenVariable(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference,
    variableRef::VariableRef
)
    return isFrozenVariableImpl(graph, inference, variableRef)
end

function isFrozenEdge(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return isFrozenEdgeImpl(graph, inference, variable, factor)
end

function isFrozenEdge(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return isFrozenEdgeImpl(graph, inference, variable, factor)
end

function freezeFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference,
    factorRef::FactorRef
)
    return freezeFactorImpl!(graph, inference, factorRef)
end

function freezeFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference,
    factorRef::FactorRef
)
    return freezeFactorImpl!(graph, inference, factorRef)
end

function freezeVariable!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference,
    variableRef::VariableRef
)
    return freezeVariableImpl!(graph, inference, variableRef)
end

function freezeVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference,
    variableRef::VariableRef
)
    return freezeVariableImpl!(graph, inference, variableRef)
end

function freezeEdge!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return freezeEdgeImpl!(graph, inference, variable, factor)
end

function freezeEdge!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return freezeEdgeImpl!(graph, inference, variable, factor)
end

function unfreezeFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference,
    factorRef::FactorRef
)
    return unfreezeFactorImpl!(graph, inference, factorRef)
end

function unfreezeFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference,
    factorRef::FactorRef
)
    return unfreezeFactorImpl!(graph, inference, factorRef)
end

function unfreezeVariable!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference,
    variableRef::VariableRef
)
    return unfreezeVariableImpl!(graph, inference, variableRef)
end

function unfreezeVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference,
    variableRef::VariableRef
)
    return unfreezeVariableImpl!(graph, inference, variableRef)
end

function unfreezeEdge!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return unfreezeEdgeImpl!(graph, inference, variable, factor)
end

function unfreezeEdge!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return unfreezeEdgeImpl!(graph, inference, variable, factor)
end

function areDampedEdges(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    return areDampedEdgesImpl(graph, inference, variable, factor)
end

function areDampedEdges(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    return areDampedEdgesImpl(graph, inference, variable, factor)
end

function dampEdges!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4
)
    return dampEdgesImpl!(graph, inference, variable, factor, prob, alpha)
end

function dampEdges!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4
)
    return dampEdgesImpl!(graph, inference, variable, factor, prob, alpha)
end

function undampEdges!(
    graph::GaussianFactorGraph,
    inference::GaussianControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    return undampEdgesImpl!(graph, inference, variable, factor)
end

function undampEdges!(
    graph::DiscreteFactorGraph,
    inference::DiscreteControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    return undampEdgesImpl!(graph, inference, variable, factor)
end

