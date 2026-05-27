"""
    ResidualSchedule

Mutable schedule state for residual belief propagation.

# Fields

- `edgeIds`: Candidate edge ids.
- `directions`: Candidate directions; `true` is factor-to-variable.
- `residuals`: Current candidate residual values.
- `updateFraction`: Fraction used when `updateCount` is omitted.
- `updateCount`: Optional fixed number of candidates to update.
- `lastUpdated`: Candidate indices updated by the most recent step.

# Notes

Use [`residualSchedule`](@ref) to create a schedule, pass it to
[`messages!`](@ref) to update the largest-residual messages, and use
[`gbp!`](@ref) with `schedule = :residual` for a complete iterative loop.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

schedule = residualSchedule(graph, inference; updateFraction = 0.25)
messages!(graph, inference, schedule)
```
"""
mutable struct ResidualSchedule
    edgeIds::Vector{Int}
    directions::Vector{Bool}
    residuals::Vector{Float64}
    updateFraction::Float64
    updateCount::Union{Nothing, Int}
    lastUpdated::Vector{Int}
end

function validateResidualSelection(updateFraction::Float64, updateCount::Union{Nothing, Int})
    if updateFraction <= 0.0 || updateFraction > 1.0
        error("updateFraction must be greater than 0 and less than or equal to 1.")
    end

    if updateCount !== nothing && updateCount < 1
        error("updateCount must be positive when provided.")
    end

    return nothing
end

function residualCandidateCount(graph::GaussianFactorGraph)
    return 2 * length(graph.edges)
end

function residualCandidateCount(graph::DiscreteFactorGraph)
    return 2 * length(graph.edges)
end

function residualUpdateCount(schedule::ResidualSchedule)
    candidateCount = length(schedule.edgeIds)

    if schedule.updateCount !== nothing
        return min(schedule.updateCount, candidateCount)
    end

    return max(1, ceil(Int, schedule.updateFraction * candidateCount))
end

"""
    residualSchedule(
        graph::AbstractFactorGraph, inference::AbstractInference;
        updateFraction = 0.1, updateCount = nothing
    )

Create a residual schedule for belief propagation.

# Arguments

- `graph`: Gaussian factor graph.
- `inference`: Matching Gaussian inference object.

# Keywords

- `updateFraction`: Fraction of directed candidates to update per step.
- `updateCount`: Fixed number of directed candidates to update per step.

# Returns

A mutable [`ResidualSchedule`](@ref).

# Notes

Each candidate is one directed message update. If `updateCount` is omitted,
each residual step updates the largest `updateFraction` of candidate messages.
If `updateCount` is provided, it updates at most that many messages. The two
selection keywords are mutually exclusive.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

schedule = residualSchedule(graph, inference; updateCount = 4)
```
"""
function residualSchedule(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    updateFraction = nothing,
    updateCount = nothing
)
    assertInferenceMatchesGraph(graph, inference)

    if updateFraction !== nothing && updateCount !== nothing
        error("Specify either updateFraction or updateCount, not both.")
    end

    fraction = updateFraction === nothing ? 0.1 : Float64(updateFraction)
    validateResidualSelection(fraction, updateCount)

    candidateCount = residualCandidateCount(graph)
    edgeIds = Vector{Int}(undef, candidateCount)
    directions = Vector{Bool}(undef, candidateCount)

    candidateIndex = 1

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = true
        candidateIndex += 1
    end

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = false
        candidateIndex += 1
    end

    return ResidualSchedule(
        edgeIds,
        directions,
        zeros(candidateCount),
        fraction,
        updateCount,
        Int[]
    )
end

function residualSchedule(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference;
    kwargs...
)
    return residualSchedule(tree.graph, inference; kwargs...)
end

function residualSchedule(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    updateFraction = nothing,
    updateCount = nothing
)
    assertInferenceMatchesGraph(graph, inference)

    if updateFraction !== nothing && updateCount !== nothing
        error("Specify either updateFraction or updateCount, not both.")
    end

    fraction = updateFraction === nothing ? 0.1 : Float64(updateFraction)
    validateResidualSelection(fraction, updateCount)

    candidateCount = residualCandidateCount(graph)
    edgeIds = Vector{Int}(undef, candidateCount)
    directions = Vector{Bool}(undef, candidateCount)
    candidateIndex = 1

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = true
        candidateIndex += 1
    end

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = false
        candidateIndex += 1
    end

    return ResidualSchedule(
        edgeIds,
        directions,
        zeros(candidateCount),
        fraction,
        updateCount,
        Int[]
    )
end

function residualSchedule(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
    updateFraction = nothing,
    updateCount = nothing
)
    assertInferenceMatchesGraph(graph, inference)

    if updateFraction !== nothing && updateCount !== nothing
        error("Specify either updateFraction or updateCount, not both.")
    end

    fraction = updateFraction === nothing ? 0.1 : Float64(updateFraction)
    validateResidualSelection(fraction, updateCount)

    candidateCount = residualCandidateCount(graph)
    edgeIds = Vector{Int}(undef, candidateCount)
    directions = Vector{Bool}(undef, candidateCount)

    candidateIndex = 1

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = true
        candidateIndex += 1
    end

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = false
        candidateIndex += 1
    end

    return ResidualSchedule(
        edgeIds,
        directions,
        zeros(candidateCount),
        fraction,
        updateCount,
        Int[]
    )
end

function residualSchedule(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    updateFraction = nothing,
    updateCount = nothing
)
    assertInferenceMatchesGraph(graph, inference)

    if updateFraction !== nothing && updateCount !== nothing
        error("Specify either updateFraction or updateCount, not both.")
    end

    fraction = updateFraction === nothing ? 0.1 : Float64(updateFraction)
    validateResidualSelection(fraction, updateCount)

    candidateCount = residualCandidateCount(graph)
    edgeIds = Vector{Int}(undef, candidateCount)
    directions = Vector{Bool}(undef, candidateCount)

    candidateIndex = 1

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = true
        candidateIndex += 1
    end

    for edge in graph.edges
        edgeIds[candidateIndex] = edge.id
        directions[candidateIndex] = false
        candidateIndex += 1
    end

    return ResidualSchedule(
        edgeIds,
        directions,
        zeros(candidateCount),
        fraction,
        updateCount,
        Int[]
    )
end

function residualSchedule(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference;
    kwargs...
)
    return residualSchedule(tree.graph, inference; kwargs...)
end

function residualSchedule(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference;
    kwargs...
)
    return residualSchedule(tree.graph, inference; kwargs...)
end

function momentMessageResidual(current::GaussianMessage, proposed::GaussianMessage)
    return max(
        norm(proposed.mean - current.mean),
        norm(proposed.covariance - current.covariance)
    )
end

function canonicalMessageResidual(current::CanonicalMessage, proposed::CanonicalMessage)
    return max(
        norm(proposed.information - current.information),
        norm(proposed.precision - current.precision)
    )
end

function discreteMessageResidual(current::Vector{Float64}, proposed::Vector{Float64})
    return norm(proposed - current)
end

function quadraticMessageResidual(current::QuadraticMessage, proposed::QuadraticMessage)
    return max(
        norm(proposed.h - current.h),
        norm(proposed.J - current.J)
    )
end

function candidateResidual(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    edgeId::Int,
    factorToVariable::Bool
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return 0.0
        end

        proposed = copyGaussianMessage(inference.factorToVariable[edgeId])
        factorToVariableMessage!(proposed, graph, inference, edgeId)

        return momentMessageResidual(inference.factorToVariable[edgeId], proposed)
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return 0.0
    end

    proposed = copyGaussianMessage(inference.variableToFactor[edgeId])
    variableToFactorMessage!(proposed, graph, inference, edgeId)

    return momentMessageResidual(inference.variableToFactor[edgeId], proposed)
end

function candidateResidual(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    edgeId::Int,
    factorToVariable::Bool
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return 0.0
        end

        proposed = copyCanonicalMessage(inference.factorToVariable[edgeId])
        factorToVariableMessage!(proposed, graph, inference, edgeId)

        return canonicalMessageResidual(inference.factorToVariable[edgeId], proposed)
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return 0.0
    end

    proposed = copyCanonicalMessage(inference.variableToFactor[edgeId])
    variableToFactorMessage!(proposed, graph, inference, edgeId)

    return canonicalMessageResidual(inference.variableToFactor[edgeId], proposed)
end

function candidateResidual(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    edgeId::Int,
    factorToVariable::Bool
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return 0.0
        end

        proposed = copy(inference.factorToVariable[edgeId])
        factorToVariableMessage!(proposed, graph, inference, edge)

        return discreteMessageResidual(inference.factorToVariable[edgeId], proposed)
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return 0.0
    end

    proposed = copy(inference.variableToFactor[edgeId])
    variableToFactorMessage!(proposed, graph, inference, edge)

    return discreteMessageResidual(inference.variableToFactor[edgeId], proposed)
end

function candidateResidual(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    edgeId::Int,
    factorToVariable::Bool
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return 0.0
        end

        proposed = copy(inference.factorToVariable[edgeId])
        factorToVariableMessage!(proposed, graph, inference, edge)

        return discreteMessageResidual(inference.factorToVariable[edgeId], proposed)
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return 0.0
    end

    proposed = copy(inference.variableToFactor[edgeId])
    variableToFactorMessage!(proposed, graph, inference, edge)

    return discreteMessageResidual(inference.variableToFactor[edgeId], proposed)
end

function candidateResidual(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    edgeId::Int,
    factorToVariable::Bool
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return 0.0
        end

        proposed = copyQuadraticMessage(inference.factorToVariable[edgeId])
        factorToVariableMessage!(proposed, graph, inference, edgeId)

        return quadraticMessageResidual(inference.factorToVariable[edgeId], proposed)
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return 0.0
    end

    proposed = copyQuadraticMessage(inference.variableToFactor[edgeId])
    variableToFactorMessage!(proposed, graph, inference, edgeId)

    return quadraticMessageResidual(inference.variableToFactor[edgeId], proposed)
end

function refreshResiduals!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    schedule::ResidualSchedule
)
    assertInferenceMatchesGraph(graph, inference)

    if length(schedule.edgeIds) != residualCandidateCount(graph)
        error("ResidualSchedule does not match the GaussianFactorGraph.")
    end

    for candidateIndex in eachindex(schedule.edgeIds)
        schedule.residuals[candidateIndex] = candidateResidual(
            graph,
            inference,
            schedule.edgeIds[candidateIndex],
            schedule.directions[candidateIndex]
        )
    end

    return schedule
end

function refreshResiduals!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    schedule::ResidualSchedule
)
    assertInferenceMatchesGraph(graph, inference)

    if length(schedule.edgeIds) != residualCandidateCount(graph)
        error("ResidualSchedule does not match the DiscreteFactorGraph.")
    end

    for candidateIndex in eachindex(schedule.edgeIds)
        schedule.residuals[candidateIndex] = candidateResidual(
            graph,
            inference,
            schedule.edgeIds[candidateIndex],
            schedule.directions[candidateIndex]
        )
    end

    return schedule
end

function refreshResiduals!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    schedule::ResidualSchedule
)
    assertInferenceMatchesGraph(graph, inference)

    if length(schedule.edgeIds) != residualCandidateCount(graph)
        error("ResidualSchedule does not match the DiscreteFactorGraph.")
    end

    for candidateIndex in eachindex(schedule.edgeIds)
        schedule.residuals[candidateIndex] = candidateResidual(
            graph,
            inference,
            schedule.edgeIds[candidateIndex],
            schedule.directions[candidateIndex]
        )
    end

    return schedule
end

function refreshResiduals!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    schedule::ResidualSchedule
)
    assertInferenceMatchesGraph(graph, inference)

    if length(schedule.edgeIds) != residualCandidateCount(graph)
        error("ResidualSchedule does not match the GaussianFactorGraph.")
    end

    for candidateIndex in eachindex(schedule.edgeIds)
        schedule.residuals[candidateIndex] = candidateResidual(
            graph,
            inference,
            schedule.edgeIds[candidateIndex],
            schedule.directions[candidateIndex]
        )
    end

    return schedule
end

function selectedResidualCandidates(schedule::ResidualSchedule)
    updateCount = residualUpdateCount(schedule)
    ordering = sortperm(schedule.residuals; rev = true)

    return ordering[1:updateCount]
end

function updateResidualCandidate!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    edgeId::Int,
    factorToVariable::Bool;
    damping::Bool,
    prob::Float64,
    alpha::Float64,
    rng
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return nothing
        end

        previous = copyGaussianMessage(inference.factorToVariable[edgeId])
        factorToVariableMessage!(inference.factorToVariable[edgeId], graph, inference, edgeId)

        shouldDamp = damping
        edgeProb = prob
        edgeAlpha = alpha

        if inference.dampedEdges[edgeId]
            shouldDamp = true
            edgeProb = inference.edgeDampingProb[edgeId]
            edgeAlpha = inference.edgeDampingAlpha[edgeId]
        end

        if shouldDamp && rand(rng) <= edgeProb
            dampMessage!(
                inference.factorToVariable[edgeId],
                previous,
                inference.factorToVariable[edgeId],
                edgeAlpha
            )
        end

        return nothing
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return nothing
    end

    variableToFactorMessage!(inference.variableToFactor[edgeId], graph, inference, edgeId)

    return nothing
end

function updateResidualCandidate!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    edgeId::Int,
    factorToVariable::Bool;
    damping::Bool,
    prob::Float64,
    alpha::Float64,
    rng
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return nothing
        end

        previous = copyCanonicalMessage(inference.factorToVariable[edgeId])
        factorToVariableMessage!(inference.factorToVariable[edgeId], graph, inference, edgeId)

        shouldDamp = damping
        edgeProb = prob
        edgeAlpha = alpha

        if inference.dampedEdges[edgeId]
            shouldDamp = true
            edgeProb = inference.edgeDampingProb[edgeId]
            edgeAlpha = inference.edgeDampingAlpha[edgeId]
        end

        if shouldDamp && rand(rng) <= edgeProb
            dampMessage!(
                inference.factorToVariable[edgeId],
                previous,
                inference.factorToVariable[edgeId],
                edgeAlpha
            )
        end

        return nothing
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return nothing
    end

    variableToFactorMessage!(inference.variableToFactor[edgeId], graph, inference, edgeId)

    return nothing
end

function updateResidualCandidate!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    edgeId::Int,
    factorToVariable::Bool;
    damping::Bool,
    prob::Float64,
    alpha::Float64,
    rng
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return nothing
        end

        previous = copy(inference.factorToVariable[edgeId])
        factorToVariableMessage!(inference.factorToVariable[edgeId], graph, inference, edge)

        shouldDamp = damping
        edgeProb = prob
        edgeAlpha = alpha

        if inference.dampedEdges[edgeId]
            shouldDamp = true
            edgeProb = inference.edgeDampingProb[edgeId]
            edgeAlpha = inference.edgeDampingAlpha[edgeId]
        end

        if shouldDamp && rand(rng) <= edgeProb
            dampMessage!(inference.factorToVariable[edgeId], previous, edgeAlpha)
        end

        return nothing
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return nothing
    end

    variableToFactorMessage!(inference.variableToFactor[edgeId], graph, inference, edge)

    return nothing
end

function updateResidualCandidate!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    edgeId::Int,
    factorToVariable::Bool;
    damping::Bool,
    prob::Float64,
    alpha::Float64,
    rng
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return nothing
        end

        previous = copy(inference.factorToVariable[edgeId])
        factorToVariableMessage!(inference.factorToVariable[edgeId], graph, inference, edge)

        shouldDamp = damping
        edgeProb = prob
        edgeAlpha = alpha

        if inference.dampedEdges[edgeId]
            shouldDamp = true
            edgeProb = inference.edgeDampingProb[edgeId]
            edgeAlpha = inference.edgeDampingAlpha[edgeId]
        end

        if shouldDamp && rand(rng) <= edgeProb
            dampCostMessage!(inference.factorToVariable[edgeId], previous, edgeAlpha)
        end

        return nothing
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return nothing
    end

    variableToFactorMessage!(inference.variableToFactor[edgeId], graph, inference, edge)

    return nothing
end

function updateResidualCandidate!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    edgeId::Int,
    factorToVariable::Bool;
    damping::Bool,
    prob::Float64,
    alpha::Float64,
    rng
)
    edge = graph.edges[edgeId]

    if factorToVariable
        if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
            return nothing
        end

        previous = copyQuadraticMessage(inference.factorToVariable[edgeId])
        factorToVariableMessage!(inference.factorToVariable[edgeId], graph, inference, edgeId)

        shouldDamp = damping
        edgeProb = prob
        edgeAlpha = alpha

        if inference.dampedEdges[edgeId]
            shouldDamp = true
            edgeProb = inference.edgeDampingProb[edgeId]
            edgeAlpha = inference.edgeDampingAlpha[edgeId]
        end

        if shouldDamp && rand(rng) <= edgeProb
            dampMessage!(
                inference.factorToVariable[edgeId],
                previous,
                inference.factorToVariable[edgeId],
                edgeAlpha
            )
        end

        return nothing
    end

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return nothing
    end

    variableToFactorMessage!(inference.variableToFactor[edgeId], graph, inference, edgeId)

    return nothing
end

function validateResidualDamping(graph::GaussianFactorGraph, prob::Float64, alpha::Float64)
    return validateDamping(prob, alpha)
end

function validateResidualDamping(graph::DiscreteFactorGraph, prob::Float64, alpha::Float64)
    return validateDiscreteDamping(prob, alpha)
end

"""
    residualStep!(
        graph::AbstractFactorGraph, inference::AbstractInference, schedule::ResidualSchedule;
        damping = false, prob = 0.6, alpha = 0.4, rng = Random.default_rng()
    )

Recompute residuals for all directed messages and update the selected
largest-residual batch in place.

# Arguments

- `graph`: Gaussian factor graph.
- `inference`: Matching Gaussian inference object.
- `schedule`: Residual schedule created for the same graph and inference object.

# Keywords

- `damping`: Enable damping for Gaussian residual updates.
- `prob`: Damping probability.
- `alpha`: Damping mixing weight.
- `rng`: Random number generator used by damping.

# Returns

The updated `schedule`.

# Notes

The schedule's `lastUpdated` field stores the selected candidate indices from
the most recent step. Returns the same `schedule`. This is equivalent to
calling [`messages!`](@ref) with the schedule as the third argument.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

schedule = residualSchedule(graph, inference)
residualStep!(graph, inference, schedule)
```
"""
function residualStep!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    schedule::ResidualSchedule;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    rng = Random.default_rng()
)
    assertInferenceMatchesGraph(graph, inference)
    validateResidualDamping(graph, prob, alpha)
    refreshResiduals!(graph, inference, schedule)

    empty!(schedule.lastUpdated)

    for candidateIndex in selectedResidualCandidates(schedule)
        updateResidualCandidate!(
            graph,
            inference,
            schedule.edgeIds[candidateIndex],
            schedule.directions[candidateIndex];
            damping = damping,
            prob = prob,
            alpha = alpha,
            rng = rng
        )
        push!(schedule.lastUpdated, candidateIndex)
    end

    return schedule
end

function residualStep!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    schedule::ResidualSchedule;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    rng = Random.default_rng()
)
    assertInferenceMatchesGraph(graph, inference)
    validateResidualDamping(graph, prob, alpha)
    refreshResiduals!(graph, inference, schedule)

    empty!(schedule.lastUpdated)

    for candidateIndex in selectedResidualCandidates(schedule)
        updateResidualCandidate!(
            graph,
            inference,
            schedule.edgeIds[candidateIndex],
            schedule.directions[candidateIndex];
            damping = damping,
            prob = prob,
            alpha = alpha,
            rng = rng
        )
        push!(schedule.lastUpdated, candidateIndex)
    end

    return schedule
end

"""
    messages!(
        graph::AbstractFactorGraph,
        inference::AbstractInference,
        schedule::ResidualSchedule;
        damping = false,
        prob = 0.6,
        alpha = 0.4,
        rng = Random.default_rng()
    )

Run one residual-scheduled message update step without recomputing marginals.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.
- `schedule`: Residual schedule created for the same graph and inference object.

# Keywords

- `damping`: Enable damping for residual updates that support it.
- `prob`: Damping probability.
- `alpha`: Damping mixing weight.
- `rng`: Random number generator used by damping.

# Returns

The updated `schedule`.

# Notes

This is the low-level residual scheduling entry point. It recomputes residuals
for all directed messages, updates the selected largest-residual batch, and
stores the selected candidate indices in `schedule.lastUpdated`. Call
[`marginals!`](@ref) afterwards when you want marginals to reflect the latest
messages.
"""
function messages!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    schedule::ResidualSchedule;
    broadcast::Bool = false,
    kwargs...
)
    if broadcast
        error("broadcast = true is not supported with ResidualSchedule.")
    end

    return residualStep!(graph, inference, schedule; kwargs...)
end

function messages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    schedule::ResidualSchedule;
    kwargs...
)
    return residualStep!(graph, inference, schedule; kwargs...)
end

function messages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    schedule::ResidualSchedule;
    kwargs...
)
    return residualStep!(graph, inference, schedule; kwargs...)
end

function residualStep!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::ResidualSchedule;
    kwargs...
)
    return residualStep!(tree.graph, inference, schedule; kwargs...)
end

function residualStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference,
    schedule::ResidualSchedule;
    kwargs...
)
    return residualStep!(tree.graph, inference, schedule; kwargs...)
end

function residualStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference,
    schedule::ResidualSchedule;
    kwargs...
)
    return residualStep!(tree.graph, inference, schedule; kwargs...)
end

function messages!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference,
    schedule::ResidualSchedule;
    broadcast::Bool = false,
    kwargs...
)
    return messages!(tree.graph, inference, schedule; broadcast = broadcast, kwargs...)
end

function messages!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference,
    schedule::ResidualSchedule;
    kwargs...
)
    return messages!(tree.graph, inference, schedule; kwargs...)
end

function messages!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference,
    schedule::ResidualSchedule;
    kwargs...
)
    return messages!(tree.graph, inference, schedule; kwargs...)
end

function runResidualGBP!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    iterations::Int = 10,
    tolerance = nothing,
    updateFraction = nothing,
    updateCount = nothing,
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    rng = Random.default_rng()
)
    assertInferenceMatchesGraph(graph, inference)

    if iterations < 0
        error("Number of iterations must be nonnegative.")
    end

    validateTolerance(tolerance)
    validateResidualDamping(graph, prob, alpha)

    schedule = residualSchedule(
        graph,
        inference;
        updateFraction = updateFraction,
        updateCount = updateCount
    )

    for _ in 1:iterations
        refreshResiduals!(graph, inference, schedule)

        if tolerance !== nothing && maximum(schedule.residuals) <= tolerance
            break
        end

        empty!(schedule.lastUpdated)

        for candidateIndex in selectedResidualCandidates(schedule)
            updateResidualCandidate!(
                graph,
                inference,
                schedule.edgeIds[candidateIndex],
                schedule.directions[candidateIndex];
                damping = damping,
                prob = prob,
                alpha = alpha,
                rng = rng
            )
            push!(schedule.lastUpdated, candidateIndex)
        end

        if inference isa GaussianMinSumInference
            estimates!(graph, inference)
        else
            marginals!(graph, inference)
        end
    end

    return schedule
end

function runResidualGBP!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
    iterations::Int = 10,
    tolerance = nothing,
    updateFraction = nothing,
    updateCount = nothing,
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    rng = Random.default_rng()
)
    assertInferenceMatchesGraph(graph, inference)

    if iterations < 0
        error("Number of iterations must be nonnegative.")
    end

    validateTolerance(tolerance)
    validateResidualDamping(graph, prob, alpha)

    schedule = residualSchedule(
        graph,
        inference;
        updateFraction = updateFraction,
        updateCount = updateCount
    )

    for _ in 1:iterations
        refreshResiduals!(graph, inference, schedule)

        if tolerance !== nothing && maximum(schedule.residuals) <= tolerance
            break
        end

        empty!(schedule.lastUpdated)

        for candidateIndex in selectedResidualCandidates(schedule)
            updateResidualCandidate!(
                graph,
                inference,
                schedule.edgeIds[candidateIndex],
                schedule.directions[candidateIndex];
                damping = damping,
                prob = prob,
                alpha = alpha,
                rng = rng
            )
            push!(schedule.lastUpdated, candidateIndex)
        end

        marginals!(graph, inference)
    end

    return schedule
end

function runResidualGBP!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    iterations::Int = 10,
    tolerance = nothing,
    updateFraction = nothing,
    updateCount = nothing,
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    rng = Random.default_rng()
)
    assertInferenceMatchesGraph(graph, inference)

    if iterations < 0
        error("Number of iterations must be nonnegative.")
    end

    validateTolerance(tolerance)
    validateResidualDamping(graph, prob, alpha)

    schedule = residualSchedule(
        graph,
        inference;
        updateFraction = updateFraction,
        updateCount = updateCount
    )

    for _ in 1:iterations
        refreshResiduals!(graph, inference, schedule)

        if tolerance !== nothing && maximum(schedule.residuals) <= tolerance
            break
        end

        empty!(schedule.lastUpdated)

        for candidateIndex in selectedResidualCandidates(schedule)
            updateResidualCandidate!(
                graph,
                inference,
                schedule.edgeIds[candidateIndex],
                schedule.directions[candidateIndex];
                damping = damping,
                prob = prob,
                alpha = alpha,
                rng = rng
            )
            push!(schedule.lastUpdated, candidateIndex)
        end

        estimates!(graph, inference)
    end

    return schedule
end

function runResidualGBP!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference;
    kwargs...
)
    return runResidualGBP!(tree.graph, inference; kwargs...)
end

function runResidualGBP!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference;
    kwargs...
)
    return runResidualGBP!(tree.graph, inference; kwargs...)
end

function runResidualGBP!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference;
    kwargs...
)
    return runResidualGBP!(tree.graph, inference; kwargs...)
end
