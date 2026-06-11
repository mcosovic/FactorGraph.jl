function graphFigureInferenceHighlights(
    graph::AbstractFactorGraph,
    inference::AbstractInference,
    residual,
    variance,
    highlight::AbstractVector
)
    figureHighlight = graphFigureVarianceHighlights(graph, inference, variance)
    append!(figureHighlight, graphFigureResidualHighlights(graph, inference, residual))
    append!(figureHighlight, highlight)

    return figureHighlight
end

function graphFigureMarginalTooltip(
    graph::AbstractFactorGraph,
    inference::Union{Nothing, AbstractInference},
    variableIndex::Int,
    detail::Symbol
)
    return Pair{String, Vector{String}}[]
end

function graphFigureMarginalTooltip(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    variableIndex::Int,
    detail::Symbol
)
    variable = graph.variables[variableIndex]
    mean = marginalMean(graph, inference, variable.id)
    covariance = marginalCovariance(graph, inference, variable.id)
    lines = graphFigureGaussianMarginalTooltipLines(variable, mean, covariance, detail)

    return ["Marginal" => lines]
end

function graphFigureMarginalTooltip(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    variableIndex::Int,
    detail::Symbol
)
    variable = graph.variables[variableIndex]
    probabilities = [
        marginalProbability(graph, inference, variable.id, state)
        for state in variable.states
    ]
    lines = graphFigureDiscreteMarginalTooltipLines(variable, probabilities, detail)

    return ["Marginal" => lines]
end

function graphFigureGaussianMarginalTooltipLines(
    variable::GaussianVariable,
    mean,
    covariance,
    detail::Symbol
)
    if variable.dimension == 1
        return [
            "mean: $(graphFigureValueText(graphFigureScalarMarginalValue(mean)))",
            "variance: $(graphFigureValueText(covariance[1, 1]))"
        ]
    end

    if detail == :full
        return [
            graphFigureFieldText("mean", mean),
            graphFigureFieldText("covariance", covariance)
        ]
    end

    return [
        "mean size: $(graphFigureSizeText(mean))",
        "covariance size: $(graphFigureSizeText(covariance))"
    ]
end

function graphFigureDiscreteMarginalTooltipLines(
    variable::DiscreteVariable,
    probabilities::AbstractVector{Float64},
    detail::Symbol
)
    confidence, stateIndex = findmax(probabilities)
    lines = String[]

    if detail == :full
        push!(lines, graphFigureFieldText("probability", probabilities))
    end

    push!(lines, "most likely: $(graphFigureValueText(variable.states[stateIndex]))")
    push!(lines, "confidence: $(graphFigureValueText(confidence))")

    return lines
end

function graphFigureScalarMarginalValue(value::AbstractVector)
    return only(value)
end

function graphFigureScalarMarginalValue(value)
    return value
end

function graphFigureResidualHighlights(
    graph::AbstractFactorGraph,
    inference::AbstractInference,
    residual
)
    if residual === nothing
        return NamedTuple[]
    end

    assertInferenceMatchesGraph(graph, inference)

    candidateValues = graphFigureResidualCandidates(graph, inference)
    selected = graphFigureSelectedResidualCandidates(candidateValues, residual)
    residuals = graphFigureEdgeResidualValues(graph, selected)

    return graphFigureResidualEdgeHighlights(graph, residuals)
end

function graphFigureEdgeResidualValues(
    graph::AbstractFactorGraph,
    candidates::AbstractVector
)
    residuals = zeros(Float64, length(graph.edges))

    for candidate in candidates
        residuals[candidate.edge] = max(residuals[candidate.edge], candidate.value)
    end

    return residuals
end

function graphFigureResidualCandidates(
    graph::AbstractFactorGraph,
    inference::AbstractInference
)
    candidates = NamedTuple[]

    for edge in graph.edges
        push!(
            candidates,
            (
                edge = edge.id,
                value = candidateResidual(graph, inference, edge.id, true)
            )
        )
    end

    for edge in graph.edges
        push!(
            candidates,
            (
                edge = edge.id,
                value = candidateResidual(graph, inference, edge.id, false)
            )
        )
    end

    return candidates
end

function graphFigureEdgeResidual(
    graph::AbstractFactorGraph,
    inference::AbstractInference,
    edge::Edge
)
    return max(
        candidateResidual(graph, inference, edge.id, true),
        candidateResidual(graph, inference, edge.id, false)
    )
end

function graphFigureSelectedResidualCandidates(
    candidates::AbstractVector,
    residual::Symbol
)
    if residual != :all
        error("Graph figure residual must be nothing, :all, a positive Int, or a fraction.")
    end

    return [candidate for candidate in candidates if candidate.value > 0.0]
end

function graphFigureSelectedResidualCandidates(
    candidates::AbstractVector,
    residual::Bool
)
    error("Graph figure residual must be nothing, :all, a positive Int, or a fraction.")
end

function graphFigureSelectedResidualCandidates(
    candidates::AbstractVector,
    residual::Integer
)
    if residual < 1
        error("Graph figure residual count must be positive.")
    end

    ordering = sortperm([candidate.value for candidate in candidates]; rev = true)
    count = min(Int(residual), length(ordering))

    return [
        candidates[index] for index in ordering[1:count] if candidates[index].value > 0.0
    ]
end

function graphFigureSelectedResidualCandidates(
    candidates::AbstractVector,
    residual::Real
)
    fraction = Float64(residual)
    if fraction <= 0.0 || fraction > 1.0
        error("Graph figure residual fraction must be greater than 0 and at most 1.")
    end

    count = max(1, ceil(Int, fraction * length(candidates)))

    return graphFigureSelectedResidualCandidates(candidates, count)
end

function graphFigureResidualEdgeHighlights(
    graph::AbstractFactorGraph,
    residualValues::AbstractVector{Float64}
)
    maxResidual = maximum(residualValues; init = 0.0)
    if maxResidual <= 0.0
        return NamedTuple[]
    end

    highlights = NamedTuple[]
    for edge in graph.edges
        residual = residualValues[edge.id]
        if residual > 0.0
            scale = residual / maxResidual
            push!(
                highlights,
                (
                    edge = edge.id,
                    haloStroke = "#dc2626",
                    haloStrokeWidth = 10.0,
                    haloStrokeOpacity = 0.08 + 0.35 * scale
                )
            )
        end
    end

    return highlights
end

function graphFigureVarianceHighlights(
    graph::AbstractFactorGraph,
    inference::AbstractInference,
    variance
)
    if variance === nothing
        return NamedTuple[]
    end

    error("Graph figure variance view is only available for Gaussian sum-product inference.")
end

function graphFigureVarianceHighlights(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    variance
)
    if variance === nothing
        return NamedTuple[]
    end

    assertInferenceMatchesGraph(graph, inference)

    candidates = graphFigureVarianceCandidates(graph, inference)
    selected = graphFigureSelectedVarianceCandidates(candidates, variance)
    values = graphFigureVariableVarianceValues(graph, selected)

    return graphFigureVarianceNodeHighlights(graph, values)
end

function graphFigureVarianceCandidates(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference
)
    candidates = NamedTuple[]

    for variableIndex in eachindex(graph.variables)
        push!(
            candidates,
            (
                variable = variableIndex,
                value = graphFigureVariableVariance(graph, inference, variableIndex)
            )
        )
    end

    return candidates
end

function graphFigureVariableVariance(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    variableIndex::Int
)
    covariance = marginalCovariance(graph, inference, graph.variables[variableIndex].id)

    return maximum(diag(covariance); init = 0.0)
end

function graphFigureVariableVarianceValues(
    graph::GaussianFactorGraph,
    candidates::AbstractVector
)
    values = zeros(Float64, length(graph.variables))

    for candidate in candidates
        values[candidate.variable] = candidate.value
    end

    return values
end

function graphFigureSelectedVarianceCandidates(
    candidates::AbstractVector,
    variance::Symbol
)
    if variance != :all
        error("Graph figure variance must be nothing, :all, a positive Int, or a fraction.")
    end

    return [candidate for candidate in candidates if candidate.value > 0.0]
end

function graphFigureSelectedVarianceCandidates(
    candidates::AbstractVector,
    variance::Bool
)
    error("Graph figure variance must be nothing, :all, a positive Int, or a fraction.")
end

function graphFigureSelectedVarianceCandidates(
    candidates::AbstractVector,
    variance::Integer
)
    if variance < 1
        error("Graph figure variance count must be positive.")
    end

    ordering = sortperm([candidate.value for candidate in candidates]; rev = true)
    count = min(Int(variance), length(ordering))

    return [
        candidates[index] for index in ordering[1:count] if candidates[index].value > 0.0
    ]
end

function graphFigureSelectedVarianceCandidates(
    candidates::AbstractVector,
    variance::Real
)
    fraction = Float64(variance)
    if fraction <= 0.0 || fraction > 1.0
        error("Graph figure variance fraction must be greater than 0 and at most 1.")
    end

    count = max(1, ceil(Int, fraction * length(candidates)))

    return graphFigureSelectedVarianceCandidates(candidates, count)
end

function graphFigureVarianceNodeHighlights(
    graph::GaussianFactorGraph,
    varianceValues::AbstractVector{Float64}
)
    maxVariance = maximum(varianceValues; init = 0.0)
    if maxVariance <= 0.0
        return NamedTuple[]
    end

    highlights = NamedTuple[]
    for variableIndex in eachindex(graph.variables)
        value = varianceValues[variableIndex]
        if value > 0.0
            scale = value / maxVariance
            push!(
                highlights,
                (
                    variable = graph.variables[variableIndex].id,
                    incidentEdges = false,
                    haloStroke = "#f59e0b",
                    haloStrokeWidth = 10.0,
                    haloStrokeOpacity = 0.08 + 0.35 * scale
                )
            )
        end
    end

    return highlights
end
