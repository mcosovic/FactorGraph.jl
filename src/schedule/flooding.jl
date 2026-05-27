"""
    FloodingSchedule

Stateless flooding message schedule.

# Notes

Use [`floodingSchedule`](@ref) to create one for low-level [`messages!`](@ref) calls.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

schedule = floodingSchedule(graph, inference)
messages!(graph, inference, schedule)
```
"""
struct FloodingSchedule end

"""
    floodingSchedule(graph::AbstractFactorGraph, inference::AbstractInference)

Create a flooding schedule for low-level message updates.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Returns

A [`FloodingSchedule`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

schedule = floodingSchedule(graph, inference)
messages!(graph, inference, schedule)
```
"""
function floodingSchedule(
    graph::GaussianFactorGraph,
    inference::GaussianInference
)
    assertInferenceMatchesGraph(graph, inference)
    return FloodingSchedule()
end

function floodingSchedule(
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference}
)
    assertInferenceMatchesGraph(graph, inference)
    return FloodingSchedule()
end

function floodingSchedule(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference
)
    return floodingSchedule(tree.graph, inference)
end

function floodingSchedule(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference}
)
    return floodingSchedule(tree.graph, inference)
end

function messages!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    schedule::FloodingSchedule;
    kwargs...
)
    return messages!(graph, inference; schedule = :flooding, kwargs...)
end

function messages!(
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    schedule::FloodingSchedule;
    kwargs...
)
    return messages!(graph, inference; schedule = :flooding, kwargs...)
end

function messages!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::FloodingSchedule;
    kwargs...
)
    return messages!(tree.graph, inference, schedule; kwargs...)
end

function messages!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    schedule::FloodingSchedule;
    kwargs...
)
    return messages!(tree.graph, inference, schedule; kwargs...)
end

function ensureFloodingBuffers!(inference::GaussianMomentInference)
    if length(inference.nextVariableToFactor) == length(inference.variableToFactor) &&
            length(inference.nextFactorToVariable) == length(inference.factorToVariable)
        return inference
    end

    inference.nextVariableToFactor = copyGaussianMessages(inference.variableToFactor)
    inference.nextFactorToVariable = copyGaussianMessages(inference.factorToVariable)

    return inference
end

function ensureFloodingBuffers!(inference::GaussianCanonicalInference)
    if length(inference.nextVariableToFactor) == length(inference.variableToFactor) &&
            length(inference.nextFactorToVariable) == length(inference.factorToVariable)
        return inference
    end

    inference.nextVariableToFactor = copyCanonicalMessages(inference.variableToFactor)
    inference.nextFactorToVariable = copyCanonicalMessages(inference.factorToVariable)

    return inference
end

function ensureFloodingBuffers!(inference::GaussianMinSumInference)
    if length(inference.nextVariableToFactor) == length(inference.variableToFactor) &&
            length(inference.nextFactorToVariable) == length(inference.factorToVariable)
        return inference
    end

    inference.nextVariableToFactor = copyQuadraticMessages(inference.variableToFactor)
    inference.nextFactorToVariable = copyQuadraticMessages(inference.factorToVariable)

    return inference
end

function ensureFloodingBuffers!(inference::DiscreteSumProductInference)
    if length(inference.nextVariableToFactor) == length(inference.variableToFactor) &&
            length(inference.nextFactorToVariable) == length(inference.factorToVariable)
        return inference
    end

    inference.nextVariableToFactor = copyDiscreteMessages(inference.variableToFactor)
    inference.nextFactorToVariable = copyDiscreteMessages(inference.factorToVariable)

    return inference
end

function ensureFloodingBuffers!(inference::DiscreteMinSumInference)
    if length(inference.nextVariableToFactor) == length(inference.variableToFactor) &&
            length(inference.nextFactorToVariable) == length(inference.factorToVariable)
        return inference
    end

    inference.nextVariableToFactor = copyDiscreteMessages(inference.variableToFactor)
    inference.nextFactorToVariable = copyDiscreteMessages(inference.factorToVariable)

    return inference
end