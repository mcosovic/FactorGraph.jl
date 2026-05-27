"""
    SequentialSchedule

Stateless sequential message schedule.

# Notes

Use [`sequentialSchedule`](@ref) to create one for low-level [`messages!`](@ref) calls.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

schedule = sequentialSchedule(graph, inference)
messages!(graph, inference, schedule)
```
"""
struct SequentialSchedule end

"""
    sequentialSchedule(graph::AbstractFactorGraph, inference::AbstractInference)

Create a sequential schedule for low-level message updates.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Returns

A [`SequentialSchedule`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

schedule = sequentialSchedule(graph, inference)
messages!(graph, inference, schedule)
```
"""
function sequentialSchedule(
    graph::GaussianFactorGraph,
    inference::GaussianInference
)
    assertInferenceMatchesGraph(graph, inference)
    return SequentialSchedule()
end

function sequentialSchedule(
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference}
)
    assertInferenceMatchesGraph(graph, inference)
    return SequentialSchedule()
end

function sequentialSchedule(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference
)
    return sequentialSchedule(tree.graph, inference)
end

function sequentialSchedule(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference}
)
    return sequentialSchedule(tree.graph, inference)
end

function messages!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    schedule::SequentialSchedule;
    kwargs...
)
    return messages!(graph, inference; schedule = :sequential, kwargs...)
end

function messages!(
    graph::DiscreteFactorGraph,
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    schedule::SequentialSchedule;
    kwargs...
)
    return messages!(graph, inference; schedule = :sequential, kwargs...)
end

function messages!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::SequentialSchedule;
    kwargs...
)
    return messages!(tree.graph, inference, schedule; kwargs...)
end

function messages!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::Union{DiscreteSumProductInference, DiscreteMinSumInference},
    schedule::SequentialSchedule;
    kwargs...
)
    return messages!(tree.graph, inference, schedule; kwargs...)
end