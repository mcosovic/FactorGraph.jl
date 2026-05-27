"""
    ForwardBackwardSchedule

Mutable cursor used by forward/backward tree inference.

# Fields

- `forwardIndex`: Next position in the forward edge order.
- `backwardIndex`: Next position in the backward edge order.

# Notes

Use [`reset!`](@ref) to restart a schedule or a tree's default step cursor.

The schedule stores the next edge position for the forward pass and the next
edge position for the backward pass. It lets callers step through tree inference
manually with [`forwardStep!`](@ref) and [`backwardStep!`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)
inference = moment(graph)

schedule = forwardBackwardSchedule(tree)
forwardBackward!(tree, inference, schedule)
```
"""
mutable struct ForwardBackwardSchedule
    forwardIndex::Int
    backwardIndex::Int

    function ForwardBackwardSchedule()
        return new(1, 1)
    end
end

"""
    forwardBackwardSchedule(tree::TreeFactorGraph)

Create a new forward/backward schedule for a tree factor graph.

# Arguments

- `tree`: Tree factor graph.

# Returns

A [`ForwardBackwardSchedule`](@ref).

# Notes

The returned schedule starts before the first forward and backward messages.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)
inference = moment(graph)

schedule = forwardBackwardSchedule(tree)
forwardBackward!(tree, inference, schedule)
```
"""
function forwardBackwardSchedule(tree::TreeFactorGraph)
    return ForwardBackwardSchedule()
end

"""
    reset!(schedule::ForwardBackwardSchedule)

Reset a [`ForwardBackwardSchedule`](@ref) to the first forward and backward
message.

# Arguments

- `schedule`: Schedule to reset.

# Returns

The same `schedule`.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)

schedule = forwardBackwardSchedule(tree)
reset!(schedule)
```
"""
function reset!(schedule::ForwardBackwardSchedule)
    schedule.forwardIndex = 1
    schedule.backwardIndex = 1

    return schedule
end

function nextForwardEdge!(
    schedule::ForwardBackwardSchedule,
    tree::TreeFactorGraph
)
    if schedule.forwardIndex > length(tree.forwardOrder)
        return nothing
    end

    edgeId = tree.forwardOrder[schedule.forwardIndex]
    schedule.forwardIndex += 1

    return edgeId
end

function nextBackwardEdge!(
    schedule::ForwardBackwardSchedule,
    tree::TreeFactorGraph
)
    if schedule.backwardIndex > length(tree.backwardOrder)
        return nothing
    end

    edgeId = tree.backwardOrder[schedule.backwardIndex]
    schedule.backwardIndex += 1

    return edgeId
end