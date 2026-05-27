"""
    reset!(tree::TreeFactorGraph)

Reset a [`TreeFactorGraph`](@ref)'s default forward/backward step cursor.

# Arguments

- `tree`: Tree factor graph.

# Returns

The same `tree`.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)

reset!(tree)
```
"""
function reset!(tree::TreeFactorGraph)
    tree.forwardIndex = 1
    tree.backwardIndex = 1

    return tree
end

"""
    refresh!(tree::TreeFactorGraph; root = nothing)

Recompute the orientation and message order of a [`TreeFactorGraph`](@ref).

# Arguments

- `tree`: Tree factor graph to refresh.

# Keywords

- `root`: Optional replacement root variable. If omitted, the previous root is reused.

# Returns

The same `tree`.

# Notes

Use this after graph topology changes when the underlying graph is still a
tree. If `root` is omitted, the previous root variable is reused. Returns the
same `tree` object.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)

refresh!(tree; root = :x1)
```
"""
function refresh!(
    tree::TreeFactorGraph;
    root::Union{Nothing, VariableRef} = nothing
)
    currentRoot = graphVariables(tree.graph)[tree.rootVariableIndex].id
    refreshed = treeFactorGraph(
        tree.graph;
        root = root === nothing ? currentRoot : root
    )

    tree.rootVariableIndex = refreshed.rootVariableIndex
    tree.variableParentEdge = refreshed.variableParentEdge
    tree.factorParentEdge = refreshed.factorParentEdge
    tree.forwardOrder = refreshed.forwardOrder
    tree.backwardOrder = refreshed.backwardOrder
    tree.forwardIndex = refreshed.forwardIndex
    tree.backwardIndex = refreshed.backwardIndex

    return tree
end

function variableIndex(
    tree::TreeFactorGraph,
    variableRef::VariableRef
)
    return variableIndex(tree.graph, variableRef)
end

function variableDimension(
    tree::TreeFactorGraph,
    variableRef::VariableRef
)
    return variableDimension(tree.graph, variableRef)
end

function componentIndex(
    tree::TreeFactorGraph{GaussianFactorGraph},
    variableRef::VariableRef,
    component::ComponentRef
)
    return componentIndex(tree.graph, variableRef, component)
end

function componentValue(
    tree::TreeFactorGraph{GaussianFactorGraph},
    variableRef::VariableRef,
    index::Int
)
    return componentValue(tree.graph, variableRef, index)
end

function stateIndex(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    variableRef::VariableRef,
    state::StateRef
)
    return stateIndex(tree.graph, variableRef, state)
end

function stateValue(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    variableRef::VariableRef,
    index::Int
)
    return stateValue(tree.graph, variableRef, index)
end

function factorIndex(tree::TreeFactorGraph, factorRef::FactorRef)
    return factorIndex(tree.graph, factorRef)
end

function edgeIndex(tree::TreeFactorGraph; kwargs...)
    return edgeIndex(tree.graph; kwargs...)
end

function edgeIndices(tree::TreeFactorGraph; kwargs...)
    return edgeIndices(tree.graph; kwargs...)
end

function addVariable!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    args...;
    kwargs...
)
    return addVariable!(tree.graph, args...; kwargs...)
end

function addFactor!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    args...;
    kwargs...
)
    added = addFactor!(tree.graph, args...; kwargs...)
    refresh!(tree)

    return added
end

function updateFactor!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    args...;
    kwargs...
)
    return updateFactor!(tree.graph, args...; kwargs...)
end

function updateFactor!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    args...;
    kwargs...
)
    return updateFactor!(tree.graph, args...; kwargs...)
end

function moment(tree::TreeFactorGraph{GaussianFactorGraph}; kwargs...)
    return moment(tree.graph; kwargs...)
end

function canonical(tree::TreeFactorGraph{GaussianFactorGraph}; kwargs...)
    return canonical(tree.graph; kwargs...)
end

function minsum(tree::TreeFactorGraph{GaussianFactorGraph}; kwargs...)
    return minsum(tree.graph; kwargs...)
end

function sumproduct(tree::TreeFactorGraph{DiscreteFactorGraph}; kwargs...)
    return sumproduct(tree.graph; kwargs...)
end

function minsum(tree::TreeFactorGraph{DiscreteFactorGraph}; kwargs...)
    return minsum(tree.graph; kwargs...)
end

function coefficientBlock(tree::TreeFactorGraph{GaussianFactorGraph}; kwargs...)
    return coefficientBlock(tree.graph; kwargs...)
end

function coefficientBlocks(tree::TreeFactorGraph{GaussianFactorGraph}, factorRef::FactorRef)
    return coefficientBlocks(tree.graph, factorRef)
end

function marginalMean(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference,
    variableRef::VariableRef
)
    return marginalMean(tree.graph, inference, variableRef)
end

function marginalCovariance(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference,
    variableRef::VariableRef
)
    return marginalCovariance(tree.graph, inference, variableRef)
end

function isFrozenFactor(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference,
    factorRef::FactorRef
)
    return isFrozenFactor(tree.graph, inference, factorRef)
end

function isFrozenVariable(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference,
    variableRef::VariableRef
)
    return isFrozenVariable(tree.graph, inference, variableRef)
end

function isFrozenEdge(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return isFrozenEdge(tree.graph, inference; variable = variable, factor = factor)
end

function freezeFactor!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference,
    factorRef::FactorRef
)
    return freezeFactor!(tree.graph, inference, factorRef)
end

function freezeVariable!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference,
    variableRef::VariableRef
)
    return freezeVariable!(tree.graph, inference, variableRef)
end

function freezeEdge!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return freezeEdge!(tree.graph, inference; variable = variable, factor = factor)
end

function unfreezeFactor!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference,
    factorRef::FactorRef
)
    return unfreezeFactor!(tree.graph, inference, factorRef)
end

function unfreezeVariable!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference,
    variableRef::VariableRef
)
    return unfreezeVariable!(tree.graph, inference, variableRef)
end

function unfreezeEdge!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return unfreezeEdge!(tree.graph, inference; variable = variable, factor = factor)
end

function isDampedEdge(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return isDampedEdge(tree.graph, inference; variable = variable, factor = factor)
end

function dampEdges!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing,
    prob::Float64 = 1.0,
    alpha::Float64 = 0.5
)
    return dampEdges!(
        tree.graph,
        inference;
        variable = variable,
        factor = factor,
        prob = prob,
        alpha = alpha
    )
end

function undampEdges!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    return undampEdges!(tree.graph, inference; variable = variable, factor = factor)
end

function isFrozenFactor(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    factorRef::FactorRef
)
    return isFrozenFactor(tree.graph, inference, factorRef)
end

function isFrozenVariable(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    variableRef::VariableRef
)
    return isFrozenVariable(tree.graph, inference, variableRef)
end

function isFrozenEdge(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return isFrozenEdge(tree.graph, inference; variable = variable, factor = factor)
end

function freezeFactor!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    factorRef::FactorRef
)
    return freezeFactor!(tree.graph, inference, factorRef)
end

function freezeVariable!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    variableRef::VariableRef
)
    return freezeVariable!(tree.graph, inference, variableRef)
end

function freezeEdge!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return freezeEdge!(tree.graph, inference; variable = variable, factor = factor)
end

function unfreezeFactor!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    factorRef::FactorRef
)
    return unfreezeFactor!(tree.graph, inference, factorRef)
end

function unfreezeVariable!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    variableRef::VariableRef
)
    return unfreezeVariable!(tree.graph, inference, variableRef)
end

function unfreezeEdge!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return unfreezeEdge!(tree.graph, inference; variable = variable, factor = factor)
end

function isDampedEdge(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::VariableRef,
    factor::FactorRef
)
    return isDampedEdge(tree.graph, inference; variable = variable, factor = factor)
end

function dampEdges!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing,
    prob::Float64 = 1.0,
    alpha::Float64 = 0.5
)
    return dampEdges!(
        tree.graph,
        inference;
        variable = variable,
        factor = factor,
        prob = prob,
        alpha = alpha
    )
end

function undampEdges!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    return undampEdges!(tree.graph, inference; variable = variable, factor = factor)
end

function marginalProbability(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference,
    variableRef::VariableRef,
    state::StateRef
)
    return marginalProbability(tree.graph, inference, variableRef, state)
end

function addVariable!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    args...;
    kwargs...
)
    return addVariable!(tree.graph, inference, args...; kwargs...)
end

function addFactor!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    args...;
    kwargs...
)
    added = addFactor!(tree.graph, inference, args...; kwargs...)
    refresh!(tree)

    return added
end

function addVariable!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    args...;
    kwargs...
)
    return addVariable!(tree.graph, inference, args...; kwargs...)
end

function addFactor!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    args...;
    kwargs...
)
    added = addFactor!(tree.graph, inference, args...; kwargs...)
    refresh!(tree)

    return added
end

function updateFactor!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    args...;
    kwargs...
)
    return updateFactor!(tree.graph, inference, args...; kwargs...)
end

function updateFactor!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    args...;
    kwargs...
)
    return updateFactor!(tree.graph, inference, args...; kwargs...)
end

function updateForwardMessage!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    edgeId::Int
)
    assertInferenceMatchesGraph(tree.graph, inference)

    graph = tree.graph
    edge = graph.edges[edgeId]

    if tree.variableParentEdge[edge.variableIndex] == edgeId
        variableToFactorMessage!(
            inference.variableToFactor[edgeId],
            graph,
            inference,
            edgeId
        )
    elseif tree.factorParentEdge[edge.factorIndex] == edgeId
        factorToVariableMessage!(
            inference.factorToVariable[edgeId],
            graph,
            inference,
            edgeId
        )
    else
        error("Edge $edgeId is not oriented in the tree factor graph.")
    end

    return edgeId
end

function updateBackwardMessage!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    edgeId::Int
)
    assertInferenceMatchesGraph(tree.graph, inference)

    graph = tree.graph
    edge = graph.edges[edgeId]

    if tree.variableParentEdge[edge.variableIndex] == edgeId
        factorToVariableMessage!(
            inference.factorToVariable[edgeId],
            graph,
            inference,
            edgeId
        )
    elseif tree.factorParentEdge[edge.factorIndex] == edgeId
        variableToFactorMessage!(
            inference.variableToFactor[edgeId],
            graph,
            inference,
            edgeId
        )
    else
        error("Edge $edgeId is not oriented in the tree factor graph.")
    end

    return edgeId
end

function updateForwardMessage!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference,
    edgeId::Int
)
    assertInferenceMatchesGraph(tree.graph, inference)

    graph = tree.graph
    edge = graph.edges[edgeId]

    if tree.variableParentEdge[edge.variableIndex] == edgeId
        variableToFactorMessage!(
            inference.variableToFactor[edgeId],
            graph,
            inference,
            edge
        )
    elseif tree.factorParentEdge[edge.factorIndex] == edgeId
        factorToVariableMessage!(
            inference.factorToVariable[edgeId],
            graph,
            inference,
            edge
        )
    else
        error("Edge $edgeId is not oriented in the tree factor graph.")
    end

    return edgeId
end

function updateForwardMessage!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference,
    edgeId::Int
)
    assertInferenceMatchesGraph(tree.graph, inference)

    graph = tree.graph
    edge = graph.edges[edgeId]

    if tree.variableParentEdge[edge.variableIndex] == edgeId
        variableToFactorMessage!(inference.variableToFactor[edgeId], graph, inference, edge)
    elseif tree.factorParentEdge[edge.factorIndex] == edgeId
        factorToVariableMessage!(inference.factorToVariable[edgeId], graph, inference, edge)
    else
        error("Edge $edgeId is not oriented in the tree factor graph.")
    end

    return edgeId
end

function updateBackwardMessage!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference,
    edgeId::Int
)
    assertInferenceMatchesGraph(tree.graph, inference)

    graph = tree.graph
    edge = graph.edges[edgeId]

    if tree.variableParentEdge[edge.variableIndex] == edgeId
        factorToVariableMessage!(
            inference.factorToVariable[edgeId],
            graph,
            inference,
            edge
        )
    elseif tree.factorParentEdge[edge.factorIndex] == edgeId
        variableToFactorMessage!(
            inference.variableToFactor[edgeId],
            graph,
            inference,
            edge
        )
    else
        error("Edge $edgeId is not oriented in the tree factor graph.")
    end

    return edgeId
end

function updateBackwardMessage!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference,
    edgeId::Int
)
    assertInferenceMatchesGraph(tree.graph, inference)

    graph = tree.graph
    edge = graph.edges[edgeId]

    if tree.variableParentEdge[edge.variableIndex] == edgeId
        factorToVariableMessage!(inference.factorToVariable[edgeId], graph, inference, edge)
    elseif tree.factorParentEdge[edge.factorIndex] == edgeId
        variableToFactorMessage!(inference.variableToFactor[edgeId], graph, inference, edge)
    else
        error("Edge $edgeId is not oriented in the tree factor graph.")
    end

    return edgeId
end

"""
    forwardStep!(tree::TreeFactorGraph, inference::AbstractInference)
    forwardStep!(
        tree::TreeFactorGraph, inference::AbstractInference, schedule::ForwardBackwardSchedule
    )
    forwardStep!(
        tree::TreeFactorGraph, inference::AbstractInference;
        variable::VariableRef, factor::FactorRef
    )
    forwardStep!(
        tree::TreeFactorGraph, inference::AbstractInference, edgeId::Int
    )

Perform one forward tree message update.

# Arguments

- `tree`: Tree factor graph.
- `inference`: Matching inference object.
- `schedule`: Optional forward/backward schedule.
- `edgeId`: Explicit edge id to update.

# Keywords

- `variable`: Variable id or label for selected-edge update.
- `factor`: Factor index or label for selected-edge update.

# Returns

The updated edge id, or `nothing` when the scheduled pass is complete.

# Notes

Without an explicit edge, the next edge in the forward schedule is updated and
the step cursor is advanced. Returns the updated edge id, or `nothing` when the
forward pass is complete.

The edge order is determined by [`treeFactorGraph`](@ref). A forward pass sends
messages from leaves toward the selected root.

When `schedule` is omitted, the tree's default step cursor is advanced. Use
[`reset!`](@ref) on the tree before reusing that cursor.

Passing an edge id, or selecting an edge by `variable` and `factor`, updates
only that selected forward message and does not advance any cursor.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)
inference = moment(graph)

schedule = forwardBackwardSchedule(tree)
forwardStep!(tree, inference, schedule)
```
"""
function forwardStep!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::ForwardBackwardSchedule
)
    edgeId = nextForwardEdge!(schedule, tree)

    if edgeId === nothing
        return nothing
    end

    return updateForwardMessage!(tree, inference, edgeId)
end

function forwardStep!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    if variable !== nothing || factor !== nothing
        if variable === nothing || factor === nothing
            error("Specify both variable and factor for a selected forward step.")
        end

        edgeId = edgeIndex(tree.graph; variable = variable, factor = factor)

        return forwardStep!(tree, inference, edgeId)
    end

    if tree.forwardIndex > length(tree.forwardOrder)
        return nothing
    end

    edgeId = tree.forwardOrder[tree.forwardIndex]
    tree.forwardIndex += 1

    return updateForwardMessage!(tree, inference, edgeId)
end

function forwardStep!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    edgeId::Int
)
    return updateForwardMessage!(tree, inference, edgeId)
end

function forwardStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    schedule::ForwardBackwardSchedule
)
    edgeId = nextForwardEdge!(schedule, tree)

    if edgeId === nothing
        return nothing
    end

    return updateForwardMessage!(tree, inference, edgeId)
end

function forwardStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    if variable !== nothing || factor !== nothing
        if variable === nothing || factor === nothing
            error("Specify both variable and factor for a selected forward step.")
        end

        edgeId = edgeIndex(tree.graph; variable = variable, factor = factor)

        return forwardStep!(tree, inference, edgeId)
    end

    if tree.forwardIndex > length(tree.forwardOrder)
        return nothing
    end

    edgeId = tree.forwardOrder[tree.forwardIndex]
    tree.forwardIndex += 1

    return updateForwardMessage!(tree, inference, edgeId)
end

function forwardStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    edgeId::Int
)
    return updateForwardMessage!(tree, inference, edgeId)
end

"""
    backwardStep!(tree::TreeFactorGraph, inference::AbstractInference)
    backwardStep!(
        tree::TreeFactorGraph, inference::AbstractInference, schedule::ForwardBackwardSchedule
    )
    backwardStep!(
        tree::TreeFactorGraph, inference::AbstractInference;
        variable::VariableRef, factor::FactorRef
    )
    backwardStep!(
        tree::TreeFactorGraph, inference::AbstractInference, edgeId::Int
    )

Perform one backward tree message update.

# Arguments

- `tree`: Tree factor graph.
- `inference`: Matching inference object.
- `schedule`: Optional forward/backward schedule.
- `edgeId`: Explicit edge id to update.

# Keywords

- `variable`: Variable id or label for selected-edge update.
- `factor`: Factor index or label for selected-edge update.

# Returns

The updated edge id, or `nothing` when the scheduled pass is complete.

# Notes

Without an explicit edge, the next edge in the backward schedule is updated and
the step cursor is advanced. Returns the updated edge id, or `nothing` when the
backward pass is complete.

A backward pass sends messages from the selected root toward the leaves.

When `schedule` is omitted, the tree's default step cursor is advanced. Use
[`reset!`](@ref) on the tree before reusing that cursor.

Passing an edge id, or selecting an edge by `variable` and `factor`, updates
only that selected backward message and does not advance any cursor.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)
inference = moment(graph)

schedule = forwardBackwardSchedule(tree)
backwardStep!(tree, inference, schedule)
```
"""
function backwardStep!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::ForwardBackwardSchedule
)
    edgeId = nextBackwardEdge!(schedule, tree)

    if edgeId === nothing
        return nothing
    end

    return updateBackwardMessage!(tree, inference, edgeId)
end

function backwardStep!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    edgeId::Int
)
    return updateBackwardMessage!(tree, inference, edgeId)
end

function backwardStep!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    if variable !== nothing || factor !== nothing
        if variable === nothing || factor === nothing
            error("Specify both variable and factor for a selected backward step.")
        end

        edgeId = edgeIndex(tree.graph; variable = variable, factor = factor)

        return backwardStep!(tree, inference, edgeId)
    end

    if tree.backwardIndex > length(tree.backwardOrder)
        return nothing
    end

    edgeId = tree.backwardOrder[tree.backwardIndex]
    tree.backwardIndex += 1

    return updateBackwardMessage!(tree, inference, edgeId)
end

function backwardStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    schedule::ForwardBackwardSchedule
)
    edgeId = nextBackwardEdge!(schedule, tree)

    if edgeId === nothing
        return nothing
    end

    return updateBackwardMessage!(tree, inference, edgeId)
end

function backwardStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    edgeId::Int
)
    return updateBackwardMessage!(tree, inference, edgeId)
end

function backwardStep!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    if variable !== nothing || factor !== nothing
        if variable === nothing || factor === nothing
            error("Specify both variable and factor for a selected backward step.")
        end

        edgeId = edgeIndex(tree.graph; variable = variable, factor = factor)

        return backwardStep!(tree, inference, edgeId)
    end

    if tree.backwardIndex > length(tree.backwardOrder)
        return nothing
    end

    edgeId = tree.backwardOrder[tree.backwardIndex]
    tree.backwardIndex += 1

    return updateBackwardMessage!(tree, inference, edgeId)
end

function forward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference
)
    return forward!(tree, inference, forwardBackwardSchedule(tree))
end

"""
    forward!(
        tree::TreeFactorGraph, inference::AbstractInference, schedule::ForwardBackwardSchedule
    )

Run forward tree message updates until the forward pass is complete.

# Arguments

- `tree`: Tree factor graph.
- `inference`: Matching inference object.
- `schedule`: Optional forward/backward schedule.

# Returns

The schedule after it has advanced past the final forward edge.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)
inference = moment(graph)

schedule = forwardBackwardSchedule(tree)
forward!(tree, inference, schedule)
```
"""
function forward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::ForwardBackwardSchedule
)
    while forwardStep!(tree, inference, schedule) !== nothing
    end

    return schedule
end

function forward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference
)
    return forward!(tree, inference, forwardBackwardSchedule(tree))
end

function forward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    schedule::ForwardBackwardSchedule
)
    while forwardStep!(tree, inference, schedule) !== nothing
    end

    return schedule
end

function backward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference
)
    return backward!(tree, inference, forwardBackwardSchedule(tree))
end

"""
    backward!(
        tree::TreeFactorGraph, inference::AbstractInference, schedule::ForwardBackwardSchedule
    )

Run backward tree message updates until the backward pass is complete.

# Arguments

- `tree`: Tree factor graph.
- `inference`: Matching inference object.
- `schedule`: Optional forward/backward schedule.

# Returns

The schedule after it has advanced past the final backward edge.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
tree = treeFactorGraph(graph; root = :x1)
inference = moment(graph)

schedule = forwardBackwardSchedule(tree)
backward!(tree, inference, schedule)
```
"""
function backward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::ForwardBackwardSchedule
)
    while backwardStep!(tree, inference, schedule) !== nothing
    end

    return schedule
end

function backward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference
)
    return backward!(tree, inference, forwardBackwardSchedule(tree))
end

function backward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteControlledInference,
    schedule::ForwardBackwardSchedule
)
    while backwardStep!(tree, inference, schedule) !== nothing
    end

    return schedule
end

function forwardBackward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference
)
    return forwardBackward!(tree, inference, forwardBackwardSchedule(tree))
end

"""
    forwardBackward!(
        tree::TreeFactorGraph, inference::AbstractInference, schedule::ForwardBackwardSchedule
    )

Run a complete forward pass followed by a complete backward pass.

# Arguments

- `tree`: Tree factor graph.
- `inference`: Matching inference object.
- `schedule`: Optional forward/backward schedule.

# Returns

The schedule after both passes complete.

# Notes

On a valid Gaussian tree factor graph, this computes exact marginals in one
forward/backward sweep for moment and canonical inference. For Gaussian
min-sum inference, the same API computes MAP estimates instead.

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
function forwardBackward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianInference,
    schedule::ForwardBackwardSchedule
)
    forward!(tree, inference, schedule)
    backward!(tree, inference, schedule)
    marginals!(tree.graph, inference)

    return schedule
end

function forwardBackward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMinSumInference
)
    return forwardBackward!(tree, inference, forwardBackwardSchedule(tree))
end

function forwardBackward!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMinSumInference,
    schedule::ForwardBackwardSchedule
)
    forward!(tree, inference, schedule)
    backward!(tree, inference, schedule)
    estimates!(tree.graph, inference)

    return schedule
end

function forwardBackward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference
)
    return forwardBackward!(tree, inference, forwardBackwardSchedule(tree))
end

function forwardBackward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference,
    schedule::ForwardBackwardSchedule
)
    forward!(tree, inference, schedule)
    backward!(tree, inference, schedule)
    marginals!(tree.graph, inference)

    return schedule
end

function forwardBackward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference
)
    return forwardBackward!(tree, inference, forwardBackwardSchedule(tree))
end

function forwardBackward!(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference,
    schedule::ForwardBackwardSchedule
)
    forward!(tree, inference, schedule)
    backward!(tree, inference, schedule)
    estimates!(tree.graph, inference)

    return schedule
end
