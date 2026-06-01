"""
    AbstractFactorGraph

Abstract supertype for factor graph containers.

# Notes

[`GaussianFactorGraph`](@ref), [`DiscreteFactorGraph`](@ref), and
[`TreeFactorGraph`](@ref) are subtypes of `AbstractFactorGraph`.
"""
abstract type AbstractFactorGraph end

"""
    VariableId

Accepted variable identifiers: `Int` or `Symbol`.

# Notes

Use an id when you want a compact programmatic reference to a variable, for example
`:x1` or `1`.
"""
const VariableId = Union{Int, Symbol}

"""
    VariableRef

Accepted variable references: `Int`, `Symbol`, or `String`.

# Notes

Symbols and integer ids reference variable ids; strings reference variable labels. Most
public functions accept either form, so both `:x1` and `"x1"` can refer to the same
variable when its id is `:x1` and its label is `"x1"`.
"""
const VariableRef = Union{Int, Symbol, String}

"""
    FactorRef

Accepted factor references: `Int`, `Symbol`, or `String`.

# Notes

Integer references select factor indices. Symbols and strings reference factor
labels, so both `:f1` and `"f1"` can refer to the same factor when its label is
`"f1"`.
"""
const FactorRef = Union{Int, Symbol, String}

"""
    ComponentRef

Accepted Gaussian component references: `Int`, `Symbol`, or `String`.

# Notes

Component references are local to a [`GaussianVariable`](@ref). They name the
entries of the variable's continuous state vector.
"""
const ComponentRef = Union{Int, Symbol, String}

"""
    StateRef

Accepted discrete state references: `Int`, `Symbol`, or `String`.

# Notes

State references are local to a [`DiscreteVariable`](@ref). They define the
order used by factor table dimensions.
"""
const StateRef = Union{Int, Symbol, String}

function defaultNodeLabel(id::Union{Int, Symbol, String})
    return replace(string(id), "_" => "")
end

function defaultFactorLabel(id::Int)
    return "f$id"
end

"""
    Edge

Internal graph edge connecting one factor node to one variable node.
"""
struct Edge
    id::Int
    factorIndex::Int
    variableIndex::Int
end

function asVector(data)
    if data isa Number
        return [Float64(data)]
    elseif data isa AbstractVector
        return Vector{Float64}(data)
    else
        error("Input must be a scalar or a vector.")
    end
end

function symmetricPart(matrix::AbstractMatrix)
    matrixFloat = Matrix{Float64}(matrix)

    return 0.5 * (matrixFloat + matrixFloat')
end

function copySymmetricPart!(destination::AbstractMatrix, matrix::AbstractMatrix)
    if size(destination) != size(matrix)
        error("Destination and source matrix dimensions must match.")
    end

    for row in axes(destination, 1)
        for column in axes(destination, 2)
            destination[row, column] = 0.5 * (matrix[row, column] + matrix[column, row])
        end
    end

    return destination
end

function variableIndex(
    referenceIndex::Dict{VariableRef, Int},
    variable::VariableRef
)
    if !haskey(referenceIndex, variable)
        error("Variable reference $variable is not defined.")
    end

    return referenceIndex[variable]
end

function variableRef(variable)
    error("Variable references must be Int, Symbol, String, or Variable.")
end

variableRef(variable::VariableRef) = variable

function variableKind(variables::AbstractVector)
    return string(nameof(eltype(variables)))
end

function validateVariables(variables::AbstractVector)
    if isempty(variables)
        error("At least one $(variableKind(variables)) node must be defined.")
    end

    ids = [variable.id for variable in variables]
    labels = [variable.label for variable in variables]

    if length(unique(ids)) != length(ids)
        error("$(variableKind(variables)) ids must be unique.")
    end

    if length(unique(labels)) != length(labels)
        error("$(variableKind(variables)) labels must be unique.")
    end

    return nothing
end

function validateVariablesIfAny(variables::AbstractVector)
    if isempty(variables)
        return nothing
    end

    return validateVariables(variables)
end

function buildRefIndex(variables::AbstractVector)
    referenceIndex = Dict{VariableRef, Int}()

    for (index, variable) in pairs(variables)
        if haskey(referenceIndex, variable.id)
            error("Duplicate $(variableKind(variables)) id $(variable.id).")
        end

        if haskey(referenceIndex, variable.label)
            error("Duplicate $(variableKind(variables)) label $(variable.label).")
        end

        referenceIndex[variable.id] = index
        referenceIndex[variable.label] = index
    end

    return referenceIndex
end

function factorIndex(factors::AbstractVector, factorLabel::Symbol)
    return factorIndex(factors, string(factorLabel))
end

function factorKind(factors::AbstractVector)
    return string(nameof(eltype(factors)))
end

function factorIndex(factors::AbstractVector, factorRef::Int)
    if factorRef < 1 || factorRef > length(factors)
        error("$(factorKind(factors)) reference $factorRef is not defined.")
    end

    return factorRef
end

function factorIndex(factors::AbstractVector, factorLabel::String)
    for (index, factorData) in pairs(factors)
        if factorData.label == factorLabel
            return index
        end
    end

    error("$(factorKind(factors)) label $factorLabel is not defined.")
end

function hasFactorLabel(graph::AbstractFactorGraph, label::String)
    for factorData in graph.factors
        if factorData.label == label
            return true
        end
    end

    return false
end

function hasVariableId(graph::AbstractFactorGraph, id::VariableId)
    return haskey(graph.referenceIndex, id)
end

function hasVariableLabel(graph::AbstractFactorGraph, label::String)
    return haskey(graph.referenceIndex, label)
end

function bumpTopologyVersion!(graph::AbstractFactorGraph)
    graph.topologyVersion += 1

    return graph
end

"""
    edgeIndex(graph::AbstractFactorGraph; variable::VariableRef, factor::FactorRef)

Return the edge id connecting `variable` and `factor`.

# Arguments

- `graph`: Gaussian or discrete factor graph.

# Keywords

- `variable`: Variable id or label.
- `factor`: Factor index or label.

# Returns

The one-based edge id.

# Notes

Throws an error if the selected variable and factor are not connected.

# Example

```julia
x1 = DiscreteVariable(:x1, 2)
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])

edgeIndex(graph; variable = :x1, factor = "f1")
```
"""
function edgeIndex(
    graph::AbstractFactorGraph;
    variable::VariableRef,
    factor::FactorRef
)
    variableIdx = variableIndex(graph.referenceIndex, variable)
    factorIdx = factorIndex(graph, factor)

    for edgeId in graph.variableEdges[variableIdx]
        edge = graph.edges[edgeId]

        if edge.factorIndex == factorIdx
            return edgeId
        end
    end

    variableData = graph.variables[variableIdx]
    factorData = graph.factors[factorIdx]

    error("Variable $(variableData.label) is not connected to factor $(factorData.label).")
end

"""
    edgeIndices(graph::AbstractFactorGraph; variable = nothing, factor = nothing)

Return edge ids connected to a variable, a factor, or their specific pair.

# Arguments

- `graph`: Gaussian or discrete factor graph.

# Keywords

- `variable`: Optional variable id or label.
- `factor`: Optional factor index or label.

# Returns

A vector of edge ids.

# Notes

Exactly one of `variable` or `factor` selects all adjacent edges for that node.
Passing both selects the single edge between them.

# Example

```julia
x1 = DiscreteVariable(:x1, 2)
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])

edgeIndices(graph; variable = :x1)
```
"""
function edgeIndices(
    graph::AbstractFactorGraph;
    variable = nothing,
    factor = nothing
)
    if variable === nothing && factor === nothing
        error("Specify variable, factor, or both.")
    end

    if variable !== nothing && factor !== nothing
        return [edgeIndex(graph; variable = variable, factor = factor)]
    end

    if variable !== nothing
        if !(variable isa VariableRef)
            error("Variable reference must be Int, Symbol, or String.")
        end

        variableIdx = variableIndex(graph.referenceIndex, variable)

        return copy(graph.variableEdges[variableIdx])
    end

    if !(factor isa FactorRef)
        error("Factor reference must be Int, Symbol, or String.")
    end

    factorIdx = factorIndex(graph, factor)

    return copy(graph.factorEdges[factorIdx])
end

function validateFactorsIfAny(
    factors::AbstractVector,
    variables::AbstractVector,
    referenceIndex::Dict{VariableRef, Int}
)
    if isempty(factors)
        return nothing
    end

    return validateFactors(factors, variables, referenceIndex)
end

function buildEdges(
    variableList::AbstractVector,
    factorList::AbstractVector,
    referenceIndex::Dict{VariableRef, Int}
)
    edges = Edge[]
    factorEdges = [Int[] for _ in factorList]
    variableEdges = [Int[] for _ in variableList]

    edgeId = 0

    for (variableIndexValue, _) in pairs(variableList)
        for (factorIndexValue, factorData) in pairs(factorList)
            for variableRef in factorData.variables
                vIndex = variableIndex(referenceIndex, variableRef)

                if vIndex != variableIndexValue
                    continue
                end

                edgeId += 1

                push!(edges, Edge(edgeId, factorIndexValue, vIndex))
                push!(factorEdges[factorIndexValue], edgeId)
                push!(variableEdges[vIndex], edgeId)
            end
        end
    end

    return edges, factorEdges, variableEdges
end

function addVariableNode!(graph::AbstractFactorGraph, variable)
    if hasVariableId(graph, variable.id)
        error("Variable id $(variable.id) is already defined.")
    end

    if hasVariableLabel(graph, variable.label)
        error("Variable label $(variable.label) is already defined.")
    end

    variableIndexValue = length(graph.variables) + 1

    push!(graph.variables, variable)
    graph.referenceIndex[variable.id] = variableIndexValue
    graph.referenceIndex[variable.label] = variableIndexValue
    push!(graph.variableEdges, Int[])
    bumpTopologyVersion!(graph)

    return variable
end

function initializingUnaryFactor(graph::AbstractFactorGraph, variableIndexValue::Int)
    return initializingUnaryFactor(
        graph.factors,
        graph.variables,
        graph.referenceIndex,
        variableIndexValue
    )
end

"""
    TreeFactorGraph{G <: AbstractFactorGraph}

Tree-oriented view of a factor graph.

`TreeFactorGraph` wraps an existing graph and stores only the tree root,
parent-edge orientation, and traversal orders. It shares the underlying graph;
it does not copy variables, factors, edges, or messages.

# Fields

- `graph`: Underlying Gaussian or discrete factor graph.
- `rootVariableIndex`: Internal root variable index.
- `variableParentEdge`: Parent edge for each variable node.
- `factorParentEdge`: Parent edge for each factor node.
- `forwardOrder`: Leaf-to-root edge order.
- `backwardOrder`: Root-to-leaf edge order.
- `forwardIndex`: Cursor used by [`forwardStep!`](@ref).
- `backwardIndex`: Cursor used by [`backwardStep!`](@ref).

# Notes

The type parameter `G` is the wrapped graph type. For example,
`TreeFactorGraph{DiscreteFactorGraph}` is a tree view over a
[`DiscreteFactorGraph`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
x2 = GaussianVariable(:x2, 1)

f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")
f2 = GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.2; label = "f2")

graph = factorGraph([x1, x2], [f1, f2])
tree = treeFactorGraph(graph; root = :x1)
```
"""
mutable struct TreeFactorGraph{G <: AbstractFactorGraph} <: AbstractFactorGraph
    graph::G
    rootVariableIndex::Int
    variableParentEdge::Vector{Union{Nothing, Int}}
    factorParentEdge::Vector{Union{Nothing, Int}}
    forwardOrder::Vector{Int}
    backwardOrder::Vector{Int}
    forwardIndex::Int
    backwardIndex::Int
end

function treeFactorGraphImpl(
    graph::G;
    root::Union{Nothing, VariableRef} = nothing
) where {G <: AbstractFactorGraph}
    rootVariableIndex =
        root === nothing ? 1 : variableIndex(graph.referenceIndex, root)

    variableCount = length(graph.variables)
    factorCount = length(graph.factors)
    nodeCount = variableCount + factorCount

    if length(graph.edges) != nodeCount - 1
        error(
            "A tree factor graph with $nodeCount nodes must have " *
            "$(nodeCount - 1) edges, but got $(length(graph.edges))."
        )
    end

    visitedVariables = falses(variableCount)
    visitedFactors = falses(factorCount)
    variableParentEdge = Union{Nothing, Int}[nothing for _ in 1:variableCount]
    factorParentEdge = Union{Nothing, Int}[nothing for _ in 1:factorCount]
    rootToLeavesOrder = Int[]

    queue = Tuple{Bool, Int}[(true, rootVariableIndex)]
    visitedVariables[rootVariableIndex] = true
    queueIndex = 1

    while queueIndex <= length(queue)
        isVariable, nodeIndex = queue[queueIndex]
        queueIndex += 1

        if isVariable
            for edgeId in graph.variableEdges[nodeIndex]
                edge = graph.edges[edgeId]

                if visitedFactors[edge.factorIndex]
                    continue
                end

                visitedFactors[edge.factorIndex] = true
                factorParentEdge[edge.factorIndex] = edgeId
                push!(rootToLeavesOrder, edgeId)
                push!(queue, (false, edge.factorIndex))
            end
        else
            for edgeId in graph.factorEdges[nodeIndex]
                edge = graph.edges[edgeId]

                if visitedVariables[edge.variableIndex]
                    continue
                end

                visitedVariables[edge.variableIndex] = true
                variableParentEdge[edge.variableIndex] = edgeId
                push!(rootToLeavesOrder, edgeId)
                push!(queue, (true, edge.variableIndex))
            end
        end
    end

    if !all(visitedVariables) || !all(visitedFactors)
        error("A tree factor graph must be connected.")
    end

    return TreeFactorGraph(
        graph,
        rootVariableIndex,
        variableParentEdge,
        factorParentEdge,
        reverse(rootToLeavesOrder),
        rootToLeavesOrder,
        1,
        1
    )
end

function formattedRefs(references::Vector{VariableRef})
    return "[" * join(string.(references), ", ") * "]"
end
