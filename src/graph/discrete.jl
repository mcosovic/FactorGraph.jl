function asDiscreteTable(data)
    if !(data isa AbstractArray)
        error("DiscreteFactor table must be an array.")
    end

    table = Array{Float64}(data)

    if isempty(table)
        error("DiscreteFactor table must be nonempty.")
    end

    if any(table .< 0.0)
        error("DiscreteFactor table entries must be nonnegative.")
    end

    return table
end

function asDiscreteProbability(data, cardinality::Int, label::String)
    probability = asVector(data)

    if length(probability) != cardinality
        error(
            "Probability dimension is not valid for DiscreteVariable $label. " *
            "Expected $cardinality, but got $(length(probability))."
        )
    end

    if any(probability .< 0.0)
        error("Probability entries must be nonnegative for DiscreteVariable $label.")
    end

    total = sum(probability)

    if total <= 0.0
        error("Probability entries must have positive sum for DiscreteVariable $label.")
    end

    return probability ./ total
end

function defaultStateRefs(cardinality::Int)
    return StateRef[1:cardinality...]
end

function asStateRefs(states)
    if !(states isa AbstractVector)
        error("DiscreteVariable states must be a vector.")
    end

    if isempty(states)
        error("DiscreteVariable states must be nonempty.")
    end

    stateRefs = StateRef[]

    for state in states
        if !(state isa StateRef)
            error("DiscreteVariable states must be Int, Symbol, or String.")
        end

        if state isa String && isempty(state)
            error("DiscreteVariable state must be nonempty.")
        end

        if state in stateRefs
            error("DiscreteVariable state $state is duplicated.")
        end

        push!(stateRefs, state)
    end

    return stateRefs
end

"""
    DiscreteVariable(
        id::VariableId, cardinality::Int;
        label = defaultNodeLabel(id), states = 1:cardinality, probability = nothing
    )

Create a discrete variable node.

# Arguments

- `id`: Symbolic or integer ID used in factors.
- `cardinality`: Number of states.

# Keywords

- `label`: Human-readable label that can also be used for lookup.
- `states`: State references for entries of the finite state set.
- `probability`: Optional initial probability vector.

# Notes

If `label` is omitted, the variable ID is converted to a compact default label.
For example, `:x_1` becomes `"x1"`. Explicit labels are kept unchanged.

The state order defines the corresponding dimension order in connected factor tables.

If `probability` is provided, it defines the variable's initial belief for
future discrete inference objects. A unary [`DiscreteFactor`](@ref) with
`initialize = true` can override this initial belief for the connected variable.
The initial belief is used only to initialize belief-propagation messages. It is
not a replacement for a unary factor when the prior should be part of the model.
If `probability` is omitted, [`sumproduct`](@ref) uses a uniform initial message
and [`minsum`](@ref) uses the corresponding uniform initial cost.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1")
x2 = DiscreteVariable(:x2, 2; states = [:off, :on], probability = [0.8, 0.2])
```
"""
struct DiscreteVariable
    id::VariableId
    cardinality::Int
    states::Vector{StateRef}
    label::String
    probability::Union{Nothing, Vector{Float64}}

    function DiscreteVariable(
        id::VariableId,
        cardinality::Int;
        label::String = defaultNodeLabel(id),
        states = defaultStateRefs(cardinality),
        probability = nothing
    )
        if cardinality <= 0
            error("DiscreteVariable cardinality must be positive.")
        end

        if isempty(label)
            error("DiscreteVariable label must be nonempty.")
        end

        stateRefs = asStateRefs(states)

        if length(stateRefs) != cardinality
            error(
                "State dimension is not valid for DiscreteVariable $label. " *
                "Expected $cardinality, but got $(length(stateRefs))."
            )
        end

        probabilityVector =
            probability === nothing ?
            nothing :
            asDiscreteProbability(probability, cardinality, label)

        return new(id, cardinality, stateRefs, label, probabilityVector)
    end
end

"""
    DiscreteFactor(
        variable::Union{VariableRef, DiscreteVariable}..., table;
        label = "", initialize = false
    )

Create a discrete factor node.

# Arguments

- `variable...`: Connected variable references or DiscreteVariable nodes.
- `table`: Nonnegative factor table.

# Keywords

- `label`: Human-readable label. If omitted inside a graph, a label is assigned.
- `initialize`: Whether a unary factor initializes the connected variable belief.

# Notes

The table stores a nonnegative potential over the connected variables in the
same order as `variables`. The table does not need to be normalized.
Its dimensions are validated against the connected variable cardinalities when
the factor is added to a [`DiscreteFactorGraph`](@ref).

If `initialize = true`, the factor must be unary. Its table is used as the
initial belief for the connected variable when a compatible discrete inference
object is created or extended.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1")
x2 = DiscreteVariable(:x2, 2; label = "x2")

f1 = DiscreteFactor(x1, [0.8, 0.2]; label = "f1", initialize = true)
f2 = DiscreteFactor(x1, x2, [0.9 0.1; 0.2 0.8]; label = "f2")
```
"""
struct DiscreteFactor
    id::Int
    variables::Vector{VariableRef}
    table::Array{Float64}
    label::String
    initialize::Bool

    function DiscreteFactor(
        id::Int,
        variables::AbstractVector{<:VariableRef},
        table;
        label::String = "",
        initialize::Bool = false
    )
        if id < 0
            error("DiscreteFactor ID must be nonnegative.")
        end

        if isempty(variables)
            error("A DiscreteFactor must be connected to at least one variable.")
        end

        if isempty(label) && id > 0
            label = defaultFactorLabel(id)
        end

        return new(
            id,
            VariableRef[variables...],
            asDiscreteTable(table),
            label,
            initialize
        )
    end
end

variableRef(variable::DiscreteVariable) = variable.id

function discreteVariableRef(variable)
    error("DiscreteVariable references must be Int, Symbol, String, or DiscreteVariable.")
end

discreteVariableRef(variable::VariableRef) = variable
discreteVariableRef(variable::DiscreteVariable) = variable.id

function DiscreteFactor(args...; label::String = "", initialize::Bool = false)
    if length(args) < 2
        error("A DiscreteFactor must contain at least one variable and a table.")
    end

    variableArguments = args[1:end - 1]
    table = args[end]
    variables = VariableRef[]

    for variable in variableArguments
        push!(variables, discreteVariableRef(variable))
    end

    return DiscreteFactor(
        0,
        variables,
        table;
        label = label,
        initialize = initialize
    )
end

"""
    DiscreteFactorGraph(
        variables::AbstractVector{DiscreteVariable}, factors::AbstractVector{DiscreteFactor}
    )
    DiscreteFactorGraph()

Construct a discrete factor graph directly from variables and factors, or construct an empty one.

# Arguments

- `variables`: Discrete variable nodes.
- `factors`: Discrete factor nodes.

# Fields

- `variables`: Discrete variable nodes in internal order.
- `factors`: Discrete factor nodes in internal order.
- `referenceIndex`: Lookup from variable IDs and labels to internal indices.
- `edges`: Factor-variable edges.
- `factorEdges`: Edge IDs adjacent to each factor.
- `variableEdges`: Edge IDs adjacent to each variable.
- `topologyVersion`: Counter used to detect stale inference states.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = DiscreteFactorGraph([x1], [f1])
```
"""
mutable struct DiscreteFactorGraph <: AbstractFactorGraph
    variables::Vector{DiscreteVariable}
    factors::Vector{DiscreteFactor}
    referenceIndex::Dict{VariableRef, Int}
    edges::Vector{Edge}
    factorEdges::Vector{Vector{Int}}
    variableEdges::Vector{Vector{Int}}
    topologyVersion::Int
end

"""
    stateIndex(graph::DiscreteFactorGraph, variable::VariableRef, state::StateRef)

Resolve a discrete state reference to its one-based index within a variable.

# Arguments

- `graph`: Discrete factor graph used to resolve `variable`.
- `variable`: Variable ID or label.
- `state`: State reference.

# Returns

The one-based state index.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])

stateIndex(graph, :x1, :on)
```
"""
function stateIndex(
    graph::DiscreteFactorGraph,
    variableRef::VariableRef,
    state::StateRef
)
    variableIdx = variableIndex(graph, variableRef)

    return _stateIndex(graph.variables[variableIdx], state)
end

"""
    stateValue(graph::DiscreteFactorGraph, variable::VariableRef, index::Int)

Return the state reference stored at a one-based index.

# Arguments

- `graph`: Discrete factor graph used to resolve `variable`.
- `variable`: Variable ID or label.
- `index`: One-based state index.

# Returns

The state reference at `index`.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])

stateValue(graph, :x1, 2)
```
"""
function stateValue(
    graph::DiscreteFactorGraph,
    variableRef::VariableRef,
    index::Int
)
    variableIdx = variableIndex(graph, variableRef)

    return _stateValue(graph.variables[variableIdx], index)
end

function variableIndex(
    graph::DiscreteFactorGraph,
    variableRef::VariableRef
)
    return variableIndex(graph.referenceIndex, variableRef)
end

function variableDimension(
    graph::DiscreteFactorGraph,
    variableRef::VariableRef
)
    return variableDimension(graph.variables, graph.referenceIndex, variableRef)
end

function factorIndex(
    graph::DiscreteFactorGraph,
    factorRef::FactorRef
)
    return factorIndex(graph.factors, factorRef)
end

"""
    factorAxis(graph::DiscreteFactorGraph; factor::FactorRef, variable::VariableRef)

Return the table axis in `DiscreteFactor` corresponding to one variable.

# Arguments

- `graph`: Discrete factor graph.

# Keywords

- `factor`: Factor index or label.
- `variable`: Variable ID or label.

# Returns

The one-based table axis for the selected variable.

# Notes

Discrete factor table axes follow the same order as `DiscreteFactor.variables`.
This helper is useful when inspecting a multi-DiscreteVariable factor table.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
x2 = DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])
f1 = DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "f1")

graph = factorGraph([x1, x2], [f1])

factorAxis(graph; factor = "f1", variable = :x2)
```
"""
function factorAxis(
    graph::DiscreteFactorGraph;
    factor::FactorRef,
    variable::VariableRef
)
    factorData = graph.factors[factorIndex(graph, factor)]
    requestedIndex = variableIndex(graph.referenceIndex, variable)

    for (axis, variableRef) in pairs(factorData.variables)
        currentIndex = variableIndex(graph.referenceIndex, variableRef)

        if currentIndex == requestedIndex
            return axis
        end
    end

    error("Variable $variable is not connected to factor $(factorData.label).")
end

"""
    factorAxes(graph::DiscreteFactorGraph, factor::FactorRef)

Return table axes for all variables connected to a DiscreteFactor.

# Arguments

- `graph`: Discrete factor graph.
- `factor`: Factor index or label.

# Returns

The axes in the same order as `DiscreteFactor.variables`.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
x2 = DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])
f1 = DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "f1")

graph = factorGraph([x1, x2], [f1])

factorAxes(graph, "f1")
```
"""
function factorAxes(graph::DiscreteFactorGraph, factorRef::FactorRef)
    factorData = graph.factors[factorIndex(graph, factorRef)]

    return collect(eachindex(factorData.variables))
end

"""
    addVariable!(
        graph::DiscreteFactorGraph, variable::DiscreteVariable
    )
    addVariable!(
        graph::DiscreteFactorGraph, id::VariableId, cardinality::Int;
        label = defaultNodeLabel(id), states = 1:cardinality, probability = nothing
    )

Add a finite-state discrete variable to an existing discrete factor graph.

# Arguments

- `graph`: Discrete factor graph to mutate.
- `variable`: Discrete variable node to add.
- `id`: Variable ID.
- `cardinality`: Number of states.

# Keywords

- `label`: Variable label.
- `states`: State references for entries of the finite state set.
- `probability`: Optional initial probability vector.

# Returns

The added [`DiscreteVariable`](@ref).

# Notes

This graph-only method changes topology and makes existing inference objects stale. For
warm-start inference, use the method that also receives an inference object.

DiscreteVariable IDs and labels must remain unique.

# Example

```julia
graph = DiscreteFactorGraph()

addVariable!(graph, DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]))
addVariable!(graph, :x2, 2; label = "x2", states = [:low, :high])
```
"""
function addVariable!(graph::DiscreteFactorGraph, variable::DiscreteVariable)
    return addVariableNode!(graph, variable)
end

function addVariable!(
    graph::DiscreteFactorGraph,
    id::VariableId,
    cardinality::Int;
    label::String = defaultNodeLabel(id),
    states = defaultStateRefs(cardinality),
    probability = nothing
)
    return addVariable!(
        graph,
        DiscreteVariable(
            id,
            cardinality;
            label = label,
            states = states,
            probability = probability
        )
    )
end

"""
    addFactor!(
        graph::DiscreteFactorGraph, factorData::DiscreteFactor
    )
    addFactor!(
        graph::DiscreteFactorGraph, variable::Union{VariableRef, DiscreteVariable}..., table;
        label = "", initialize = false
    )

Add a DiscreteFactor node to an existing discrete factor graph.

# Arguments

- `graph`: Discrete factor graph to mutate.
- `factorData`: Factor node to add.
- `variable...`: Connected variable references or DiscreteVariable nodes.
- `table`: Nonnegative factor table.

# Keywords

- `label`: Factor label.
- `initialize`: Whether a unary factor initializes the connected variable belief.

# Returns

The added [`DiscreteFactor`](@ref).

# Notes

This graph-only method changes topology and makes existing inference objects stale. For
warm-start inference, use the method that also receives an inference object.

If `label` is omitted, the DiscreteFactor is named `fN`, where `N` is the
assigned DiscreteFactor ID. Explicit labels are kept unchanged.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])

graph = DiscreteFactorGraph([x1], DiscreteFactor[])

addFactor!(graph, DiscreteFactor(:x1, [0.8, 0.2]; label = "f1"))
```
"""
function addFactor!(graph::DiscreteFactorGraph, factorData::DiscreteFactor)
    factorId = length(graph.factors) + 1
    factorLabel = isempty(factorData.label) ? defaultFactorLabel(factorId) : factorData.label

    if hasFactorLabel(graph, factorLabel)
        error("Factor label $factorLabel is already defined.")
    end

    assertNoInitializingUnaryFactorConflict(graph, factorData)

    added = DiscreteFactor(
        factorId,
        factorData.variables,
        factorData.table;
        label = factorLabel,
        initialize = factorData.initialize
    )

    validateFactor(
        added,
        graph.variables,
        graph.referenceIndex
    )

    newEdges = Edge[]
    nextEdgeId = length(graph.edges) + 1

    for variableRef in added.variables
        variableIdx = variableIndex(graph.referenceIndex, variableRef)
        push!(newEdges, Edge(nextEdgeId, factorId, variableIdx))
        nextEdgeId += 1
    end

    push!(graph.factors, added)
    push!(graph.factorEdges, Int[])

    for edge in newEdges
        push!(graph.edges, edge)
        push!(graph.factorEdges[factorId], edge.id)
        push!(graph.variableEdges[edge.variableIndex], edge.id)
    end

    bumpTopologyVersion!(graph)

    return added
end

function addFactor!(
    graph::DiscreteFactorGraph,
    args...;
    label::String = "",
    initialize::Bool = false
)
    return addFactor!(graph, DiscreteFactor(args...; label = label, initialize = initialize))
end

function validateUpdatedFactorDimensions(current::DiscreteFactor, updated::DiscreteFactor)
    if size(updated.table) != size(current.table)
        error(
            "Updated factor table size must remain $(size(current.table)), " *
            "but got $(size(updated.table))."
        )
    end

    return nothing
end

function validateUpdatedInitializingUnaryFactors(
    graph::DiscreteFactorGraph,
    factorIdx::Int,
    updated::DiscreteFactor
)
    factors = copy(graph.factors)
    factors[factorIdx] = updated

    validateInitializingUnaryFactors(factors, graph.referenceIndex)

    return nothing
end

"""
    updateFactor!(
        graph::DiscreteFactorGraph, factor::FactorRef;
        table = nothing, initialize = nothing
    )

Update a DiscreteFactor table without changing graph topology or table dimensions.

# Arguments

- `graph`: Discrete factor graph to mutate.
- `factor`: Factor index or label.

# Keywords

- `table`: Replacement table.
- `initialize`: Replacement initialization flag; omitted keeps the current flag.

# Returns

The updated [`DiscreteFactor`](@ref).

# Notes

The connected variables stay fixed. The replacement table must have the same
dimensions as the current table. This method does not change graph topology, so
existing inference objects are not stale, but use the inference-aware
[`updateFactor!`](@ref) method when a running inference state should also update
an initializing unary factor.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = DiscreteFactorGraph([x1], [f1])

updateFactor!(graph, "f1"; table = [0.6, 0.4])
```
"""
function updateFactor!(
    graph::DiscreteFactorGraph,
    factorRef::FactorRef;
    table = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    factorIdx = factorIndex(graph, factorRef)
    current = graph.factors[factorIdx]
    nextTable = table === nothing ? current.table : table
    nextInitialize = initialize === nothing ? current.initialize : initialize

    updated = DiscreteFactor(
        current.id,
        current.variables,
        nextTable;
        label = current.label,
        initialize = nextInitialize
    )

    validateUpdatedFactorDimensions(current, updated)

    validateFactor(
        updated,
        graph.variables,
        graph.referenceIndex
    )

    validateUpdatedInitializingUnaryFactors(graph, factorIdx, updated)

    graph.factors[factorIdx] = updated

    return updated
end

function updateFactor!(
    graph::DiscreteFactorGraph;
    factor::FactorRef,
    table = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    return updateFactor!(graph, factor; table = table, initialize = initialize)
end

function DiscreteFactorGraph()
    return DiscreteFactorGraph(DiscreteVariable[], DiscreteFactor[])
end

function DiscreteFactorGraph(
    variables::AbstractVector{DiscreteVariable},
    factors::AbstractVector{DiscreteFactor}
)
    variableList = Vector{DiscreteVariable}(variables)

    validateVariablesIfAny(variableList)

    referenceIndex = buildRefIndex(variableList)
    factorList = assignFactorIds(factors)

    validateFactorsIfAny(factorList, variableList, referenceIndex)

    edges, factorEdges, variableEdges = buildEdges(variableList, factorList, referenceIndex)

    return DiscreteFactorGraph(
        variableList,
        factorList,
        referenceIndex,
        edges,
        factorEdges,
        variableEdges,
        0
    )
end

function factorGraph(
    variables::AbstractVector{DiscreteVariable},
    factors::AbstractVector{DiscreteFactor};
    root::Union{Nothing, VariableRef} = nothing
)
    graph = DiscreteFactorGraph(variables, factors)

    if root === nothing
        return graph
    end

    return treeFactorGraph(graph; root = root)
end

function variableIndex(
    referenceIndex::Dict{VariableRef, Int},
    variable::DiscreteVariable
)
    return variableIndex(referenceIndex, variable.id)
end

function assignFactorIds(factors::AbstractVector{DiscreteFactor})
    factorList = DiscreteFactor[]

    for (index, factorData) in pairs(factors)
        factorLabel = isempty(factorData.label) ? defaultFactorLabel(index) : factorData.label

        push!(
            factorList,
            DiscreteFactor(
                index,
                factorData.variables,
                factorData.table;
                label = factorLabel,
                initialize = factorData.initialize
            )
        )
    end

    return factorList
end

function variableDimension(
    variables::Vector{DiscreteVariable},
    referenceIndex::Dict{VariableRef, Int},
    variableRef::VariableRef
)
    index = variableIndex(referenceIndex, variableRef)

    return variables[index].cardinality
end

function _stateIndex(variable::DiscreteVariable, state::StateRef)
    for (index, stateRef) in pairs(variable.states)
        if stateRef == state
            return index
        end
    end

    error("State reference $state is not defined for DiscreteVariable $(variable.label).")
end

function _stateValue(variable::DiscreteVariable, index::Int)
    if index < 1 || index > variable.cardinality
        error(
            "State index $index is not valid for DiscreteVariable $(variable.label). " *
            "Expected 1:$(variable.cardinality)."
        )
    end

    return variable.states[index]
end

function validateFactor(
    factorData::DiscreteFactor,
    variables::Vector{DiscreteVariable},
    referenceIndex::Dict{VariableRef, Int}
)
    resolvedIndices = Int[]

    for variableRef in factorData.variables
        push!(resolvedIndices, variableIndex(referenceIndex, variableRef))
    end

    if length(unique(resolvedIndices)) != length(resolvedIndices)
        error("Factor $(factorData.label) contains duplicate variable references.")
    end

    expectedSize = Tuple(variables[index].cardinality for index in resolvedIndices)

    if size(factorData.table) != expectedSize
        error(
            "Factor table size is not valid for factor $(factorData.label). " *
            "Expected $expectedSize, but got $(size(factorData.table))."
        )
    end

    if factorData.initialize && length(factorData.variables) != 1
        error("Only unary factors can be used for variable initialization.")
    end

    return nothing
end

function validateInitializingUnaryFactors(
    factors::Vector{DiscreteFactor},
    referenceIndex::Dict{VariableRef, Int}
)
    initializedVariables = Dict{Int, String}()

    for factorData in factors
        if !factorData.initialize
            continue
        end

        variableIdx = variableIndex(referenceIndex, only(factorData.variables))

        if haskey(initializedVariables, variableIdx)
            error(
                "Variable $(factorData.variables[1]) has more than one " *
                "initializing unary factor: $(initializedVariables[variableIdx]) " *
                "and $(factorData.label)."
            )
        end

        initializedVariables[variableIdx] = factorData.label
    end

    return nothing
end

function validateFactors(
    factors::Vector{DiscreteFactor},
    variables::Vector{DiscreteVariable},
    referenceIndex::Dict{VariableRef, Int}
)
    if isempty(factors)
        error("At least one factor node must be defined.")
    end

    labels = [factorData.label for factorData in factors]

    if length(unique(labels)) != length(labels)
        error("Factor labels must be unique.")
    end

    for factorData in factors
        validateFactor(factorData, variables, referenceIndex)
    end

    validateInitializingUnaryFactors(factors, referenceIndex)

    return nothing
end

function initializingUnaryFactor(
    factors::Vector{DiscreteFactor},
    variables::Vector{DiscreteVariable},
    referenceIndex::Dict{VariableRef, Int},
    variableIndexValue::Int
)
    found = nothing

    for factorData in factors
        if !factorData.initialize
            continue
        end

        factorVariableIndex = variableIndex(
            referenceIndex,
            only(factorData.variables)
        )

        if factorVariableIndex != variableIndexValue
            continue
        end

        if found !== nothing
            error(
                "Variable $(variables[variableIndexValue].label) has " *
                "more than one initializing unary factor."
            )
        end

        found = factorData
    end

    return found
end

function assertNoInitializingUnaryFactorConflict(
    graph::DiscreteFactorGraph,
    factorData::DiscreteFactor
)
    if !factorData.initialize
        return nothing
    end

    if length(factorData.variables) != 1
        error("Only unary factors can be used for variable initialization.")
    end

    variableRef = only(factorData.variables)
    variableIdx = variableIndex(graph.referenceIndex, variableRef)
    existing = initializingUnaryFactor(graph, variableIdx)

    if existing !== nothing
        error(
            "Variable $variableRef already has initializing unary factor " *
            "$(existing.label)."
        )
    end

    return nothing
end

function treeFactorGraph(
    graph::DiscreteFactorGraph;
    root::Union{Nothing, VariableRef} = nothing
)
    return treeFactorGraphImpl(graph; root = root)
end

function printModel(graph::DiscreteFactorGraph)
    println("Discrete factor graph")
    println("Number of DiscreteVariable nodes: ", length(graph.variables))
    println("Number of DiscreteFactor nodes: ", length(graph.factors))
    println()

    println("DiscreteVariable nodes:")
    for variable in graph.variables
        println(
            "  label = ",
            variable.label,
            ", id = ",
            variable.id,
            ", cardinality = ",
            variable.cardinality,
            ", states = ",
            formattedRefs(variable.states)
        )
    end

    println()
    println("DiscreteFactor nodes:")
    for factorData in graph.factors
        println("  label = ", factorData.label, ", id = ", factorData.id)
        println("    variables = ", formattedRefs(factorData.variables))
        println("    table size = ", size(factorData.table))
        println("    initialize = ", factorData.initialize)
    end

    return nothing
end

function printModel(tree::TreeFactorGraph{DiscreteFactorGraph})
    return printModel(tree.graph)
end
