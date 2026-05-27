function asCoefficientMatrix(data)
    if data isa Number
        return reshape([Float64(data)], 1, 1)
    elseif data isa AbstractMatrix
        return Matrix{Float64}(data)
    elseif data isa AbstractVector && length(data) == 1
        return reshape([Float64(only(data))], 1, 1)
    else
        error("Coefficient must be a scalar, a one-element vector, or a matrix.")
    end
end

function asCovarianceMatrix(data)
    if data isa Number
        return reshape([Float64(data)], 1, 1)
    elseif data isa AbstractMatrix
        return Matrix{Float64}(data)
    elseif data isa AbstractVector
        vector = Vector{Float64}(data)

        if length(vector) == 1
            return reshape(vector, 1, 1)
        end

        return Diagonal(vector) |> Matrix{Float64}
    else
        error("Covariance must be a scalar, a vector, or a matrix.")
    end
end

function defaultComponentRefs(dimension::Int)
    return ComponentRef[1:dimension...]
end

function asComponentRefs(components)
    if !(components isa AbstractVector)
        error("GaussianVariable components must be a vector.")
    end

    if isempty(components)
        error("GaussianVariable components must be nonempty.")
    end

    componentRefs = ComponentRef[]

    for component in components
        if !(component isa ComponentRef)
            error("GaussianVariable components must be Int, Symbol, or String.")
        end

        if component isa String && isempty(component)
            error("GaussianVariable component must be nonempty.")
        end

        if component in componentRefs
            error("GaussianVariable component $component is duplicated.")
        end

        push!(componentRefs, component)
    end

    return componentRefs
end

"""
    GaussianVariable(
        id::VariableId, dimension::Int;
        label = string(id), components = 1:dimension, mean = nothing, covariance = nothing
    )

Create a GaussianVariable node.

# Arguments

- `id`: Symbolic or integer identifier used in factors.
- `dimension`: State dimension of the variable.

# Keywords

- `label`: Human-readable label that can also be used for lookup.
- `components`: Names for entries of the continuous state vector.
- `mean`: Optional initial mean for Gaussian inference.
- `covariance`: Optional initial covariance for Gaussian inference.

# Notes

If `mean` and `covariance` are provided, they define the GaussianVariable initial belief
used by new inference objects. If they are omitted, [`moment`](@ref), [`canonical`](@ref),
and Gaussian [`minsum`](@ref) use their own default initial belief arguments instead. This
belief is used only to initialize Gaussian belief propagation messages. Scalar values are
accepted for one-dimensional variables.

# Example

```julia
x1 = GaussianVariable(:x1, 1; label = "x1", mean = 0.0, covariance = 1.0)
x2 = GaussianVariable(:x2, 2; mean = [0.0, 0.0], covariance = [1.0, 1.0])
x3 = GaussianVariable(:x3, 2; label = "x3", components = [:position, :velocity])
```
"""
struct GaussianVariable
    id::VariableId
    dimension::Int
    components::Vector{ComponentRef}
    label::String
    mean::Union{Nothing, Vector{Float64}}
    covariance::Union{Nothing, Matrix{Float64}}

    function GaussianVariable(
        id::VariableId,
        dimension::Int;
        label::String = string(id),
        components = defaultComponentRefs(dimension),
        mean = nothing,
        covariance = nothing
    )
        if dimension <= 0
            error("GaussianVariable dimension must be positive.")
        end

        if isempty(label)
            error("GaussianVariable label must be nonempty.")
        end

        componentRefs = asComponentRefs(components)

        if length(componentRefs) != dimension
            error(
                "Component dimension is not valid for GaussianVariable $label. " *
                "Expected $dimension, but got $(length(componentRefs))."
            )
        end

        meanVector = nothing

        if mean !== nothing
            meanVector = asVector(mean)

            if length(meanVector) != dimension
                error(
                    "Mean dimension is not valid for GaussianVariable $label. " *
                    "Expected $dimension, but got $(length(meanVector))."
                )
            end
        end

        covarianceMatrix = nothing

        if covariance !== nothing
            covarianceMatrix = symmetricPart(asCovarianceMatrix(covariance))

            if size(covarianceMatrix, 1) != dimension
                error(
                    "Covariance row dimension is not valid for GaussianVariable $label. " *
                    "Expected $dimension, but got $(size(covarianceMatrix, 1))."
                )
            end

            if size(covarianceMatrix, 2) != dimension
                error(
                    "Covariance column dimension is not valid for GaussianVariable $label. " *
                    "Expected $dimension, but got $(size(covarianceMatrix, 2))."
                )
            end

            if !isposdef(Symmetric(covarianceMatrix))
                error("Covariance must be positive definite for GaussianVariable $label.")
            end
        end

        return new(id, dimension, componentRefs, label, meanVector, covarianceMatrix)
    end
end

"""
    GaussianFactor(
        variable::Union{VariableRef, GaussianVariable}..., mean, coefficient, covariance;
        label = "", initialize = false
    )

Create a linear Gaussian factor node.

# Arguments

- `variable...`: Connected variable references or GaussianVariable nodes.
- `mean`: Observed or target value.
- `coefficient`: Linear coefficient matrix.
- `covariance`: Measurement or factor covariance.

# Keywords

- `label`: Human-readable label. If omitted inside a graph, a label is assigned.
- `initialize`: Whether a unary factor initializes the connected variable belief.

# Notes

The GaussianFactor represents a linear Gaussian relation
`mean = coefficient * x + noise`, where `mean` is the GaussianFactor's observed or
target value and the noise covariance is `covariance`. The columns of
`coefficient` are ordered by the variables listed in `variables`.

The constructor takes one or more GaussianVariable references or nodes followed by
`mean`, `coefficient`, and `covariance`.

Scalars and one-element vectors are accepted for one-dimensional factors. Coefficients for
higher-dimensional factors must be given as matrices.

If `initialize = true`, the GaussianFactor must be unary. Its Gaussian information is
used as the initial belief for the connected GaussianVariable when an inference object
is created or extended.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
x2 = GaussianVariable(:x2, 1)

f1 = GaussianFactor(:x1, 0.25, 1.0, 0.15; label = "f1", initialize = true)
f2 = GaussianFactor(x1, x2, [0.4, -0.2], [1.0 0.5; -0.4 1.2], [0.3, 0.25])
```
"""
struct GaussianFactor
    id::Int
    variables::Vector{VariableRef}
    mean::Vector{Float64}
    coefficient::Matrix{Float64}
    covariance::Matrix{Float64}
    label::String
    initialize::Bool

    function GaussianFactor(
        id::Int,
        variables::AbstractVector{<:VariableRef},
        mean,
        coefficient,
        covariance;
        label::String = "",
        initialize::Bool = false
    )
        if id < 0
            error("GaussianFactor id must be nonnegative.")
        end

        if isempty(variables)
            error("A GaussianFactor must be connected to at least one variable.")
        end

        if isempty(label) && id > 0
            label = "f_$id"
        end

        covarianceMatrix = symmetricPart(asCovarianceMatrix(covariance))

        return new(
            id,
            VariableRef[variables...],
            asVector(mean),
            asCoefficientMatrix(coefficient),
            covarianceMatrix,
            label,
            initialize
        )
    end
end

variableRef(variable::GaussianVariable) = variable.id

function gaussianVariableRef(variable)
    error("GaussianVariable references must be Int, Symbol, String, or GaussianVariable.")
end

gaussianVariableRef(variable::VariableRef) = variable
gaussianVariableRef(variable::GaussianVariable) = variable.id

function GaussianFactor(args...; label::String = "", initialize::Bool = false)
    if length(args) < 4
        error(
            "A GaussianFactor must contain at least one GaussianVariable and the " *
            "arguments mean, coefficient, and covariance."
        )
    end

    variableArguments = args[1:end - 3]
    mean = args[end - 2]
    coefficient = args[end - 1]
    covariance = args[end]

    variables = VariableRef[]

    for variable in variableArguments
        push!(variables, gaussianVariableRef(variable))
    end

    return GaussianFactor(
        0,
        variables,
        mean,
        coefficient,
        covariance;
        label = label,
        initialize = initialize
    )
end

"""
    GaussianFactorGraph(
        variables::AbstractVector{GaussianVariable}, factors::AbstractVector{GaussianFactor}
    )
    GaussianFactorGraph()

Construct a Gaussian factor graph from variables and factors, or construct an empty one.

# Arguments

- `variables`: Gaussian variable nodes.
- `factors`: Gaussian factor nodes.

# Fields

- `variables`: Gaussian variable nodes in internal order.
- `factors`: Gaussian factor nodes in internal order.
- `referenceIndex`: Lookup from variable ids and labels to internal indices.
- `edges`: Factor-variable edges.
- `factorEdges`: Edge ids adjacent to each factor.
- `variableEdges`: Edge ids adjacent to each variable.
- `topologyVersion`: Counter used to detect stale inference states.

# Notes

The graph stores edges and adjacency lists used by message passing. New factors
may be added with [`addFactor!`](@ref), and existing GaussianFactor means and
coefficients and covariances may be updated with [`updateFactor!`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
x2 = GaussianVariable(:x2, 1)

f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")
f2 = GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.2; label = "f2")

graph = GaussianFactorGraph([x1, x2], [f1, f2])
```
"""
mutable struct GaussianFactorGraph <: AbstractFactorGraph
    variables::Vector{GaussianVariable}
    factors::Vector{GaussianFactor}
    referenceIndex::Dict{VariableRef, Int}
    edges::Vector{Edge}
    factorEdges::Vector{Vector{Int}}
    variableEdges::Vector{Vector{Int}}
    topologyVersion::Int
end

"""
    componentIndex(graph::GaussianFactorGraph, variable::VariableRef, component::ComponentRef)

Resolve a Gaussian component reference to its one-based index within a variable.

# Arguments

- `graph`: Gaussian factor graph used to resolve `variable`.
- `variable`: Variable id or label.
- `component`: Component reference.

# Returns

The one-based component index.

# Example

```julia
x1 = GaussianVariable(:x1, 2; label = "state", components = [:position, :velocity])
f1 = GaussianFactor(:x1, [0.0, 0.0], [1.0 0.0; 0.0 1.0], [0.1 0.0; 0.0 0.1])

graph = factorGraph([x1], [f1])

componentIndex(graph, :x1, :velocity)
```
"""
function componentIndex(
    graph::GaussianFactorGraph,
    variableRef::VariableRef,
    component::ComponentRef
)
    variableIdx = variableIndex(graph, variableRef)

    return _componentIndex(graph.variables[variableIdx], component)
end

"""
    componentValue(graph::GaussianFactorGraph, variable::VariableRef, index::Int)

Return the component reference stored at a one-based index.

# Arguments

- `graph`: Gaussian factor graph used to resolve `variable`.
- `variable`: Variable id or label.
- `index`: One-based component index.

# Returns

The component reference at `index`.

# Example

```julia
x1 = GaussianVariable(:x1, 2; label = "state", components = [:position, :velocity])
f1 = GaussianFactor(:x1, [0.0, 0.0], [1.0 0.0; 0.0 1.0], [0.1 0.0; 0.0 0.1])

graph = factorGraph([x1], [f1])

componentValue(graph, :x1, 1)
```
"""
function componentValue(
    graph::GaussianFactorGraph,
    variableRef::VariableRef,
    index::Int
)
    variableIdx = variableIndex(graph, variableRef)

    return _componentValue(graph.variables[variableIdx], index)
end

"""
    variableIndex(graph::AbstractFactorGraph, variable::VariableRef)

Resolve a variable reference to its internal one-based variable index.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `variable`: Variable id or label.

# Returns

The one-based internal variable index.

# Example

```julia
x1 = GaussianVariable(:x1, 1; label = "state")
graph = factorGraph([x1], GaussianFactor[])

variableIndex(graph, "state")
```
"""
function variableIndex(
    graph::GaussianFactorGraph,
    variableRef::VariableRef
)
    return variableIndex(graph.referenceIndex, variableRef)
end

"""
    variableDimension(graph::AbstractFactorGraph, variable::VariableRef)

Return the dimension of a variable.

For Gaussian graphs this is the continuous state-vector dimension. For discrete
graphs this is the variable cardinality.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `variable`: Variable id or label.

# Returns

The variable dimension.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; states = [:off, :on])
graph = factorGraph([x1], DiscreteFactor[])

variableDimension(graph, :x1)
```
"""
function variableDimension(
    graph::GaussianFactorGraph,
    variableRef::VariableRef
)
    return variableDimension(graph.variables, graph.referenceIndex, variableRef)
end

function coefficientBlock(
    graph::GaussianFactorGraph,
    factorData::GaussianFactor,
    variableRef::VariableRef
)
    return coefficientBlock(graph.variables, graph.referenceIndex, factorData, variableRef)
end

function coefficientBlock(
    graph::GaussianFactorGraph;
    factor::FactorRef,
    variable::VariableRef
)
    factorData = graph.factors[factorIndex(graph, factor)]

    return coefficientBlock(graph, factorData, variable)
end

function coefficientBlocks(
    graph::GaussianFactorGraph,
    factorRef::FactorRef
)
    factorData = graph.factors[factorIndex(graph, factorRef)]

    return coefficientBlocks(graph, factorData)
end

function coefficientBlocks(graph::GaussianFactorGraph, factorData::GaussianFactor)
    return coefficientBlocks(graph.variables, graph.referenceIndex, factorData)
end

"""
    factorIndex(graph::AbstractFactorGraph, factor::FactorRef)

Resolve a factor reference to its internal one-based factor index.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `factor`: Factor index or label.

# Returns

The one-based internal factor index.

# Example

```julia
x1 = DiscreteVariable(:x1, 2)
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])

factorIndex(graph, "f1")
```
"""
function factorIndex(
    graph::GaussianFactorGraph,
    factorRef::FactorRef
)
    return factorIndex(graph.factors, factorRef)
end

"""
    addVariable!(
        graph::GaussianFactorGraph, variable::GaussianVariable
    )
    addVariable!(
        graph::GaussianFactorGraph, id::VariableId, dimension::Int;
        label = string(id), components = 1:dimension, mean = nothing, covariance = nothing
    )

Add a continuous Gaussian variable to an existing Gaussian factor graph.

# Arguments

- `graph`: Gaussian factor graph to mutate.
- `variable`: Gaussian variable node to add.
- `id`: Variable id.
- `dimension`: Dimension of the continuous state vector.

# Keywords

- `label`: Variable label.
- `components`: Component references for entries of the continuous state vector.
- `mean`: Optional initial mean.
- `covariance`: Optional initial covariance.

# Returns

The added [`GaussianVariable`](@ref).

# Notes

This is intended for incremental graph construction with [`factorGraph`](@ref).
Calling this graph-only method changes the graph topology and makes existing
inference objects for this graph stale. For warm-start inference, use the
method that also receives an inference object.

GaussianVariable ids and labels must remain unique.

# Example

```julia
graph = GaussianFactorGraph()

addVariable!(graph, GaussianVariable(:x1, 1; label = "x1"))
addVariable!(graph, :x2, 2; label = "x2", components = [:position, :velocity])
```
"""
function addVariable!(graph::GaussianFactorGraph, variable::GaussianVariable)
    return addVariableNode!(graph, variable)
end

function addVariable!(
    graph::GaussianFactorGraph,
    id::VariableId,
    dimension::Int;
    label::String = string(id),
    components = defaultComponentRefs(dimension),
    mean = nothing,
    covariance = nothing
)
    return addVariable!(
        graph,
        GaussianVariable(
            id,
            dimension;
            label = label,
            components = components,
            mean = mean,
            covariance = covariance
        )
    )
end

"""
    addFactor!(
        graph::GaussianFactorGraph, factor::GaussianFactor
    )
    addFactor!(
        graph::GaussianFactorGraph,
        variable::Union{VariableRef, GaussianVariable}..., mean, coefficient, covariance;
        label = "", initialize = false
    )

Add a GaussianFactor node to an existing Gaussian factor graph.

# Arguments

- `graph`: Gaussian factor graph to mutate.
- `factor`: Factor node to add.
- `variable...`: Connected variable references or GaussianVariable nodes for constructor-style
  insertion.
- `mean`: Observed or target value.
- `coefficient`: Linear coefficient matrix.
- `covariance`: Factor covariance.

# Keywords

- `label`: Factor label.
- `initialize`: Whether a unary factor initializes the connected variable belief.

# Returns

The added [`GaussianFactor`](@ref).

# Notes

The new GaussianFactor is validated against the graph's existing variables, assigned
the next GaussianFactor id, inserted into `graph.factors`, and connected by new
edges. Existing variables are not modified, but their adjacency lists receive
the new edge ids.

If `label` is omitted, the GaussianFactor is named `f_N`, where `N` is the assigned
GaussianFactor id. GaussianFactor labels must remain unique.

Calling this graph-only method does not resize existing inference objects.
Create a fresh [`moment`](@ref) or [`canonical`](@ref) inference object before
running message passing on the updated graph, or call
`addFactor!(graph, inference, factorData)` to extend an existing inference object
as a warm start.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
x2 = GaussianVariable(:x2, 1)

graph = GaussianFactorGraph([x1, x2], GaussianFactor[])

addFactor!(graph, GaussianFactor(:x1, 0.5, 1.0, 0.1))
addFactor!(graph, :x1, :x2, [0.4, -0.2], [1.0 0.5; -0.4 1.2], [0.3, 0.25])
```
"""
function addFactor!(graph::GaussianFactorGraph, factorData::GaussianFactor)
    factorId = length(graph.factors) + 1
    factorLabel = isempty(factorData.label) ? "f_$factorId" : factorData.label

    if hasFactorLabel(graph, factorLabel)
        error("Factor label $factorLabel is already defined.")
    end

    assertNoInitializingUnaryFactorConflict(graph, factorData)

    added = GaussianFactor(
        factorId,
        factorData.variables,
        factorData.mean,
        factorData.coefficient,
        factorData.covariance;
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
    graph::GaussianFactorGraph,
    args...;
    label::String = "",
    initialize::Bool = false
)
    return addFactor!(graph, GaussianFactor(args...; label = label, initialize = initialize))
end

function validateUpdatedFactorDimensions(current::GaussianFactor, updated::GaussianFactor)
    if length(updated.mean) != length(current.mean)
        error(
            "Updated factor mean dimension must remain $(length(current.mean)), " *
            "but got $(length(updated.mean))."
        )
    end

    if size(updated.coefficient) != size(current.coefficient)
        error(
            "Updated factor coefficient size must remain $(size(current.coefficient)), " *
            "but got $(size(updated.coefficient))."
        )
    end

    if size(updated.covariance) != size(current.covariance)
        error(
            "Updated factor covariance size must remain $(size(current.covariance)), " *
            "but got $(size(updated.covariance))."
        )
    end

    return nothing
end

function validateUpdatedInitializingUnaryFactors(
    graph::GaussianFactorGraph,
    factorIdx::Int,
    updated::GaussianFactor
)
    factors = copy(graph.factors)
    factors[factorIdx] = updated

    validateInitializingUnaryFactors(factors, graph.referenceIndex)

    return nothing
end

"""
    updateFactor!(
        graph::GaussianFactorGraph, factor::FactorRef;
        mean = nothing, coefficient = nothing, covariance = nothing, initialize = nothing
    )

Update a GaussianFactor node's mean, coefficient, and/or covariance without changing
graph topology or GaussianFactor dimensions.

# Arguments

- `graph`: Gaussian factor graph to mutate.
- `factor`: Factor index or label.

# Keywords

- `mean`: Replacement observed or target value.
- `coefficient`: Replacement coefficient matrix.
- `covariance`: Replacement covariance.
- `initialize`: Replacement initialization flag; omitted keeps the current flag.

# Returns

The updated [`GaussianFactor`](@ref).

# Notes

Use `mean` for the observed or target value stored in the GaussianFactor. Use
`coefficient` to replace the GaussianFactor coefficient matrix.

The update keeps the existing connected variables, id, and label. If
`initialize` is omitted, the existing initialization flag is kept. The updated
mean length, coefficient matrix size, and covariance matrix size must match the
old GaussianFactor. Dimensions, positive definiteness, and initializing unary GaussianFactor
conflicts are validated before the graph is modified, so a failed update leaves the old
GaussianFactor unchanged.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.25, 1.0, 0.15; label = "f1", initialize = true)

graph = GaussianFactorGraph([x1], [f1])

updateFactor!(graph, "f1"; mean = 0.6, covariance = 0.2)
```
"""
function updateFactor!(
    graph::GaussianFactorGraph,
    factorRef::FactorRef;
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    factorIdx = factorIndex(graph, factorRef)
    current = graph.factors[factorIdx]

    nextMean = mean === nothing ? current.mean : mean
    nextCoefficient = coefficient === nothing ? current.coefficient : coefficient
    nextCovariance = covariance === nothing ? current.covariance : covariance
    nextInitialize = initialize === nothing ? current.initialize : initialize

    updated = GaussianFactor(
        current.id,
        current.variables,
        nextMean,
        nextCoefficient,
        nextCovariance;
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
    graph::GaussianFactorGraph;
    factor::FactorRef,
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    return updateFactor!(
        graph,
        factor;
        mean = mean,
        coefficient = coefficient,
        covariance = covariance,
        initialize = initialize
    )
end

function GaussianFactorGraph()
    return GaussianFactorGraph(GaussianVariable[], GaussianFactor[])
end

function GaussianFactorGraph(
    variables::AbstractVector{GaussianVariable},
    factors::AbstractVector{GaussianFactor}
)
    variableList = Vector{GaussianVariable}(variables)

    validateVariablesIfAny(variableList)

    referenceIndex = buildRefIndex(variableList)
    factorList = assignFactorIds(factors)

    validateFactorsIfAny(factorList, variableList, referenceIndex)

    edges, factorEdges, variableEdges = buildEdges(variableList, factorList, referenceIndex)

    return GaussianFactorGraph(
        variableList,
        factorList,
        referenceIndex,
        edges,
        factorEdges,
        variableEdges,
        0
    )
end

"""
    factorGraph(
        variables::AbstractVector{GaussianVariable}, factors::AbstractVector{GaussianFactor};
        root = nothing
    )
    factorGraph(
        variables::AbstractVector{DiscreteVariable}, factors::AbstractVector{DiscreteFactor};
        root = nothing
    )

Construct a Gaussian or discrete factor graph from variable and factor nodes.

# Arguments

- `variables`: Variable nodes.
- `factors`: Factor nodes.

# Keywords

- `root`: Optional root variable for returning a [`TreeFactorGraph`](@ref).

# Returns

A [`GaussianFactorGraph`](@ref) or [`DiscreteFactorGraph`](@ref), or a
[`TreeFactorGraph`](@ref) view when `root` is provided.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; states = [:off, :on])
f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")

graph = factorGraph([x1], [f1])
```
"""
function factorGraph(
    variables::AbstractVector{GaussianVariable},
    factors::AbstractVector{GaussianFactor};
    root::Union{Nothing, VariableRef} = nothing
)
    graph = GaussianFactorGraph(variables, factors)

    if root === nothing
        return graph
    end

    return treeFactorGraph(graph; root = root)
end

function variableIndex(
    referenceIndex::Dict{VariableRef, Int},
    variable::GaussianVariable
)
    return variableIndex(referenceIndex, variable.id)
end

function assignFactorIds(factors::AbstractVector{GaussianFactor})
    factorList = GaussianFactor[]

    for (index, factorData) in pairs(factors)
        factorLabel = isempty(factorData.label) ? "f_$index" : factorData.label

        push!(
            factorList,
            GaussianFactor(
                index,
                factorData.variables,
                factorData.mean,
                factorData.coefficient,
                factorData.covariance;
                label = factorLabel,
                initialize = factorData.initialize
            )
        )
    end

    return factorList
end

function variableDimension(
    variables::Vector{GaussianVariable},
    referenceIndex::Dict{VariableRef, Int},
    variable::VariableRef
)
    index = variableIndex(referenceIndex, variable)

    return variables[index].dimension
end

function _componentIndex(variable::GaussianVariable, component::ComponentRef)
    for (index, componentRef) in pairs(variable.components)
        if componentRef == component
            return index
        end
    end

    error(
        "Component reference $component is not defined for " *
        "GaussianVariable $(variable.label)."
    )
end

function _componentValue(variable::GaussianVariable, index::Int)
    if index < 1 || index > variable.dimension
        error(
            "Component index $index is not valid for GaussianVariable $(variable.label). " *
            "Expected 1:$(variable.dimension)."
        )
    end

    return variable.components[index]
end

function validateFactor(
    factorData::GaussianFactor,
    variables::Vector{GaussianVariable},
    referenceIndex::Dict{VariableRef, Int}
)
    resolvedIndices = Int[]

    for variableRef in factorData.variables
        push!(resolvedIndices, variableIndex(referenceIndex, variableRef))
    end

    if length(unique(resolvedIndices)) != length(resolvedIndices)
        error("Factor $(factorData.label) contains duplicate variable references.")
    end

    meanDimension = length(factorData.mean)

    if meanDimension <= 0
        error("Factor mean vector must be nonempty for factor $(factorData.label).")
    end

    if size(factorData.coefficient, 1) != meanDimension
        error(
            "Coefficient matrix row dimension is not valid for factor $(factorData.label). " *
            "Expected $meanDimension rows from factor mean vector length, " *
            "but got $(size(factorData.coefficient, 1))."
        )
    end

    totalStateDimension = 0

    for variableRef in factorData.variables
        totalStateDimension += variableDimension(variables, referenceIndex, variableRef)
    end

    if size(factorData.coefficient, 2) != totalStateDimension
        error(
            "Coefficient matrix column dimension is not valid for factor $(factorData.label). " *
            "Expected $totalStateDimension columns from connected variable dimensions, " *
            "but got $(size(factorData.coefficient, 2))."
        )
    end

    if size(factorData.covariance, 1) != meanDimension
        error(
            "Covariance row dimension is not valid for factor $(factorData.label). " *
            "Expected $meanDimension rows, " *
            "but got $(size(factorData.covariance, 1))."
        )
    end

    if size(factorData.covariance, 2) != meanDimension
        error(
            "Covariance column dimension is not valid for factor $(factorData.label). " *
            "Expected $meanDimension columns, " *
            "but got $(size(factorData.covariance, 2))."
        )
    end

    if !isposdef(Symmetric(factorData.covariance))
        error("Covariance matrix must be positive definite for factor $(factorData.label).")
    end

    if factorData.initialize && length(factorData.variables) != 1
        error("Only unary factors can be used for variable initialization.")
    end

    return nothing
end

function validateInitializingUnaryFactors(
    factors::Vector{GaussianFactor},
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
    factors::Vector{GaussianFactor},
    variables::Vector{GaussianVariable},
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
    factors::Vector{GaussianFactor},
    variables::Vector{GaussianVariable},
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
    graph::GaussianFactorGraph,
    factorData::GaussianFactor
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

"""
    coefficientBlock(graph::GaussianFactorGraph; factor::FactorRef, variable::VariableRef)

Return the coefficient submatrix in `GaussianFactor` corresponding to one variable.

# Arguments

- `graph`: Gaussian factor graph.

# Keywords

- `factor`: Factor index or label.
- `variable`: Variable id or label.

# Returns

The coefficient block for the selected variable.

# Notes

This is useful when inspecting a multi-GaussianVariable GaussianFactor. The returned block uses
the same row dimension as `GaussianFactor.coefficient` and the column dimension of the
selected variable.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
x2 = GaussianVariable(:x2, 1)
f1 = GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.2; label = "f1")

graph = factorGraph([x1, x2], [f1])

coefficientBlock(graph; factor = "f1", variable = :x2)
```
"""
function coefficientBlock(
    variables::Vector{GaussianVariable},
    referenceIndex::Dict{VariableRef, Int},
    factorData::GaussianFactor,
    variableRef::VariableRef
)
    requestedIndex = variableIndex(referenceIndex, variableRef)
    startColumn = 1

    for reference in factorData.variables
        currentIndex = variableIndex(referenceIndex, reference)
        dimension = variables[currentIndex].dimension
        stopColumn = startColumn + dimension - 1

        if currentIndex == requestedIndex
            return factorData.coefficient[:, startColumn:stopColumn]
        end

        startColumn = stopColumn + 1
    end

    error("Variable $variableRef is not connected to factor $(factorData.label).")
end

"""
    coefficientBlocks(graph::GaussianFactorGraph, factor::FactorRef)

Return coefficient submatrices for all variables connected to a GaussianFactor.

# Arguments

- `graph`: Gaussian factor graph.
- `factor`: Factor index or label.

# Returns

The blocks in the same order as `GaussianFactor.variables`.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
x2 = GaussianVariable(:x2, 1)
f1 = GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.2; label = "f1")

graph = factorGraph([x1, x2], [f1])

coefficientBlocks(graph, "f1")
```
"""
function coefficientBlocks(
    variables::Vector{GaussianVariable},
    referenceIndex::Dict{VariableRef, Int},
    factorData::GaussianFactor
)
    blocks = Matrix{Float64}[]

    startColumn = 1

    for variableRef in factorData.variables
        index = variableIndex(referenceIndex, variableRef)
        dimension = variables[index].dimension
        stopColumn = startColumn + dimension - 1

        push!(blocks, factorData.coefficient[:, startColumn:stopColumn])

        startColumn = stopColumn + 1
    end

    return blocks
end

"""
    treeFactorGraph(graph::AbstractFactorGraph; root = nothing)

Create a tree-oriented factor graph view rooted at `root`.

# Arguments

- `graph`: Gaussian or discrete factor graph.

# Keywords

- `root`: Optional root variable id or label. If omitted, the first variable is used.

# Returns

A [`TreeFactorGraph`](@ref) view that shares the original graph.

# Notes

The graph must be connected and must have exactly `nodes - 1` edges.
The returned tree view does not copy variables, factors, or messages.

# Example

```julia
x1 = DiscreteVariable(:x1, 2; states = [:off, :on])
x2 = DiscreteVariable(:x2, 2; states = [:low, :high])

f1 = DiscreteFactor(:x1, [0.8, 0.2]; label = "f1")
f2 = DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "f2")

graph = factorGraph([x1, x2], [f1, f2])
tree = treeFactorGraph(graph; root = :x1)
```
"""
function treeFactorGraph(
    graph::GaussianFactorGraph;
    root::Union{Nothing, VariableRef} = nothing
)
    return treeFactorGraphImpl(graph; root = root)
end

"""
    printModel(graph::AbstractFactorGraph)

Print a compact summary of variables and factors.

# Arguments

- `graph`: Gaussian factor graph, discrete factor graph, or tree view.

# Notes

This helper is intended for interactive inspection.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])

printModel(graph)
```
"""
function printModel(graph::GaussianFactorGraph)
    println("Gaussian factor graph")
    println("Number of GaussianVariable nodes: ", length(graph.variables))
    println("Number of GaussianFactor nodes: ", length(graph.factors))
    println()

    println("GaussianVariable nodes:")
    for variable in graph.variables
        println(
            "  label = ",
            variable.label,
            ", id = ",
            variable.id,
            ", dimension = ",
            variable.dimension,
            ", components = ",
            formattedRefs(variable.components)
        )
    end

    println()
    println("GaussianFactor nodes:")
    for factorData in graph.factors
        println("  label = ", factorData.label, ", id = ", factorData.id)
        println("    variables = ", formattedRefs(factorData.variables))
        println("    mean dimension = ", length(factorData.mean))
        println("    coefficient size = ", size(factorData.coefficient))
        println("    covariance size = ", size(factorData.covariance))
        println("    initialize = ", factorData.initialize)
    end

    return nothing
end

function printModel(tree::TreeFactorGraph{GaussianFactorGraph})
    return printModel(tree.graph)
end
