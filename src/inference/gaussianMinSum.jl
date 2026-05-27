function zeroQuadraticMessage(dimension::Int)
    return QuadraticMessage(
        zeros(dimension, dimension),
        zeros(dimension),
        0.0
    )
end

function zeroQuadraticMessage(variable::GaussianVariable)
    return zeroQuadraticMessage(variable.dimension)
end

function quadraticMessageFromMeanCovariance(mean::Vector{Float64}, covariance::Matrix{Float64})
    precision = inverseFromPositiveDefinite(covariance)
    information = precision * mean

    return QuadraticMessage(precision, information, 0.0)
end

function initialQuadraticMessage(
    variable::GaussianVariable;
    mean = 0.0,
    covariance = 1e6
)
    messageMean =
        variable.mean === nothing ?
        defaultGaussianMean(variable, mean) :
        variable.mean

    messageCovariance =
        variable.covariance === nothing ?
        defaultGaussianCovariance(variable, covariance) :
        variable.covariance

    return quadraticMessageFromMeanCovariance(messageMean, messageCovariance)
end

function initialQuadraticMessageFromUnaryFactor(factorData::GaussianFactor)
    R = inverseFromPositiveDefinite(factorData.covariance)
    H = factorData.coefficient
    precision = symmetricPart(H' * R * H)
    information = H' * R * factorData.mean

    if !isposdef(Symmetric(precision))
        error(
            "Initializing unary GaussianFactor $(factorData.label) does not define " *
            "a positive definite precision."
        )
    end

    return QuadraticMessage(precision, information, 0.0)
end

function initialQuadraticMessage(
    graph::GaussianFactorGraph,
    variableIndexValue::Int;
    mean = 0.0,
    covariance = 1e6
)
    factorData = initializingUnaryFactor(graph, variableIndexValue)

    if factorData !== nothing
        return initialQuadraticMessageFromUnaryFactor(factorData)
    end

    return initialQuadraticMessage(
        graph.variables[variableIndexValue];
        mean = mean,
        covariance = covariance
    )
end

function copyQuadraticMessage(message::QuadraticMessage)
    return QuadraticMessage(copy(message.J), copy(message.h), message.c)
end

function copyQuadraticMessages(messages::Vector{QuadraticMessage})
    return copyMessages(messages, copyQuadraticMessage)
end

function copyMessage!(output::QuadraticMessage, source::QuadraticMessage)
    copyto!(output.J, source.J)
    copyto!(output.h, source.h)
    output.c = source.c

    return output
end

function dampMessage!(
    output::QuadraticMessage,
    previous::QuadraticMessage,
    current::QuadraticMessage,
    alpha::Float64
)
    output.J .= alpha .* previous.J .+ (1.0 - alpha) .* current.J
    output.h .= alpha .* previous.h .+ (1.0 - alpha) .* current.h
    output.c = alpha * previous.c + (1.0 - alpha) * current.c

    return output
end

"""
    minsum(graph::GaussianFactorGraph; mean = 0.0, covariance = 1e6)

Create a Gaussian min-sum MAP inference state.

# Arguments

- `graph`: Gaussian factor graph or tree view.

# Keywords

- `mean`: Default prior mean for variables without their own mean.
- `covariance`: Default prior covariance for variables without their own covariance.

# Returns

A [`GaussianMinSumInference`](@ref) object.

# Notes

Gaussian min-sum runs in negative-log quadratic cost form. Messages are
quadratic costs, and the stored estimates are MAP estimates. Marginal
covariances are not computed by this inference object. The default `mean` and
`covariance` are converted to initial quadratic messages for variables without
their own initial belief.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = minsum(graph; covariance = 1e4)
```
"""
function minsum(
    graph::GaussianFactorGraph;
    mean = 0.0,
    covariance = 1e6
)
    initial = QuadraticMessage[]

    for variableIndexValue in eachindex(graph.variables)
        push!(
            initial,
            initialQuadraticMessage(
                graph,
                variableIndexValue;
                mean = mean,
                covariance = covariance
            )
        )
    end

    variableToFactor = QuadraticMessage[]
    factorToVariable = QuadraticMessage[]
    nextVariableToFactor = QuadraticMessage[]
    nextFactorToVariable = QuadraticMessage[]

    for edge in graph.edges
        message = initial[edge.variableIndex]

        push!(variableToFactor, copyQuadraticMessage(message))
        push!(factorToVariable, copyQuadraticMessage(message))

    end

    estimate = [
        solveQuadraticEstimate(message.J, message.h)
        for message in initial
    ]

    return GaussianMinSumInference(
        variableToFactor,
        factorToVariable,
        estimate,
        nextVariableToFactor,
        nextFactorToVariable,
        falses(length(graph.edges)),
        falses(length(graph.variables)),
        falses(length(graph.factors)),
        falses(length(graph.edges)),
        fill(0.6, length(graph.edges)),
        fill(0.4, length(graph.edges)),
        mean,
        covariance,
        graph.topologyVersion
    )
end

function setMinSumInitialFromUnaryFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    factorData::GaussianFactor
)
    if !factorData.initialize
        return nothing
    end

    variableIdx = variableIndex(graph.referenceIndex, only(factorData.variables))
    message = initialQuadraticMessageFromUnaryFactor(factorData)

    inference.estimate[variableIdx] .= solveQuadraticEstimate(message.J, message.h)

    for edgeId in graph.variableEdges[variableIdx]
        copyMessage!(inference.variableToFactor[edgeId], message)
        copyMessage!(inference.factorToVariable[edgeId], message)
    end

    return nothing
end

function assertGaussianMinSumInferenceMatchesGraph(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference
)
    assertInferenceCurrent(graph, inference)

    assertVariableStorageMatchesGraph(
        graph,
        "Min-sum inference object does not match the graph variables.",
        inference.estimate,
        inference.frozenVariables
    )

    assertCommonInferenceStorageMatchesGraph(
        graph,
        inference,
        "Min-sum inference object does not match the graph."
    )

    return nothing
end

function assertInferenceMatchesGraph(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference
)
    return assertGaussianMinSumInferenceMatchesGraph(graph, inference)
end

function extendGaussianMinSumInferenceForAddedFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    firstNewEdge::Int
)
    push!(inference.frozenFactors, false)

    for edgeId in firstNewEdge:length(graph.edges)
        edge = graph.edges[edgeId]
        message = initialQuadraticMessage(
            graph,
            edge.variableIndex;
            mean = inference.defaultMean,
            covariance = inference.defaultCovariance
        )

        push!(inference.variableToFactor, copyQuadraticMessage(message))
        push!(inference.factorToVariable, copyQuadraticMessage(message))

        push!(inference.frozenEdges, false)
        push!(inference.dampedEdges, false)
        push!(inference.edgeDampingProb, 0.6)
        push!(inference.edgeDampingAlpha, 0.4)
    end

    return inference
end

function extendGaussianMinSumInferenceForAddedVariable!(
    inference::GaussianMinSumInference,
    variable::GaussianVariable;
    mean = inference.defaultMean,
    covariance = inference.defaultCovariance
)
    message = initialQuadraticMessage(
        variable;
        mean = mean,
        covariance = covariance
    )

    push!(inference.estimate, solveQuadraticEstimate(message.J, message.h))
    push!(inference.frozenVariables, false)

    return inference
end

function extendInferenceForAddedVariable!(
    inference::GaussianMinSumInference,
    variable::GaussianVariable;
    mean = inference.defaultMean,
    covariance = inference.defaultCovariance
)
    return extendGaussianMinSumInferenceForAddedVariable!(
        inference,
        variable;
        mean = mean,
        covariance = covariance
    )
end

"""
    addFactor!(
        graph::GaussianFactorGraph, inference::GaussianInference, factor::GaussianFactor
    )
    addFactor!(
        graph::GaussianFactorGraph, inference::GaussianInference,
        variable::Union{VariableRef, GaussianVariable}..., mean, coefficient, covariance;
        label = "", initialize = false
    )

Add a GaussianFactor node and extend an existing Gaussian inference state.

Existing messages and results are kept as a warm start. Moment inference extends
Gaussian messages and marginals, canonical inference extends canonical messages
and marginals, and min-sum inference extends quadratic messages and MAP
estimates. In all cases, the returned value is the added [`GaussianFactor`](@ref).

If `initialize = true`, the factor must be unary. The connected variable's
initial belief and incident messages are reset in the representation used by
the selected inference object.

# Example

```julia
x1 = GaussianVariable(:x1, 1)

graph = factorGraph([x1], GaussianFactor[])
inference = moment(graph)

addFactor!(graph, inference, :x1, 0.0, 1.0, 0.1; label = "f1")
```
"""
function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    factorData::GaussianFactor
)
    error(
        "Warm-start addFactor! supports GaussianMomentInference, " *
        "GaussianCanonicalInference, and GaussianMinSumInference."
    )
end

function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    args...;
    label::String = "",
    initialize::Bool = false
)
    return addFactor!(
        graph,
        inference,
        GaussianFactor(args...; label = label, initialize = initialize)
    )
end

function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    factorData::GaussianFactor
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)

    firstNewEdge = length(graph.edges) + 1
    added = addFactor!(graph, factorData)
    extendGaussianMinSumInferenceForAddedFactor!(graph, inference, firstNewEdge)
    setMinSumInitialFromUnaryFactor!(graph, inference, added)
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    args...;
    label::String = "",
    initialize::Bool = false
)
    return addFactor!(
        graph,
        inference,
        GaussianFactor(args...; label = label, initialize = initialize)
    )
end

"""
    updateFactor!(
        graph::GaussianFactorGraph, inference::GaussianInference, factor::FactorRef;
        mean = nothing, coefficient = nothing, covariance = nothing, initialize = nothing
    )
    updateFactor!(
        graph::GaussianFactorGraph, inference::GaussianInference;
        factor::FactorRef,
        mean = nothing, coefficient = nothing, covariance = nothing, initialize = nothing
    )

Update an existing GaussianFactor and refresh the matching Gaussian inference state.

The graph factor is updated first, then the inference object refreshes the
representation it owns: moment messages, canonical messages, or min-sum
quadratic messages and MAP estimates. When `initialize = true`, the updated
factor must be unary and its value is copied into the connected variable's
initial belief and incident messages.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

updateFactor!(graph, inference, "f1"; mean = 0.2)
```
"""
function updateFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    factorRef::FactorRef;
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    error(
        "Warm-start updateFactor! supports GaussianMomentInference, " *
        "GaussianCanonicalInference, and GaussianMinSumInference."
    )
end

function updateFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    factor::FactorRef,
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    return updateFactor!(
        graph,
        inference,
        factor;
        mean = mean,
        coefficient = coefficient,
        covariance = covariance,
        initialize = initialize
    )
end

function updateFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    factorRef::FactorRef;
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)

    updated = updateFactor!(
        graph,
        factorRef;
        mean = mean,
        coefficient = coefficient,
        covariance = covariance,
        initialize = initialize
    )

    if initialize === true
        setMinSumInitialFromUnaryFactor!(graph, inference, updated)
    end

    return updated
end

function updateFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    factor::FactorRef,
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    return updateFactor!(
        graph,
        inference,
        factor;
        mean = mean,
        coefficient = coefficient,
        covariance = covariance,
        initialize = initialize
    )
end

function unsupportedGaussianControlError()
    error(
        "Freezing and damping controls support GaussianMomentInference, " *
        "GaussianCanonicalInference, and GaussianMinSumInference."
    )
end

"""
    isFrozenFactor(
        graph::AbstractFactorGraph, inference::AbstractInference, factor::FactorRef
    )

Return whether outgoing messages from a factor are frozen.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.
- `factor`: Factor index or label.

# Returns

`true` if the selected factor is frozen.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

isFrozenFactor(graph, inference, "f1")
```
"""
function isFrozenFactor(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    factorRef::FactorRef
)
    unsupportedGaussianControlError()
end

"""
    isFrozenVariable(
        graph::AbstractFactorGraph, inference::AbstractInference, variable::VariableRef
    )

Return whether outgoing messages from a variable are frozen.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.
- `variable`: Variable id or label.

# Returns

`true` if the selected variable is frozen.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

isFrozenVariable(graph, inference, :x1)
```
"""
function isFrozenVariable(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    variableRef::VariableRef
)
    unsupportedGaussianControlError()
end

"""
    isFrozenEdge(
        graph::AbstractFactorGraph, inference::AbstractInference;
        variable::VariableRef, factor::FactorRef
    )

Return whether both message directions on an edge are frozen.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Keywords

- `variable`: Variable id or label.
- `factor`: Factor index or label.

# Returns

`true` if the selected edge is frozen.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

isFrozenEdge(graph, inference; variable = :x1, factor = "f1")
```
"""
function isFrozenEdge(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    variable::VariableRef,
    factor::FactorRef
)
    unsupportedGaussianControlError()
end

"""
    freezeFactor!(
        graph::AbstractFactorGraph, inference::AbstractInference, factor::FactorRef
    )

Freeze all outgoing factor-to-variable messages for a factor.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.
- `factor`: Factor index or label.

# Notes

While frozen, belief propagation iterations keep the selected outgoing messages
unchanged.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

freezeFactor!(graph, inference, "f1")
```
"""
function freezeFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    factorRef::FactorRef
)
    unsupportedGaussianControlError()
end

"""
    freezeVariable!(
        graph::AbstractFactorGraph, inference::AbstractInference, variable::VariableRef
    )

Freeze all outgoing variable-to-factor messages for a variable.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.
- `variable`: Variable id or label.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

freezeVariable!(graph, inference, :x1)
```
"""
function freezeVariable!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    variableRef::VariableRef
)
    unsupportedGaussianControlError()
end

"""
    freezeEdge!(
        graph::AbstractFactorGraph, inference::AbstractInference;
        variable::VariableRef, factor::FactorRef
    )

Freeze both message directions on the edge connecting `variable` and `factor`.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Keywords

- `variable`: Variable id or label.
- `factor`: Factor index or label.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

freezeEdge!(graph, inference; variable = :x1, factor = "f1")
```
"""
function freezeEdge!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    variable::VariableRef,
    factor::FactorRef
)
    unsupportedGaussianControlError()
end

"""
    unfreezeFactor!(
        graph::AbstractFactorGraph, inference::AbstractInference, factor::FactorRef
    )

Unfreeze outgoing factor-to-variable messages for a factor.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.
- `factor`: Factor index or label.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

unfreezeFactor!(graph, inference, "f1")
```
"""
function unfreezeFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    factorRef::FactorRef
)
    unsupportedGaussianControlError()
end

"""
    unfreezeVariable!(
        graph::AbstractFactorGraph, inference::AbstractInference, variable::VariableRef
    )

Unfreeze outgoing variable-to-factor messages for a variable.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.
- `variable`: Variable id or label.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

unfreezeVariable!(graph, inference, :x1)
```
"""
function unfreezeVariable!(
    graph::GaussianFactorGraph,
    inference::GaussianInference,
    variableRef::VariableRef
)
    unsupportedGaussianControlError()
end

"""
    unfreezeEdge!(
        graph::AbstractFactorGraph, inference::AbstractInference;
        variable::VariableRef, factor::FactorRef
    )

Unfreeze both message directions on the edge connecting `variable` and `factor`.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Keywords

- `variable`: Variable id or label.
- `factor`: Factor index or label.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

unfreezeEdge!(graph, inference; variable = :x1, factor = "f1")
```
"""
function unfreezeEdge!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    variable::VariableRef,
    factor::FactorRef
)
    unsupportedGaussianControlError()
end

"""
    isDampedEdge(
        graph::AbstractFactorGraph, inference::AbstractInference;
        variable = nothing, factor = nothing
    )

Return whether all selected edges are marked for damping.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Keywords

- `variable`: Optional variable id or label.
- `factor`: Optional factor index or label.

# Returns

`true` if all selected edges are damped.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

isDampedEdge(graph, inference; variable = :x1, factor = "f1")
```
"""
function isDampedEdge(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    unsupportedGaussianControlError()
end

"""
    dampEdges!(
        graph::AbstractFactorGraph, inference::AbstractInference;
        variable = nothing, factor = nothing, prob = 0.6, alpha = 0.4
    )

Enable damping for selected edges.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Keywords

- `variable`: Optional variable id or label.
- `factor`: Optional factor index or label.
- `prob`: Damping probability.
- `alpha`: Previous-message mixing weight.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
x2 = GaussianVariable(:x2, 1)
f1 = GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.1; label = "f1")

graph = factorGraph([x1, x2], [f1])
inference = moment(graph)

dampEdges!(graph, inference; variable = :x2, factor = "f1", prob = 1.0, alpha = 0.35)
```
"""
function dampEdges!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4
)
    unsupportedGaussianControlError()
end

"""
    undampEdges!(
        graph::AbstractFactorGraph, inference::AbstractInference;
        variable = nothing, factor = nothing
    )

Disable damping for selected edges.

# Arguments

- `graph`: Gaussian or discrete factor graph.
- `inference`: Matching inference object.

# Keywords

- `variable`: Optional variable id or label.
- `factor`: Optional factor index or label.

# Notes

The same selector rules as [`dampEdges!`](@ref) apply.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

undampEdges!(graph, inference; variable = :x1, factor = "f1")
```
"""
function undampEdges!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    variable::Union{Nothing, VariableRef} = nothing,
    factor::Union{Nothing, FactorRef} = nothing
)
    unsupportedGaussianControlError()
end



function localVariableRanges(graph::GaussianFactorGraph, factorData::GaussianFactor)
    ranges = UnitRange{Int}[]
    startIndex = 1

    for variableRef in factorData.variables
        variableIdx = variableIndex(graph.referenceIndex, variableRef)
        dimension = graph.variables[variableIdx].dimension
        stopIndex = startIndex + dimension - 1

        push!(ranges, startIndex:stopIndex)
        startIndex = stopIndex + 1
    end

    return ranges
end

function localIndexForEdge(
    graph::GaussianFactorGraph,
    factorData::GaussianFactor,
    edge::Edge
)
    for (localIndex, variableRef) in pairs(factorData.variables)
        variableIdx = variableIndex(graph.referenceIndex, variableRef)

        if variableIdx == edge.variableIndex
            return localIndex
        end
    end

    error("Edge $(edge.id) does not belong to GaussianFactor $(factorData.label).")
end

function factorQuadratic(graph::GaussianFactorGraph, factorData::GaussianFactor)
    R = inverseFromPositiveDefinite(factorData.covariance)
    J = symmetricPart(factorData.coefficient' * R * factorData.coefficient)
    h = factorData.coefficient' * R * factorData.mean

    return J, h
end

function variableToFactorMessage!(
    output::QuadraticMessage,
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    edgeId::Int
)
    edge = graph.edges[edgeId]

    fill!(output.J, 0.0)
    fill!(output.h, 0.0)
    output.c = 0.0

    if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
        return output
    end

    for incomingEdgeId in graph.variableEdges[edge.variableIndex]
        if incomingEdgeId == edgeId
            continue
        end

        incoming = inference.factorToVariable[incomingEdgeId]
        output.J .+= incoming.J
        output.h .+= incoming.h
        output.c += incoming.c
    end

    copySymmetricPart!(output.J, output.J)

    return output
end

function updateVariableToFactorMessagesStandard!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    output::Vector{QuadraticMessage}
)
    for edge in graph.edges
        variableToFactorMessage!(output[edge.id], graph, inference, edge.id)
    end

    return nothing
end

function updateVariableToFactorMessagesBroadcast!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    output::Vector{QuadraticMessage}
)
    for (variableIndexValue, variable) in pairs(graph.variables)
        edgeIds = graph.variableEdges[variableIndexValue]
        totalJ = zeros(variable.dimension, variable.dimension)
        totalH = zeros(variable.dimension)
        totalC = 0.0

        for edgeId in edgeIds
            incoming = inference.factorToVariable[edgeId]
            totalJ .+= incoming.J
            totalH .+= incoming.h
            totalC += incoming.c
        end

        for edgeId in edgeIds
            edge = graph.edges[edgeId]

            if inference.frozenEdges[edgeId] || inference.frozenVariables[edge.variableIndex]
                copyMessage!(output[edgeId], inference.variableToFactor[edgeId])
                continue
            end

            incoming = inference.factorToVariable[edgeId]
            output[edgeId].J .= totalJ .- incoming.J
            output[edgeId].h .= totalH .- incoming.h
            output[edgeId].c = totalC - incoming.c
            copySymmetricPart!(output[edgeId].J, output[edgeId].J)
        end
    end

    return nothing
end

function updateVariableToFactorMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    output::Vector{QuadraticMessage};
    broadcast::Bool = false
)
    if broadcast
        return updateVariableToFactorMessagesBroadcast!(graph, inference, output)
    end

    return updateVariableToFactorMessagesStandard!(graph, inference, output)
end

function factorToVariableMessage!(
    output::QuadraticMessage,
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    edgeId::Int
)
    edge = graph.edges[edgeId]

    if inference.frozenEdges[edgeId] || inference.frozenFactors[edge.factorIndex]
        return output
    end

    factorData = graph.factors[edge.factorIndex]
    J, h = factorQuadratic(graph, factorData)
    ranges = localVariableRanges(graph, factorData)
    targetLocalIndex = localIndexForEdge(graph, factorData, edge)
    targetRange = ranges[targetLocalIndex]

    for (localIndex, variableRef) in pairs(factorData.variables)
        if localIndex == targetLocalIndex
            continue
        end

        variableIdx = variableIndex(graph.referenceIndex, variableRef)
        localRange = ranges[localIndex]

        for incomingEdgeId in graph.variableEdges[variableIdx]
            incomingEdge = graph.edges[incomingEdgeId]

            if incomingEdge.factorIndex != edge.factorIndex
                continue
            end

            incoming = inference.variableToFactor[incomingEdgeId]
            J[localRange, localRange] .+= incoming.J
            h[localRange] .+= incoming.h

            break
        end
    end

    otherRanges = UnitRange{Int}[]

    for (localIndex, localRange) in pairs(ranges)
        if localIndex != targetLocalIndex
            push!(otherRanges, localRange)
        end
    end

    if isempty(otherRanges)
        copyto!(output.J, J[targetRange, targetRange])
        copyto!(output.h, h[targetRange])
        output.c = 0.0
        copySymmetricPart!(output.J, output.J)

        return output
    end

    otherIndices = reduce(vcat, [collect(localRange) for localRange in otherRanges])
    Jtt = J[targetRange, targetRange]
    Jtr = J[targetRange, otherIndices]
    Jrt = J[otherIndices, targetRange]
    Jrr = symmetricPart(J[otherIndices, otherIndices])
    ht = h[targetRange]
    hr = h[otherIndices]
    JrrInverse = pinv(Symmetric(Jrr))

    copyto!(output.J, symmetricPart(Jtt - Jtr * JrrInverse * Jrt))
    copyto!(output.h, ht - Jtr * JrrInverse * hr)
    output.c = 0.0

    return output
end

function setFactorToVariableMessageFromJoint!(
    output::QuadraticMessage,
    jointJ::Matrix{Float64},
    jointH::Vector{Float64},
    ranges::Vector{UnitRange{Int}},
    targetLocalIndex::Int
)
    targetRange = ranges[targetLocalIndex]
    otherRanges = UnitRange{Int}[]

    for (localIndex, localRange) in pairs(ranges)
        if localIndex != targetLocalIndex
            push!(otherRanges, localRange)
        end
    end

    if isempty(otherRanges)
        copyto!(output.J, jointJ[targetRange, targetRange])
        copyto!(output.h, jointH[targetRange])
        output.c = 0.0
        copySymmetricPart!(output.J, output.J)

        return output
    end

    otherIndices = reduce(vcat, [collect(localRange) for localRange in otherRanges])
    Jtt = jointJ[targetRange, targetRange]
    Jtr = jointJ[targetRange, otherIndices]
    Jrt = jointJ[otherIndices, targetRange]
    Jrr = symmetricPart(jointJ[otherIndices, otherIndices])
    ht = jointH[targetRange]
    hr = jointH[otherIndices]
    JrrInverse = pinv(Symmetric(Jrr))

    copyto!(output.J, symmetricPart(Jtt - Jtr * JrrInverse * Jrt))
    copyto!(output.h, ht - Jtr * JrrInverse * hr)
    output.c = 0.0

    return output
end

function updateFactorToVariableMessagesStandard!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    output::Vector{QuadraticMessage}
)
    for edge in graph.edges
        factorToVariableMessage!(output[edge.id], graph, inference, edge.id)
    end

    return nothing
end

function updateFactorToVariableMessagesBroadcast!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    output::Vector{QuadraticMessage}
)
    for (factorIndexValue, factorData) in pairs(graph.factors)
        J, h = factorQuadratic(graph, factorData)
        ranges = localVariableRanges(graph, factorData)
        incomingJs = Matrix{Float64}[]
        incomingHs = Vector{Float64}[]

        for (localIndex, variableRef) in pairs(factorData.variables)
            variableIdx = variableIndex(graph.referenceIndex, variableRef)
            localRange = ranges[localIndex]
            incomingJ = zeros(length(localRange), length(localRange))
            incomingH = zeros(length(localRange))

            for incomingEdgeId in graph.variableEdges[variableIdx]
                incomingEdge = graph.edges[incomingEdgeId]

                if incomingEdge.factorIndex != factorIndexValue
                    continue
                end

                incoming = inference.variableToFactor[incomingEdgeId]
                incomingJ .= incoming.J
                incomingH .= incoming.h

                break
            end

            J[localRange, localRange] .+= incomingJ
            h[localRange] .+= incomingH
            push!(incomingJs, incomingJ)
            push!(incomingHs, incomingH)
        end

        for edgeId in graph.factorEdges[factorIndexValue]
            edge = graph.edges[edgeId]

            if inference.frozenEdges[edgeId] || inference.frozenFactors[factorIndexValue]
                copyMessage!(output[edgeId], inference.factorToVariable[edgeId])
                continue
            end

            targetLocalIndex = localIndexForEdge(graph, factorData, edge)
            targetRange = ranges[targetLocalIndex]
            targetJ = copy(J)
            targetH = copy(h)
            targetJ[targetRange, targetRange] .-= incomingJs[targetLocalIndex]
            targetH[targetRange] .-= incomingHs[targetLocalIndex]
            setFactorToVariableMessageFromJoint!(
                output[edgeId],
                targetJ,
                targetH,
                ranges,
                targetLocalIndex
            )
        end
    end

    return nothing
end

function updateFactorToVariableMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    output::Vector{QuadraticMessage};
    broadcast::Bool = false
)
    if broadcast
        return updateFactorToVariableMessagesBroadcast!(graph, inference, output)
    end

    return updateFactorToVariableMessagesStandard!(graph, inference, output)
end

function applyFactorToVariableDamping!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    output::Vector{QuadraticMessage},
    previous::Vector{QuadraticMessage};
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
            dampMessage!(
                output[edge.id],
                previous[edge.id],
                output[edge.id],
                edgeAlpha
            )
        end
    end

    return output
end

"""
    factorToVariableMessages!(
        graph::GaussianFactorGraph, inference::GaussianInference;
        damping = false, prob = 0.6, alpha = 0.4, rng = Random.default_rng(),
        broadcast = false
    )

Run one GaussianFactor-to-GaussianVariable message update pass.

# Arguments

- `graph`: Gaussian factor graph.
- `inference`: Matching Gaussian inference object.

# Keywords

- `damping`: Enable global factor-to-variable message damping.
- `prob`: Damping probability.
- `alpha`: Previous-message mixing weight.
- `rng`: Random number generator used by damping.
- `broadcast`: Compute messages in grouped broadcast form.

# Notes

This updates only messages sent from GaussianFactor nodes to GaussianVariable
nodes. It does not recompute marginals or MAP estimates; call [`marginals!`](@ref)
for moment and canonical inference, or [`estimates!`](@ref) for min-sum inference.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

factorToVariableMessages!(graph, inference)
```
"""
function factorToVariableMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    unsupportedGaussianControlError()
end

function factorToVariableMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)
    validateDamping(prob, alpha)

    previous = copyQuadraticMessages(inference.factorToVariable)

    updateFactorToVariableMessages!(
        graph,
        inference,
        inference.factorToVariable;
        broadcast = broadcast
    )

    if damping || any(inference.dampedEdges)
        applyFactorToVariableDamping!(
            graph,
            inference,
            inference.factorToVariable,
            previous;
            damping = damping,
            prob = prob,
            alpha = alpha,
            rng = rng
        )
    end

    return nothing
end

"""
    variableToFactorMessages!(
        graph::GaussianFactorGraph, inference::GaussianInference;
        broadcast = false
    )

Run one GaussianVariable-to-GaussianFactor message update pass.

# Arguments

- `graph`: Gaussian factor graph.
- `inference`: Matching Gaussian inference object.

# Keywords

- `broadcast`: Compute messages in grouped broadcast form.

# Notes

This updates only messages sent from GaussianVariable nodes to GaussianFactor
nodes. It does not recompute marginals or MAP estimates.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

variableToFactorMessages!(graph, inference)
```
"""
function variableToFactorMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    broadcast::Bool = false
)
    unsupportedGaussianControlError()
end

function variableToFactorMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    broadcast::Bool = false
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)

    updateVariableToFactorMessages!(
        graph,
        inference,
        inference.variableToFactor;
        broadcast = broadcast
    )

    return nothing
end

"""
    messages!(
        graph::GaussianFactorGraph, inference::GaussianInference;
        damping = false, prob = 0.6, alpha = 0.4, rng = Random.default_rng(),
        schedule = nothing, broadcast = false
    )

Run one complete Gaussian message update step.

# Arguments

- `graph`: Gaussian factor graph.
- `inference`: Matching Gaussian inference object.

# Keywords

- `damping`: Enable global factor-to-variable message damping.
- `prob`: Damping probability.
- `alpha`: Previous-message mixing weight.
- `rng`: Random number generator used by damping.
- `schedule`: `:sequential`, `:flooding`, `:residual`, a schedule object, or `nothing`.
- `broadcast`: Compute sequential or flooding passes in grouped form.

# Notes

This updates messages stored in `inference`. Call [`marginals!`](@ref) afterwards
for moment and canonical inference, or [`estimates!`](@ref) for min-sum inference.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

messages!(graph, inference)
```
"""
function messages!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    unsupportedGaussianControlError()
end

function messages!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)
    selectedSchedule = effectiveSchedule(inference, schedule)
    validateGBPSchedule(selectedSchedule)

    if selectedSchedule == :residual
        error("Pass a ResidualSchedule object to messages! for residual scheduling.")
    end

    if selectedSchedule == :flooding
        validateDamping(prob, alpha)
        ensureFloodingBuffers!(inference)

        for edge in graph.edges
            copyMessage!(
                inference.nextFactorToVariable[edge.id],
                inference.factorToVariable[edge.id]
            )
            copyMessage!(
                inference.nextVariableToFactor[edge.id],
                inference.variableToFactor[edge.id]
            )
        end

        previous = copyQuadraticMessages(inference.factorToVariable)

        updateFactorToVariableMessages!(
            graph,
            inference,
            inference.nextFactorToVariable;
            broadcast = broadcast
        )

        if damping || any(inference.dampedEdges)
            applyFactorToVariableDamping!(
                graph,
                inference,
                inference.nextFactorToVariable,
                previous;
                damping = damping,
                prob = prob,
                alpha = alpha,
                rng = rng
            )
        end

        updateVariableToFactorMessages!(
            graph,
            inference,
            inference.nextVariableToFactor;
            broadcast = broadcast
        )

        for edge in graph.edges
            copyMessage!(
                inference.factorToVariable[edge.id],
                inference.nextFactorToVariable[edge.id]
            )
            copyMessage!(
                inference.variableToFactor[edge.id],
                inference.nextVariableToFactor[edge.id]
            )
        end

        return nothing
    end

    factorToVariableMessages!(
        graph,
        inference;
        damping = damping,
        prob = prob,
        alpha = alpha,
        broadcast = broadcast,
        rng = rng
    )
    variableToFactorMessages!(graph, inference; broadcast = broadcast)

    return nothing
end

function solveQuadraticEstimate(J::Matrix{Float64}, h::Vector{Float64})
    if isposdef(Symmetric(J))
        return Symmetric(J) \ h
    end

    return pinv(Symmetric(J)) * h
end

"""
    estimates!(graph::GaussianFactorGraph, inference::GaussianMinSumInference)

Recompute MAP estimates from current Gaussian min-sum messages.

# Arguments

- `graph`: Gaussian factor graph.
- `inference`: Matching Gaussian min-sum inference object.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = minsum(graph)

messages!(graph, inference)
estimates!(graph, inference)
```
"""
function estimates!(graph::GaussianFactorGraph, inference::GaussianMinSumInference)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)

    for (variableIdx, variable) in pairs(graph.variables)
        J = zeros(variable.dimension, variable.dimension)
        h = zeros(variable.dimension)

        for edgeId in graph.variableEdges[variableIdx]
            incoming = inference.factorToVariable[edgeId]
            J .+= incoming.J
            h .+= incoming.h
        end

        inference.estimate[variableIdx] .= solveQuadraticEstimate(symmetricPart(J), h)
    end

    return nothing
end

function stepGBP!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    messages!(
        graph,
        inference;
        damping = damping,
        prob = prob,
        alpha = alpha,
        schedule = schedule,
        broadcast = broadcast,
        rng = rng
    )
    estimates!(graph, inference)

    return nothing
end

function maxEstimateChange(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    previousEstimate::Vector{Vector{Float64}}
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)

    return maximumVectorChange(inference.estimate, previousEstimate)
end

function maxEstimateChange(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference,
    previous::GaussianMinSumInference
)
    assertGaussianMinSumInferenceMatchesGraph(graph, previous)

    return maxEstimateChange(graph, inference, previous.estimate)
end

"""
    gbp!(
        graph::GaussianFactorGraph, inference::GaussianInference;
        iterations = 10, tolerance = nothing,
        damping = false, prob = 0.6, alpha = 0.4, rng = Random.default_rng()
        schedule = nothing, updateFraction = nothing, updateCount = nothing
    )

Run Gaussian belief propagation updates in place.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching Gaussian inference object.

# Keywords

- `iterations`: Maximum number of belief propagation iterations to run.
- `tolerance`: Optional early-stopping tolerance.
- `damping`: Enable global factor-to-variable message damping.
- `prob`: Damping probability.
- `rng`: Random number generator used for damping.
- `alpha`: Damping mixing weight.
- `schedule`: `:sequential`, `:flooding`, `:residual`, or `nothing`.
- `broadcast`: Compute sequential or flooding message passes in grouped form.
- `updateFraction`: Residual-schedule fraction when `schedule = :residual`.
- `updateCount`: Residual-schedule count when `schedule = :residual`.

# Notes

Each iteration updates messages, then updates the result state owned by the
inference object: marginals for moment and canonical inference, or MAP estimates
for min-sum inference. When `schedule` is omitted, `gbp!` uses the sequential
schedule.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

gbp!(graph, inference; iterations = 20, tolerance = 1e-6)
```
"""
function gbp!(
    graph::GaussianFactorGraph,
    inference::GaussianInference;
    iterations::Int = 10,
    tolerance = nothing,
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    broadcast::Bool = false,
    updateFraction::Union{Nothing, Float64} = nothing,
    updateCount::Union{Nothing, Int} = nothing,
    rng = Random.default_rng()
)
    unsupportedGaussianControlError()
end

function gbp!(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    iterations::Int = 10,
    tolerance = nothing,
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    broadcast::Bool = false,
    updateFraction::Union{Nothing, Float64} = nothing,
    updateCount::Union{Nothing, Int} = nothing,
    rng = Random.default_rng()
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)

    if iterations < 0
        error("Number of iterations must be nonnegative.")
    end

    validateDamping(prob, alpha)
    validateTolerance(tolerance)
    selectedSchedule = effectiveSchedule(inference, schedule)
    validateGBPSchedule(selectedSchedule)

    if selectedSchedule == :residual
        if broadcast
            error("broadcast = true is not supported with schedule = :residual.")
        end

        runResidualGBP!(
            graph,
            inference;
            iterations = iterations,
            tolerance = tolerance,
            updateFraction = updateFraction,
            updateCount = updateCount,
            damping = damping,
            prob = prob,
            alpha = alpha,
            rng = rng
        )

        return nothing
    end

    for _ in 1:iterations
        previous =
            tolerance === nothing ?
            nothing :
            [copy(estimate) for estimate in inference.estimate]

        stepGBP!(
            graph,
            inference;
            damping = damping,
            prob = prob,
            alpha = alpha,
            schedule = selectedSchedule,
            broadcast = broadcast,
            rng = rng
        )

        if tolerance !== nothing &&
                maxEstimateChange(graph, inference, previous) <= tolerance
            break
        end
    end

    return nothing
end

function messages!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMinSumInference;
    kwargs...
)
    return messages!(tree.graph, inference; kwargs...)
end

function gbp!(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMinSumInference;
    kwargs...
)
    return gbp!(tree.graph, inference; kwargs...)
end
