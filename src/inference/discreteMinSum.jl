"""
    minsum(graph::DiscreteFactorGraph)

Create a discrete min-sum inference state for `graph`.

The returned [`DiscreteMinSumInference`](@ref) stores cost-vector messages and
MAP state estimates. Variable probability vectors and initializing unary
discrete factors are converted to negative-log costs.
"""
function minsum(graph::DiscreteFactorGraph)
    initial = Vector{Float64}[]

    for variableIndexValue in eachindex(graph.variables)
        push!(initial, initialCost(graph, variableIndexValue))
    end

    variableToFactor = Vector{Float64}[]
    factorToVariable = Vector{Float64}[]

    for edge in graph.edges
        message = initial[edge.variableIndex]

        push!(variableToFactor, copy(message))
        push!(factorToVariable, copy(message))
    end

    estimate = StateRef[]

    for variableIndexValue in eachindex(graph.variables)
        push!(estimate, currentEstimate(graph, initial[variableIndexValue], variableIndexValue))
    end

    return DiscreteMinSumInference(
        initial,
        copyDiscreteMessages(variableToFactor),
        factorToVariable,
        estimate,
        Vector{Float64}[],
        Vector{Float64}[],
        falses(length(graph.edges)),
        falses(length(graph.variables)),
        falses(length(graph.factors)),
        falses(length(graph.edges)),
        fill(0.6, length(graph.edges)),
        fill(0.4, length(graph.edges)),
        graph.topologyVersion
    )
end

function assertDiscreteMinSumInferenceMatchesGraph(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference
)
    assertVariableStorageMatchesGraph(
        graph,
        "DiscreteMinSumInference does not match the DiscreteFactorGraph.",
        inference.initial,
        inference.estimate
    )

    assertCommonInferenceStorageMatchesGraph(
        graph,
        inference,
        "DiscreteMinSumInference does not match the DiscreteFactorGraph."
    )
    assertInferenceCurrent(graph, inference)

    return nothing
end

function assertInferenceMatchesGraph(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference
)
    return assertDiscreteMinSumInferenceMatchesGraph(graph, inference)
end

function extendDiscreteMinSumInferenceForAddedVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    variableIndexValue::Int
)
    message = initialCost(graph, variableIndexValue)

    push!(inference.initial, message)
    push!(inference.estimate, currentEstimate(graph, message, variableIndexValue))
    push!(inference.frozenVariables, false)

    return inference
end

function extendDiscreteMinSumInferenceForAddedFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    firstNewEdge::Int
)
    push!(inference.frozenFactors, false)

    for edgeId in firstNewEdge:length(graph.edges)
        edge = graph.edges[edgeId]
        message = inference.initial[edge.variableIndex]

        push!(inference.variableToFactor, copy(message))
        push!(inference.factorToVariable, copy(message))
        push!(inference.frozenEdges, false)
        push!(inference.dampedEdges, false)
        push!(inference.edgeDampingProb, 0.6)
        push!(inference.edgeDampingAlpha, 0.4)
    end

    return inference
end

function setDiscreteMinSumInitialFromUnaryFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    factorData::DiscreteFactor
)
    if !factorData.initialize
        return nothing
    end

    variableIdx = variableIndex(graph.referenceIndex, only(factorData.variables))
    message = probabilityToCost(initialMessageFromUnaryFactor(factorData))

    copyto!(inference.initial[variableIdx], message)
    inference.estimate[variableIdx] = currentEstimate(graph, message, variableIdx)

    for edgeId in graph.variableEdges[variableIdx]
        copyto!(inference.variableToFactor[edgeId], message)
    end

    return nothing
end

"""
    addVariable!(
        graph::DiscreteFactorGraph, inference::DiscreteMinSumInference,
        variable::DiscreteVariable
    )

Add a discrete variable to a graph and extend a matching min-sum inference
state in place.

The inference object must be current with the graph topology. The added
variable receives its initial min-sum cost message and MAP estimate.
"""
function addVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    variable::DiscreteVariable
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)

    added = addVariable!(graph, variable)
    extendDiscreteMinSumInferenceForAddedVariable!(
        graph,
        inference,
        length(graph.variables)
    )
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    id::VariableId,
    cardinality::Int;
    label::String = string(id),
    states = defaultStateRefs(cardinality),
    probability = nothing
)
    return addVariable!(
        graph,
        inference,
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
        graph::DiscreteFactorGraph, inference::DiscreteMinSumInference,
        factorData::DiscreteFactor
    )

Add a discrete factor to a graph and extend a matching min-sum inference state
in place.

New edge messages are initialized from the connected variable costs. If the
factor is an initializing unary factor, it updates the connected variable's
initial cost and MAP estimate.
"""
function addFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    factorData::DiscreteFactor
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)

    firstNewEdge = length(graph.edges) + 1
    added = addFactor!(graph, factorData)
    extendDiscreteMinSumInferenceForAddedFactor!(graph, inference, firstNewEdge)
    setDiscreteMinSumInitialFromUnaryFactor!(graph, inference, added)
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    args...;
    label::String = "",
    initialize::Bool = false
)
    return addFactor!(
        graph,
        inference,
        DiscreteFactor(args...; label = label, initialize = initialize)
    )
end

"""
    factorToVariableMessages!(
        graph::DiscreteFactorGraph, inference::DiscreteMinSumInference;
        damping = false, prob = 0.6, alpha = 0.4, rng = Random.default_rng()
    )

Update all factor-to-variable min-sum cost messages in place.

When `damping` is enabled, each updated message is mixed with its previous value
with probability `prob` and damping weight `alpha`. Per-edge damping settings
configured with [`dampEdges!`](@ref) are also honored.
"""
function factorToVariableMessages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    rng = Random.default_rng()
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)
    ensureMessageDimensions(graph, inference)
    validateDiscreteDamping(prob, alpha)

    hasDamping = damping || any(inference.dampedEdges)
    previousFactorToVariable =
        hasDamping ? copyDiscreteMessages(inference.factorToVariable) : Vector{Float64}[]

    updateFactorToVariableMessages!(graph, inference, inference.factorToVariable)

    if hasDamping
        applyFactorToVariableDamping!(
            graph,
            inference,
            inference.factorToVariable,
            previousFactorToVariable;
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
        graph::DiscreteFactorGraph, inference::DiscreteMinSumInference
    )

Update all variable-to-factor min-sum cost messages in place.
"""
function variableToFactorMessages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)
    ensureMessageDimensions(graph, inference)

    updateVariableToFactorMessages!(graph, inference, inference.variableToFactor)

    return nothing
end

"""
    messages!(
        graph::DiscreteFactorGraph, inference::DiscreteMinSumInference;
        damping = false, prob = 0.6, alpha = 0.4, schedule = nothing,
        updateFraction = nothing, updateCount = nothing, rng = Random.default_rng()
    )

Run one discrete min-sum message update sweep.

The default schedule updates factor-to-variable messages and then
variable-to-factor messages. Pass `schedule = :flooding` for simultaneous
message commits, or pass `schedule = :residual` with `updateFraction` or
`updateCount` for one residual-scheduled step.
"""
function messages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    updateFraction = nothing,
    updateCount = nothing,
    rng = Random.default_rng()
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)
    ensureMessageDimensions(graph, inference)
    validateDiscreteDamping(prob, alpha)
    selectedSchedule = effectiveSchedule(inference, schedule)
    validateGBPSchedule(selectedSchedule)

    if selectedSchedule == :residual
        residual = residualSchedule(
            graph,
            inference;
            updateFraction = updateFraction,
            updateCount = updateCount
        )
        return messages!(
            graph,
            inference,
            residual;
            damping = damping,
            prob = prob,
            alpha = alpha,
            rng = rng
        )
    end

    if selectedSchedule == :flooding
        ensureFloodingBuffers!(inference)
        updateFactorToVariableMessages!(graph, inference, inference.nextFactorToVariable)

        if damping || any(inference.dampedEdges)
            applyFactorToVariableDamping!(
                graph,
                inference,
                inference.nextFactorToVariable,
                inference.factorToVariable;
                damping = damping,
                prob = prob,
                alpha = alpha,
                rng = rng
            )
        end

        updateVariableToFactorMessages!(graph, inference, inference.nextVariableToFactor)
        commitFactorToVariableMessages!(inference)
        commitVariableToFactorMessages!(inference)

        return nothing
    end

    factorToVariableMessages!(
        graph,
        inference;
        damping = damping,
        prob = prob,
        alpha = alpha,
        rng = rng
    )
    variableToFactorMessages!(graph, inference)

    return nothing
end

"""
    estimates!(graph::DiscreteFactorGraph, inference::DiscreteMinSumInference)

Recompute all discrete min-sum MAP estimates from the current incoming
factor-to-variable cost messages.
"""
function estimates!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)

    for variableIndexValue in eachindex(graph.variables)
        costs = copy(inference.initial[variableIndexValue])
        hasIncomingMessage = false
        fill!(costs, 0.0)

        for edgeId in graph.variableEdges[variableIndexValue]
            hasIncomingMessage = true
            costs .+= inference.factorToVariable[edgeId]
        end

        if !hasIncomingMessage
            copyto!(costs, inference.initial[variableIndexValue])
        end

        inference.estimate[variableIndexValue] =
            currentEstimate(graph, costs, variableIndexValue)
    end

    return nothing
end

function stepGBP!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    rng = Random.default_rng()
)
    messages!(
        graph,
        inference;
        damping = damping,
        prob = prob,
        alpha = alpha,
        schedule = schedule,
        rng = rng
    )
    estimates!(graph, inference)

    return nothing
end

"""
    gbp!(
        graph::DiscreteFactorGraph, inference::DiscreteMinSumInference;
        iterations = 10, tolerance = nothing,
        damping = false, prob = 0.6, alpha = 0.4, schedule = nothing,
        updateFraction = nothing, updateCount = nothing, rng = Random.default_rng()
    )

Run iterative discrete min-sum belief propagation in place.

Each iteration updates messages and recomputes MAP estimates. If `tolerance` is
provided, iteration stops once the maximum estimate change is at most that
tolerance.
"""
function gbp!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    iterations::Int = 10,
    tolerance = nothing,
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    updateFraction::Union{Nothing, Float64} = nothing,
    updateCount::Union{Nothing, Int} = nothing,
    rng = Random.default_rng()
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)

    if iterations < 0
        error("Number of iterations must be nonnegative.")
    end

    validateDiscreteDamping(prob, alpha)
    validateTolerance(tolerance)
    selectedSchedule = effectiveSchedule(inference, schedule)
    validateGBPSchedule(selectedSchedule)

    if selectedSchedule == :residual
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
        previous = tolerance === nothing ? nothing : copy(inference.estimate)

        stepGBP!(
            graph,
            inference;
            damping = damping,
            prob = prob,
            alpha = alpha,
            schedule = selectedSchedule,
            rng = rng
        )

        if tolerance !== nothing && maxEstimateChange(graph, inference, previous) <= tolerance
            break
        end
    end

    return nothing
end

"""
    maxEstimateChange(
        graph::AbstractFactorGraph, inference::AbstractMinSumInference,
        previous::Union{Vector, AbstractMinSumInference}
    )

Return the largest change between current and previous MAP estimates.

# Arguments

- `graph`: Matching factor graph or tree view.
- `inference`: Current Gaussian or discrete min-sum inference object.
- `previous`: Previous estimate snapshot or previous matching min-sum inference object.

# Returns

For Gaussian min-sum, the largest Euclidean change between corresponding
numeric MAP estimate vectors. For discrete min-sum, `0.0` if all MAP state
references match, otherwise `1.0`.

# Example

```julia
oldEstimates = deepcopy(inference.estimate)
gbp!(graph, inference; iterations = 1)

maxEstimateChange(graph, inference, oldEstimates)
```
"""
function maxEstimateChange(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    previousEstimate::Vector{StateRef}
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)

    if length(inference.estimate) != length(previousEstimate)
        error("Discrete min-sum estimate vectors do not have the same length.")
    end

    maximumChange = 0.0

    for variableIndexValue in eachindex(inference.estimate)
        maximumChange = max(
            maximumChange,
            inference.estimate[variableIndexValue] == previousEstimate[variableIndexValue] ?
            0.0 :
            1.0
        )
    end

    return maximumChange
end

function maxEstimateChange(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    previous::DiscreteMinSumInference
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, previous)

    return maxEstimateChange(graph, inference, previous.estimate)
end


"""
    updateFactor!(
        graph::DiscreteFactorGraph, inference::DiscreteMinSumInference,
        factorRef::FactorRef; table = nothing, initialize = nothing
    )

Update a discrete factor and refresh a matching min-sum inference state.

If an initializing unary factor is updated, the connected variable's initial
cost and MAP estimate are refreshed.
"""
function updateFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference,
    factorRef::FactorRef;
    table = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)

    updated = updateFactor!(
        graph,
        factorRef;
        table = table,
        initialize = initialize
    )
    setDiscreteMinSumInitialFromUnaryFactor!(graph, inference, updated)

    return updated
end

function updateFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    factor::FactorRef,
    table = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    return updateFactor!(
        graph,
        inference,
        factor;
        table = table,
        initialize = initialize
    )
end
