"""
    sumproduct(graph::DiscreteFactorGraph)

Create a discrete sum-product inference state for `graph`.

The returned [`DiscreteSumProductInference`](@ref) stores probability-vector
messages and marginals. Discrete variables with an initializing unary discrete
factor use that factor as their initial belief; other variables use
their own probability vector or a uniform default.
"""
function sumproduct(graph::DiscreteFactorGraph)
    initial = Vector{Float64}[]

    for variableIndexValue in eachindex(graph.variables)
        push!(initial, initialMessage(graph, variableIndexValue))
    end

    variableToFactor = Vector{Float64}[]
    factorToVariable = Vector{Float64}[]

    for edge in graph.edges
        message = initial[edge.variableIndex]

        push!(variableToFactor, copy(message))
        push!(factorToVariable, copy(message))
    end

    return DiscreteSumProductInference(
        initial,
        copyDiscreteMessages(variableToFactor),
        factorToVariable,
        copyDiscreteMessages(initial),
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

function assertDiscreteSumProductInferenceMatchesGraph(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference
)
    assertVariableStorageMatchesGraph(
        graph,
        "DiscreteSumProductInference does not match the DiscreteFactorGraph.",
        inference.initial,
        inference.marginal
    )

    assertCommonInferenceStorageMatchesGraph(
        graph,
        inference,
        "DiscreteSumProductInference does not match the DiscreteFactorGraph."
    )
    assertInferenceCurrent(graph, inference)

    return nothing
end

function assertInferenceMatchesGraph(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference
)
    return assertDiscreteSumProductInferenceMatchesGraph(graph, inference)
end

function extendDiscreteSumProductInferenceForAddedVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    variableIndexValue::Int
)
    message = initialMessage(graph, variableIndexValue)

    push!(inference.initial, message)
    push!(inference.marginal, copy(message))
    push!(inference.frozenVariables, false)

    return inference
end

function extendDiscreteSumProductInferenceForAddedFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
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

function setSumProductInitialFromUnaryFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    factorData::DiscreteFactor
)
    if !factorData.initialize
        return nothing
    end

    variableIdx = variableIndex(graph.referenceIndex, only(factorData.variables))
    message = initialMessageFromUnaryFactor(factorData)

    copyto!(inference.initial[variableIdx], message)
    copyto!(inference.marginal[variableIdx], message)

    for edgeId in graph.variableEdges[variableIdx]
        copyto!(inference.variableToFactor[edgeId], message)
    end

    return nothing
end

"""
    addVariable!(
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference,
        variable::DiscreteVariable
    )

Add a discrete variable to a graph and extend a matching sum-product inference
state in place.

The inference object must be current with the graph topology. The added
variable receives its initial sum-product message and marginal.
"""
function addVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    variable::DiscreteVariable
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)

    added = addVariable!(graph, variable)
    extendDiscreteSumProductInferenceForAddedVariable!(
        graph,
        inference,
        length(graph.variables)
    )
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addVariable!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
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
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference,
        factorData::DiscreteFactor
    )

Add a discrete factor to a graph and extend a matching sum-product inference
state in place.

New edge messages are initialized from the connected variable beliefs. If the
factor is an initializing unary factor, it updates the connected variable's
initial belief and incident variable-to-factor messages.
"""
function addFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    factorData::DiscreteFactor
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)

    firstNewEdge = length(graph.edges) + 1
    added = addFactor!(graph, factorData)
    extendDiscreteSumProductInferenceForAddedFactor!(graph, inference, firstNewEdge)
    setSumProductInitialFromUnaryFactor!(graph, inference, added)
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
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
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference;
        damping = false, prob = 0.6, alpha = 0.4, rng = Random.default_rng()
    )

Update all factor-to-variable sum-product messages in place.

When `damping` is enabled, each updated message is mixed with its previous value
with probability `prob` and damping weight `alpha`. Per-edge damping settings
configured with [`dampEdges!`](@ref) are also honored.
"""
function factorToVariableMessages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    rng = Random.default_rng()
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)
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
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference
    )

Update all variable-to-factor sum-product messages in place.
"""
function variableToFactorMessages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)
    ensureMessageDimensions(graph, inference)

    updateVariableToFactorMessages!(graph, inference, inference.variableToFactor)

    return nothing
end

"""
    messages!(
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference;
        damping = false, prob = 0.6, alpha = 0.4, schedule = nothing,
        updateFraction = nothing, updateCount = nothing, rng = Random.default_rng()
    )

Run one discrete sum-product message update sweep.

The default schedule updates factor-to-variable messages and then
variable-to-factor messages. Pass `schedule = :flooding` for simultaneous
message commits, or pass `schedule = :residual` with `updateFraction` or
`updateCount` for one residual-scheduled step.
"""
function messages!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    updateFraction = nothing,
    updateCount = nothing,
    rng = Random.default_rng()
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)
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
    marginals!(graph::DiscreteFactorGraph, inference::DiscreteSumProductInference)

Recompute all discrete sum-product marginals from the current incoming
factor-to-variable messages.
"""
function marginals!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)

    for variableIndexValue in eachindex(graph.variables)
        marginal = inference.marginal[variableIndexValue]
        hasIncomingMessage = false
        fill!(marginal, 1.0)

        for edgeId in graph.variableEdges[variableIndexValue]
            hasIncomingMessage = true
            marginal .*= inference.factorToVariable[edgeId]
        end

        if !hasIncomingMessage
            copyto!(marginal, inference.initial[variableIndexValue])
        end

        normalizeDiscreteVector!(marginal, "DiscreteVariable marginal")
    end

    return nothing
end

function stepGBP!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
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
    marginals!(graph, inference)

    return nothing
end

"""
    gbp!(
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference;
        iterations = 10, tolerance = nothing,
        damping = false, prob = 0.6, alpha = 0.4, schedule = nothing,
        updateFraction = nothing, updateCount = nothing, rng = Random.default_rng()
    )

Run iterative discrete sum-product belief propagation in place.

Each iteration updates messages and recomputes marginals. If `tolerance` is
provided, iteration stops once the maximum marginal change is at most that
tolerance.
"""
function gbp!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
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
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)

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
        previous =
            tolerance === nothing ?
            nothing :
            copyDiscreteMessages(inference.marginal)

        stepGBP!(
            graph,
            inference;
            damping = damping,
            prob = prob,
            alpha = alpha,
            schedule = selectedSchedule,
            rng = rng
        )

        if tolerance !== nothing && maxMarginalChange(graph, inference, previous) <= tolerance
            break
        end
    end

    return nothing
end

"""
    marginalProbability(
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference,
        variableRef::VariableRef, state::StateRef
    )

Return the stored marginal probability for one variable state.
"""
function marginalProbability(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    variableRef::VariableRef,
    state::StateRef
)
    variableIdx = variableIndex(graph, variableRef)
    stateIdx = stateIndex(graph, variableRef, state)

    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)

    return inference.marginal[variableIdx][stateIdx]
end


"""
    updateFactor!(
        graph::DiscreteFactorGraph, inference::DiscreteSumProductInference,
        factorRef::FactorRef; table = nothing, initialize = nothing
    )

Update a discrete factor and refresh a matching sum-product inference state.

If an initializing unary factor is updated, the connected variable's initial
belief and incident variable-to-factor messages are refreshed.
"""
function updateFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference,
    factorRef::FactorRef;
    table = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    assertDiscreteSumProductInferenceMatchesGraph(graph, inference)
    validateDiscreteInferenceFactorUpdate(
        graph,
        factorRef;
        table = table,
        initialize = initialize
    )

    updated = updateFactor!(
        graph,
        factorRef;
        table = table,
        initialize = initialize
    )
    setSumProductInitialFromUnaryFactor!(graph, inference, updated)

    return updated
end

function updateFactor!(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
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
