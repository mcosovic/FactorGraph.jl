function inverseSPDCanonical(matrix::AbstractMatrix)
    if size(matrix, 1) != size(matrix, 2)
        error("Matrix must be square.")
    end

    symmetricMatrix = symmetricPart(matrix)

    if !isposdef(Symmetric(symmetricMatrix))
        error("Matrix must be positive definite.")
    end

    factorization = cholesky(Symmetric(symmetricMatrix))
    dimension = size(symmetricMatrix, 1)

    inverseMatrix = Matrix(factorization \ Matrix{Float64}(I, dimension, dimension))

    return symmetricPart(inverseMatrix)
end

function solveSPDCanonical(matrix::AbstractMatrix, rhs::AbstractVecOrMat)
    symmetricMatrix = symmetricPart(matrix)

    if !isposdef(Symmetric(symmetricMatrix))
        error("Matrix must be positive definite.")
    end

    return Symmetric(symmetricMatrix) \ rhs
end

mutable struct CanonicalMessage
    information::Vector{Float64}
    precision::Matrix{Float64}

    function CanonicalMessage(information, precision)
        informationVector = asVector(information)
        precisionMatrix = symmetricPart(asCovarianceMatrix(precision))

        if size(precisionMatrix, 1) != length(informationVector)
            error("Message precision row dimension must match information dimension.")
        end

        if size(precisionMatrix, 2) != length(informationVector)
            error("Message precision column dimension must match information dimension.")
        end

        return new(informationVector, precisionMatrix)
    end
end

mutable struct CanonicalMarginal
    information::Vector{Float64}
    precision::Matrix{Float64}
    mean::Vector{Float64}
    covariance::Matrix{Float64}
end

function setCanonicalMessage!(
    message::CanonicalMessage,
    information,
    precision
)
    informationVector = asVector(information)
    precisionMatrix = asCovarianceMatrix(precision)

    if length(informationVector) != length(message.information)
        error("Information dimension does not match existing message dimension.")
    end

    if size(precisionMatrix, 1) != length(message.information)
        error("Precision row dimension does not match existing message dimension.")
    end

    if size(precisionMatrix, 2) != length(message.information)
        error("Precision column dimension does not match existing message dimension.")
    end

    copyto!(message.information, informationVector)
    copySymmetricPart!(message.precision, precisionMatrix)

    return message
end

function copyCanonicalMessage(message::CanonicalMessage)
    return CanonicalMessage(
        copy(message.information),
        copy(message.precision)
    )
end

function copyCanonicalMessages(messages::Vector{CanonicalMessage})
    return copyMessages(messages, copyCanonicalMessage)
end

function copyMessage!(
    output::CanonicalMessage,
    message::CanonicalMessage
)
    copyto!(output.information, message.information)
    copyto!(output.precision, message.precision)

    return output
end

function dampMessage!(
    output::CanonicalMessage,
    previous::CanonicalMessage,
    current::CanonicalMessage,
    alpha::Float64
)
    oneMinusAlpha = 1.0 - alpha

    @. output.information = alpha * previous.information + oneMinusAlpha * current.information
    @. output.precision = alpha * previous.precision + oneMinusAlpha * current.precision

    copySymmetricPart!(output.precision, output.precision)

    return output
end

function canonicalMarginal(message::CanonicalMessage)
    precision = symmetricPart(message.precision)

    if !isposdef(Symmetric(precision))
        error("Marginal precision must be positive definite to compute mean and covariance.")
    end

    covariance = inverseSPDCanonical(precision)
    mean = covariance * message.information

    return CanonicalMarginal(
        copy(message.information),
        precision,
        mean,
        covariance
    )
end

function setCanonicalMarginal!(
    marginal::CanonicalMarginal,
    information,
    precision
)
    informationVector = asVector(information)
    precisionMatrix = asCovarianceMatrix(precision)

    if length(informationVector) != length(marginal.information)
        error("Information dimension does not match existing marginal dimension.")
    end

    if size(precisionMatrix, 1) != length(marginal.information)
        error("Precision row dimension does not match existing marginal dimension.")
    end

    if size(precisionMatrix, 2) != length(marginal.information)
        error("Precision column dimension does not match existing marginal dimension.")
    end

    copySymmetricPart!(marginal.precision, precisionMatrix)

    if !isposdef(Symmetric(marginal.precision))
        error("Marginal precision must be positive definite to compute mean and covariance.")
    end

    covarianceMatrix = inverseSPDCanonical(marginal.precision)

    copyto!(marginal.information, informationVector)
    copyto!(marginal.covariance, covarianceMatrix)
    mul!(marginal.mean, marginal.covariance, marginal.information)

    return marginal
end

function copyCanonicalMarginal(marginal::CanonicalMarginal)
    return CanonicalMarginal(
        copy(marginal.information),
        copy(marginal.precision),
        copy(marginal.mean),
        copy(marginal.covariance)
    )
end

"""
    GaussianCanonicalInference

Inference state for canonical-form Gaussian belief propagation.

# Fields

- `initial`: Initial variable beliefs.
- `variableToFactor`: Variable-to-factor canonical messages.
- `factorToVariable`: Factor-to-variable canonical messages.
- `marginal`: Stored canonical marginals.
- `nextVariableToFactor`: Flooding buffer for variable-to-factor messages.
- `nextFactorToVariable`: Flooding buffer for factor-to-variable messages.
- `frozenEdges`: Edge freeze flags.
- `frozenVariables`: Variable freeze flags.
- `frozenFactors`: Factor freeze flags.
- `dampedEdges`: Edge damping flags.
- `edgeDampingProb`: Per-edge damping probabilities.
- `edgeDampingAlpha`: Per-edge damping mixing weights.
- `defaultMean`: Default mean used for warm-start variable additions.
- `defaultCovariance`: Default covariance used for warm-start variable additions.
- `graphVersion`: Topology version used to detect stale inference states.

# Notes

Canonical form represents Gaussian messages as information vectors and precision
matrices. It is often convenient for algebraic message products. Create
instances with [`canonical`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = canonical(graph)
```
"""
mutable struct GaussianCanonicalInference <: GaussianSumProductInference
    initial::Vector{CanonicalMessage}
    variableToFactor::Vector{CanonicalMessage}
    factorToVariable::Vector{CanonicalMessage}
    marginal::Vector{CanonicalMarginal}
    nextVariableToFactor::Vector{CanonicalMessage}
    nextFactorToVariable::Vector{CanonicalMessage}
    frozenEdges::Vector{Bool}
    frozenVariables::Vector{Bool}
    frozenFactors::Vector{Bool}
    dampedEdges::Vector{Bool}
    edgeDampingProb::Vector{Float64}
    edgeDampingAlpha::Vector{Float64}
    defaultMean
    defaultCovariance
    graphVersion::Int
end

function defaultMeanCanonical(variable::GaussianVariable, mean)
    return defaultGaussianMean(variable, mean)
end

function defaultCovarianceCanonical(variable::GaussianVariable, covariance)
    return defaultGaussianCovariance(variable, covariance)
end

function initialCanonicalMessage(
    variable::GaussianVariable;
    mean = 0.0,
    covariance = 1e6
)
    messageMean =
        variable.mean === nothing ?
        defaultMeanCanonical(variable, mean) :
        variable.mean

    messageCovariance =
        variable.covariance === nothing ?
        defaultCovarianceCanonical(variable, covariance) :
        variable.covariance

    precision = inverseSPDCanonical(messageCovariance)
    information = precision * messageMean

    return CanonicalMessage(information, precision)
end

function initialCanonicalMessageFromUnaryFactor(factorData::GaussianFactor)
    R = inverseSPDCanonical(factorData.covariance)
    H = factorData.coefficient
    precision = symmetricPart(H' * R * H)
    information = H' * R * factorData.mean

    if !isposdef(Symmetric(precision))
        error(
            "Initializing unary GaussianFactor $(factorData.label) does not define " *
            "a positive definite precision."
        )
    end

    return CanonicalMessage(information, precision)
end

function initialCanonicalMessage(
    graph::GaussianFactorGraph,
    variableIndexValue::Int;
    mean = 0.0,
    covariance = 1e6
)
    factorData = initializingUnaryFactor(graph, variableIndexValue)

    if factorData !== nothing
        return initialCanonicalMessageFromUnaryFactor(factorData)
    end

    return initialCanonicalMessage(
        graph.variables[variableIndexValue];
        mean = mean,
        covariance = covariance
    )
end

"""
    canonical(graph::GaussianFactorGraph; mean = 0.0, covariance = 1e6)

Create canonical-form Gaussian belief propagation inference state.

# Arguments

- `graph`: Gaussian factor graph or tree view.

# Keywords

- `mean`: Default prior mean for variables without their own mean.
- `covariance`: Default prior covariance for variables without their own covariance.

# Returns

A [`GaussianCanonicalInference`](@ref) object.

# Notes

The default `mean` and `covariance` are converted to information/precision
messages for variables without explicit priors. Scalar defaults are expanded to
the dimension of each variable. These defaults are stored and reused when new
variables are added later through the warm-start [`addVariable!`](@ref) API.
A unary Gaussian factor created with `initialize = true` overrides both the
Gaussian variable prior and inference default for its connected variable.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = canonical(graph; covariance = 1e4)
```
"""
function canonical(
    graph::GaussianFactorGraph;
    mean = 0.0,
    covariance = 1e6
)
    initial = CanonicalMessage[]

    for variableIndexValue in eachindex(graph.variables)
        push!(
            initial,
            initialCanonicalMessage(
                graph,
                variableIndexValue;
                mean = mean,
                covariance = covariance
            )
        )
    end

    variableToFactor = CanonicalMessage[]
    factorToVariable = CanonicalMessage[]
    nextVariableToFactor = CanonicalMessage[]
    nextFactorToVariable = CanonicalMessage[]

    for edge in graph.edges
        message = initial[edge.variableIndex]

        push!(variableToFactor, copyCanonicalMessage(message))
        push!(factorToVariable, copyCanonicalMessage(message))

    end

    marginal = CanonicalMarginal[]

    for message in initial
        push!(marginal, canonicalMarginal(message))
    end

    return GaussianCanonicalInference(
        initial,
        variableToFactor,
        factorToVariable,
        marginal,
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

function setCanonicalInitialFromUnaryFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    factorData::GaussianFactor
)
    if !factorData.initialize
        return nothing
    end

    variableIdx = variableIndex(graph.referenceIndex, only(factorData.variables))
    message = initialCanonicalMessageFromUnaryFactor(factorData)
    marginal = canonicalMarginal(message)

    copyMessage!(inference.initial[variableIdx], message)
    setCanonicalMarginal!(
        inference.marginal[variableIdx],
        marginal.information,
        marginal.precision
    )

    for edgeId in graph.variableEdges[variableIdx]
        copyMessage!(inference.variableToFactor[edgeId], message)
        copyMessage!(inference.factorToVariable[edgeId], message)
    end

    return nothing
end

function assertGaussianCanonicalInferenceMatchesGraph(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference
)
    assertInferenceCurrent(graph, inference)

    assertVariableStorageMatchesGraph(
        graph,
        "Canonical inference object does not match the graph variables.",
        inference.initial,
        inference.marginal,
        inference.frozenVariables
    )

    assertCommonInferenceStorageMatchesGraph(
        graph,
        inference,
        "Canonical inference object does not match the graph."
    )

    return nothing
end

function assertInferenceMatchesGraph(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference
)
    return assertGaussianCanonicalInferenceMatchesGraph(graph, inference)
end

function extendGaussianCanonicalInferenceForAddedFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    firstNewEdge::Int
)
    push!(inference.frozenFactors, false)

    for edgeId in firstNewEdge:length(graph.edges)
        edge = graph.edges[edgeId]
        message = inference.initial[edge.variableIndex]

        push!(inference.variableToFactor, copyCanonicalMessage(message))
        push!(inference.factorToVariable, copyCanonicalMessage(message))

        push!(inference.frozenEdges, false)
        push!(inference.dampedEdges, false)
        push!(inference.edgeDampingProb, 0.6)
        push!(inference.edgeDampingAlpha, 0.4)
    end

    return inference
end

function extendGaussianCanonicalInferenceForAddedVariable!(
    inference::GaussianCanonicalInference,
    variable::GaussianVariable;
    mean = inference.defaultMean,
    covariance = inference.defaultCovariance
)
    message = initialCanonicalMessage(
        variable;
        mean = mean,
        covariance = covariance
    )

    push!(inference.initial, message)
    push!(inference.marginal, canonicalMarginal(message))
    push!(inference.frozenVariables, false)

    return inference
end

function extendInferenceForAddedVariable!(
    inference::GaussianCanonicalInference,
    variable::GaussianVariable;
    mean = inference.defaultMean,
    covariance = inference.defaultCovariance
)
    return extendGaussianCanonicalInferenceForAddedVariable!(
        inference,
        variable;
        mean = mean,
        covariance = covariance
    )
end

function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    factorData::GaussianFactor
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, inference)

    firstNewEdge = length(graph.edges) + 1
    added = addFactor!(graph, factorData)
    extendGaussianCanonicalInferenceForAddedFactor!(graph, inference, firstNewEdge)
    setCanonicalInitialFromUnaryFactor!(graph, inference, added)
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
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

function updateFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    factorRef::FactorRef;
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, inference)

    updated = updateFactor!(
        graph,
        factorRef;
        mean = mean,
        coefficient = coefficient,
        covariance = covariance,
        initialize = initialize
    )

    if initialize === true
        setCanonicalInitialFromUnaryFactor!(graph, inference, updated)
    end

    return updated
end

function updateFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
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


function applyFactorToVariableDamping!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    output::Vector{CanonicalMessage},
    previous::Vector{CanonicalMessage};
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

    return nothing
end

function edgeVariableCanonical(graph::GaussianFactorGraph, edgeId::Int)
    edge = graph.edges[edgeId]

    return graph.variables[edge.variableIndex]
end

function coefficientBlockViewCanonical(
    graph::GaussianFactorGraph,
    factorData::GaussianFactor,
    variableRef
)
    requestedIndex = variableIndex(graph.referenceIndex, variableRef)
    startColumn = 1

    for reference in factorData.variables
        currentIndex = variableIndex(graph.referenceIndex, reference)
        dimension = graph.variables[currentIndex].dimension
        stopColumn = startColumn + dimension - 1

        if currentIndex == requestedIndex
            return @view factorData.coefficient[:, startColumn:stopColumn]
        end

        startColumn = stopColumn + 1
    end

    error("GaussianVariable $variableRef is not connected to GaussianFactor $(factorData.label).")
end

function blockDiagonalCanonical(blocks::Vector{Matrix{Float64}})
    totalDimension = 0

    for block in blocks
        totalDimension += size(block, 1)
    end

    matrix = zeros(totalDimension, totalDimension)

    startIndex = 1

    for block in blocks
        stopIndex = startIndex + size(block, 1) - 1

        matrix[startIndex:stopIndex, startIndex:stopIndex] .= block

        startIndex = stopIndex + 1
    end

    return matrix
end

function factorToVariableMessage!(
    output::CanonicalMessage,
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    edgeId::Int
)
    edge = graph.edges[edgeId]
    factorData = graph.factors[edge.factorIndex]
    targetVariable = graph.variables[edge.variableIndex]

    R = inverseSPDCanonical(factorData.covariance)

    Htarget = coefficientBlockViewCanonical(graph, factorData, targetVariable.id)

    targetPrecision = Htarget' * R * Htarget
    targetInformation = Htarget' * R * factorData.mean

    otherBlocks = Matrix{Float64}[]
    otherPrecisions = Matrix{Float64}[]
    otherInformations = Vector{Float64}[]

    for otherEdgeId in graph.factorEdges[edge.factorIndex]
        if otherEdgeId == edgeId
            continue
        end

        otherVariable = edgeVariableCanonical(graph, otherEdgeId)
        otherMessage = inference.variableToFactor[otherEdgeId]

        Hother = coefficientBlockViewCanonical(graph, factorData, otherVariable.id)

        push!(otherBlocks, Matrix(Hother))
        push!(otherPrecisions, otherMessage.precision)
        push!(otherInformations, otherMessage.information)
    end

    if isempty(otherBlocks)
        return setCanonicalMessage!(
            output,
            targetInformation,
            symmetricPart(targetPrecision)
        )
    end

    HotherStack = hcat(otherBlocks...)

    incomingPrecision = blockDiagonalCanonical(otherPrecisions)
    incomingInformation = vcat(otherInformations...)

    otherPrecision =
        HotherStack' * R * HotherStack +
        incomingPrecision

    otherInformation =
        HotherStack' * R * factorData.mean +
        incomingInformation

    crossPrecision = Htarget' * R * HotherStack

    solvedCross = solveSPDCanonical(otherPrecision, crossPrecision')
    solvedInformation = solveSPDCanonical(otherPrecision, otherInformation)

    messagePrecision =
        targetPrecision -
        crossPrecision * solvedCross

    messageInformation =
        targetInformation -
        crossPrecision * solvedInformation

    return setCanonicalMessage!(
        output,
        messageInformation,
        symmetricPart(messagePrecision)
    )
end

function factorEdgeRangesCanonical(
    graph::GaussianFactorGraph,
    factorData::GaussianFactor,
    edgeIds::Vector{Int}
)
    rangeByVariable = Dict{Int, UnitRange{Int}}()
    startColumn = 1

    for variableRef in factorData.variables
        variableIdx = variableIndex(graph.referenceIndex, variableRef)
        dimension = graph.variables[variableIdx].dimension
        stopColumn = startColumn + dimension - 1

        rangeByVariable[variableIdx] = startColumn:stopColumn

        startColumn = stopColumn + 1
    end

    ranges = UnitRange{Int}[]

    for edgeId in edgeIds
        edge = graph.edges[edgeId]
        push!(ranges, rangeByVariable[edge.variableIndex])
    end

    return ranges
end

function concatenateRangesCanonical(ranges::Vector{UnitRange{Int}})
    indices = Int[]

    for range in ranges
        append!(indices, collect(range))
    end

    return indices
end

function updateFactorToVariableMessagesStandard!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    output::Vector{CanonicalMessage}
)
    for edge in graph.edges
        if inference.frozenEdges[edge.id] || inference.frozenFactors[edge.factorIndex]
            copyMessage!(output[edge.id], inference.factorToVariable[edge.id])
            continue
        end

        factorToVariableMessage!(
            output[edge.id],
            graph,
            inference,
            edge.id
        )
    end

    return nothing
end

function updateFactorToVariableMessagesBroadcast!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    output::Vector{CanonicalMessage}
)
    for factorIndex in eachindex(graph.factors)
        factorData = graph.factors[factorIndex]
        edgeIds = graph.factorEdges[factorIndex]

        if inference.frozenFactors[factorIndex]
            copyFactorToVariableMessages!(
                graph,
                output,
                inference.factorToVariable,
                factorIndex
            )
            continue
        end

        R = inverseSPDCanonical(factorData.covariance)
        H = factorData.coefficient

        factorPrecision = symmetricPart(H' * R * H)
        factorInformation = H' * R * factorData.mean

        ranges = factorEdgeRangesCanonical(graph, factorData, edgeIds)

        incomingPrecisions = Matrix{Float64}[]
        incomingInformations = Vector{Float64}[]

        for edgeId in edgeIds
            message = inference.variableToFactor[edgeId]

            push!(incomingPrecisions, message.precision)
            push!(incomingInformations, message.information)
        end

        for localIndex in eachindex(edgeIds)
            edgeId = edgeIds[localIndex]
            targetRange = ranges[localIndex]

            if inference.frozenEdges[edgeId]
                copyMessage!(output[edgeId], inference.factorToVariable[edgeId])
                continue
            end

            otherRanges = UnitRange{Int}[]

            for candidateIndex in eachindex(ranges)
                if candidateIndex != localIndex
                    push!(otherRanges, ranges[candidateIndex])
                end
            end

            Jii = factorPrecision[targetRange, targetRange]
            hi = factorInformation[targetRange]

            if isempty(otherRanges)
                setCanonicalMessage!(
                    output[edgeId],
                    hi,
                    Jii
                )

                continue
            end

            otherIndex = concatenateRangesCanonical(otherRanges)

            Jio = factorPrecision[targetRange, otherIndex]
            Joo = factorPrecision[otherIndex, otherIndex]
            ho = factorInformation[otherIndex]

            incomingOtherPrecisions = Matrix{Float64}[]
            incomingOtherInformations = Vector{Float64}[]

            for candidateIndex in eachindex(edgeIds)
                if candidateIndex == localIndex
                    continue
                end

                push!(incomingOtherPrecisions, incomingPrecisions[candidateIndex])
                push!(incomingOtherInformations, incomingInformations[candidateIndex])
            end

            Joo += blockDiagonalCanonical(incomingOtherPrecisions)
            ho += vcat(incomingOtherInformations...)

            solvedCross = solveSPDCanonical(Joo, Jio')
            solvedInformation = solveSPDCanonical(Joo, ho)

            messagePrecision = Jii - Jio * solvedCross
            messageInformation = hi - Jio * solvedInformation

            setCanonicalMessage!(
                output[edgeId],
                messageInformation,
                symmetricPart(messagePrecision)
            )
        end
    end

    return nothing
end

function updateFactorToVariableMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    output::Vector{CanonicalMessage};
    broadcast::Bool = false
)
    if broadcast
        return updateFactorToVariableMessagesBroadcast!(graph, inference, output)
    end

    return updateFactorToVariableMessagesStandard!(graph, inference, output)
end

function variableToFactorMessage!(
    output::CanonicalMessage,
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    edgeId::Int
)
    edge = graph.edges[edgeId]
    variableIndex = edge.variableIndex

    edgeIds = graph.variableEdges[variableIndex]
    hasIncomingMessage = false

    fill!(output.information, 0.0)
    fill!(output.precision, 0.0)

    for otherEdgeId in edgeIds
        if otherEdgeId == edgeId
            continue
        end

        output.information .+= inference.factorToVariable[otherEdgeId].information
        output.precision .+= inference.factorToVariable[otherEdgeId].precision

        hasIncomingMessage = true
    end

    if !hasIncomingMessage
        return setCanonicalMessage!(
            output,
            inference.initial[variableIndex].information,
            inference.initial[variableIndex].precision
        )
    end

    output.precision .= symmetricPart(output.precision)

    return output
end

function updateVariableToFactorMessagesStandard!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    output::Vector{CanonicalMessage}
)
    for edge in graph.edges
        if inference.frozenEdges[edge.id] || inference.frozenVariables[edge.variableIndex]
            copyMessage!(output[edge.id], inference.variableToFactor[edge.id])
            continue
        end

        variableToFactorMessage!(
            output[edge.id],
            graph,
            inference,
            edge.id
        )
    end

    return nothing
end

function updateVariableToFactorMessagesBroadcast!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    output::Vector{CanonicalMessage}
)
    for variableIndex in eachindex(graph.variables)
        edgeIds = graph.variableEdges[variableIndex]
        variable = graph.variables[variableIndex]
        dimension = variable.dimension

        if inference.frozenVariables[variableIndex]
            copyVariableToFactorMessages!(
                graph,
                output,
                inference.variableToFactor,
                variableIndex
            )
            continue
        end

        totalPrecision = zeros(dimension, dimension)
        totalInformation = zeros(dimension)

        for edgeId in edgeIds
            totalPrecision += inference.factorToVariable[edgeId].precision
            totalInformation += inference.factorToVariable[edgeId].information
        end

        for edgeId in edgeIds
            if inference.frozenEdges[edgeId]
                copyMessage!(output[edgeId], inference.variableToFactor[edgeId])
                continue
            end

            if length(edgeIds) == 1
                setCanonicalMessage!(
                    output[edgeId],
                    inference.initial[variableIndex].information,
                    inference.initial[variableIndex].precision
                )

                continue
            end

            copyto!(output[edgeId].information, totalInformation)
            output[edgeId].information .-= inference.factorToVariable[edgeId].information

            copyto!(output[edgeId].precision, totalPrecision)
            output[edgeId].precision .-= inference.factorToVariable[edgeId].precision
            copySymmetricPart!(output[edgeId].precision, output[edgeId].precision)
        end
    end

    return nothing
end

function updateVariableToFactorMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    output::Vector{CanonicalMessage};
    broadcast::Bool = false
)
    if broadcast
        return updateVariableToFactorMessagesBroadcast!(graph, inference, output)
    end

    return updateVariableToFactorMessagesStandard!(graph, inference, output)
end

function marginals!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, inference)

    for variableIndex in eachindex(graph.variables)
        edgeIds = graph.variableEdges[variableIndex]

        if isempty(edgeIds)
            setCanonicalMarginal!(
                inference.marginal[variableIndex],
                inference.initial[variableIndex].information,
                inference.initial[variableIndex].precision
            )

            continue
        end

        fill!(inference.marginal[variableIndex].information, 0.0)
        fill!(inference.marginal[variableIndex].precision, 0.0)

        for edgeId in edgeIds
            inference.marginal[variableIndex].information .+=
                inference.factorToVariable[edgeId].information

            inference.marginal[variableIndex].precision .+=
                inference.factorToVariable[edgeId].precision
        end

        setCanonicalMarginal!(
            inference.marginal[variableIndex],
            inference.marginal[variableIndex].information,
            symmetricPart(inference.marginal[variableIndex].precision)
        )
    end

    return nothing
end

function commitFactorToVariableMessages!(inference::GaussianCanonicalInference)
    inference.factorToVariable, inference.nextFactorToVariable =
        inference.nextFactorToVariable, inference.factorToVariable

    return nothing
end

function commitVariableToFactorMessages!(inference::GaussianCanonicalInference)
    inference.variableToFactor, inference.nextVariableToFactor =
        inference.nextVariableToFactor, inference.variableToFactor

    return nothing
end

function factorToVariableMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, inference)
    validateDamping(prob, alpha)

    hasDamping = damping || any(inference.dampedEdges)
    previousFactorToVariable =
        hasDamping ? copyCanonicalMessages(inference.factorToVariable) : CanonicalMessage[]

    updateFactorToVariableMessages!(
        graph,
        inference,
        inference.factorToVariable;
        broadcast = broadcast
    )

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

function variableToFactorMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
    broadcast::Bool = false
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, inference)

    updateVariableToFactorMessages!(
        graph,
        inference,
        inference.variableToFactor;
        broadcast = broadcast
    )

    return nothing
end

function messages!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    updateFraction = nothing,
    updateCount = nothing,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    assertGaussianCanonicalInferenceMatchesGraph(graph, inference)
    validateDamping(prob, alpha)
    selectedSchedule = effectiveSchedule(inference, schedule)
    validateGBPSchedule(selectedSchedule)

    if selectedSchedule == :residual
        if broadcast
            error("broadcast = true is not supported with schedule = :residual.")
        end

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
                inference.factorToVariable;
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

        commitFactorToVariableMessages!(inference)
        commitVariableToFactorMessages!(inference)

        return nothing
    end

    hasDamping = damping || any(inference.dampedEdges)
    previousFactorToVariable =
        hasDamping ? copyCanonicalMessages(inference.factorToVariable) : CanonicalMessage[]

    updateFactorToVariableMessages!(
        graph,
        inference,
        inference.factorToVariable;
        broadcast = broadcast
    )

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

    updateVariableToFactorMessages!(
        graph,
        inference,
        inference.variableToFactor;
        broadcast = broadcast
    )

    return nothing
end

function stepGBP!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
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
    marginals!(graph, inference)

    return nothing
end

function gbp!(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
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
    assertGaussianCanonicalInferenceMatchesGraph(graph, inference)

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
            [copyCanonicalMarginal(marginal) for marginal in inference.marginal]

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
                maxMarginalChange(graph, inference, previous) <= tolerance
            break
        end
    end

    return nothing
end

function marginalMean(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    variableRef::VariableRef
)
    return marginal(graph, inference, variableRef).mean
end

function marginalCovariance(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference,
    variableRef::VariableRef
)
    return marginal(graph, inference, variableRef).covariance
end
