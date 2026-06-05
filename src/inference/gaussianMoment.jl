function inverseFromPositiveDefinite(matrix::AbstractMatrix)
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

mutable struct GaussianMessage
    mean::Vector{Float64}
    covariance::Matrix{Float64}
    precision::Matrix{Float64}
    information::Vector{Float64}

    function GaussianMessage(mean, covariance)
        meanVector = asVector(mean)
        covarianceMatrix = symmetricPart(asCovarianceMatrix(covariance))

        if size(covarianceMatrix, 1) != length(meanVector)
            error("Message covariance row dimension must match mean dimension.")
        end

        if size(covarianceMatrix, 2) != length(meanVector)
            error("Message covariance column dimension must match mean dimension.")
        end

        if !isposdef(Symmetric(covarianceMatrix))
            error("Message covariance must be positive definite.")
        end

        precisionMatrix = inverseFromPositiveDefinite(covarianceMatrix)
        informationVector = precisionMatrix * meanVector

        return new(
            meanVector,
            covarianceMatrix,
            precisionMatrix,
            informationVector
        )
    end
end

function setGaussianMessage!(
    message::GaussianMessage,
    mean,
    covariance
)
    meanVector = asVector(mean)
    covarianceMatrix = symmetricPart(asCovarianceMatrix(covariance))

    if length(meanVector) != length(message.mean)
        error("Mean dimension does not match existing message dimension.")
    end

    if size(covarianceMatrix, 1) != length(message.mean)
        error("Covariance row dimension does not match existing message dimension.")
    end

    if size(covarianceMatrix, 2) != length(message.mean)
        error("Covariance column dimension does not match existing message dimension.")
    end

    if !isposdef(Symmetric(covarianceMatrix))
        error("Message covariance must be positive definite.")
    end

    precisionMatrix = inverseFromPositiveDefinite(covarianceMatrix)

    copyto!(message.mean, meanVector)
    copyto!(message.covariance, covarianceMatrix)
    copyto!(message.precision, precisionMatrix)
    mul!(message.information, message.precision, message.mean)

    return message
end

function copyGaussianMessage(message::GaussianMessage)
    return GaussianMessage(
        copy(message.mean),
        copy(message.covariance)
    )
end

function copyGaussianMessages(messages::Vector{GaussianMessage})
    return copyMessages(messages, copyGaussianMessage)
end

function copyMessage!(
    output::GaussianMessage,
    message::GaussianMessage
)
    copyto!(output.mean, message.mean)
    copyto!(output.covariance, message.covariance)
    copyto!(output.precision, message.precision)
    copyto!(output.information, message.information)

    return output
end

function dampMessage!(
    output::GaussianMessage,
    previous::GaussianMessage,
    current::GaussianMessage,
    alpha::Float64
)
    oneMinusAlpha = 1.0 - alpha

    @. output.mean = alpha * previous.mean + oneMinusAlpha * current.mean
    @. output.covariance = alpha * previous.covariance + oneMinusAlpha * current.covariance

    precisionMatrix = inverseFromPositiveDefinite(output.covariance)
    copyto!(output.precision, precisionMatrix)
    mul!(output.information, output.precision, output.mean)

    return output
end

"""
    GaussianMomentInference

Inference state for moment-form Gaussian belief propagation.

# Fields

- `initial`: Initial variable beliefs.
- `variableToFactor`: Variable-to-factor messages.
- `factorToVariable`: Factor-to-variable messages.
- `marginal`: Stored variable marginals.
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

Moment form stores Gaussian messages directly as means and covariances, with
precision and information cached for efficient multiplication of messages.
Create instances with [`moment`](@ref).

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)
```
"""
mutable struct GaussianMomentInference <: GaussianSumProductInference
    initial::Vector{GaussianMessage}
    variableToFactor::Vector{GaussianMessage}
    factorToVariable::Vector{GaussianMessage}
    marginal::Vector{GaussianMessage}
    nextVariableToFactor::Vector{GaussianMessage}
    nextFactorToVariable::Vector{GaussianMessage}
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

function defaultMean(variable::GaussianVariable, mean)
    return defaultGaussianMean(variable, mean)
end

function defaultCovariance(variable::GaussianVariable, covariance)
    return defaultGaussianCovariance(variable, covariance)
end

function initialMessage(
    variable::GaussianVariable;
    mean = 0.0,
    covariance = 1e6
)
    messageMean =
        variable.mean === nothing ?
        defaultMean(variable, mean) :
        variable.mean

    messageCovariance =
        variable.covariance === nothing ?
        defaultCovariance(variable, covariance) :
        variable.covariance

    return GaussianMessage(messageMean, messageCovariance)
end

function initialMessageFromUnaryFactor(factorData::GaussianFactor)
    R = inverseFromPositiveDefinite(factorData.covariance)
    H = factorData.coefficient
    precision = symmetricPart(H' * R * H)
    information = H' * R * factorData.mean
    covariance = inverseFromPositiveDefinite(precision)
    mean = covariance * information

    return GaussianMessage(mean, covariance)
end

function initialMessage(
    graph::GaussianFactorGraph,
    variableIndexValue::Int;
    mean = 0.0,
    covariance = 1e6
)
    factorData = initializingUnaryFactor(graph, variableIndexValue)

    if factorData !== nothing
        return initialMessageFromUnaryFactor(factorData)
    end

    return initialMessage(
        graph.variables[variableIndexValue];
        mean = mean,
        covariance = covariance
    )
end

"""
    moment(graph::GaussianFactorGraph; mean = 0.0, covariance = 1e6)

Create moment-form Gaussian belief propagation inference state.

# Arguments

- `graph`: Gaussian factor graph or tree view.

# Keywords

- `mean`: Default prior mean for variables without their own mean.
- `covariance`: Default prior covariance for variables without their own covariance.

# Returns

A [`GaussianMomentInference`](@ref) object.

# Notes

The returned object contains all messages and marginals needed by [`gbp!`](@ref).
Variables that were constructed without explicit priors use the provided
default `mean` and `covariance`. Scalar defaults are expanded to the dimension
of each variable. These defaults are stored and reused when new variables are
added later through the warm-start [`addVariable!`](@ref) API. A unary Gaussian
factor created with `initialize = true` overrides both the Gaussian variable
prior and inference default for its connected variable.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph; covariance = 1e4)
```
"""
function moment(
    graph::GaussianFactorGraph;
    mean = 0.0,
    covariance = 1e6
)
    initial = GaussianMessage[]

    for variableIndexValue in eachindex(graph.variables)
        push!(
            initial,
            initialMessage(
                graph,
                variableIndexValue;
                mean = mean,
                covariance = covariance
            )
        )
    end

    variableToFactor = GaussianMessage[]
    factorToVariable = GaussianMessage[]
    nextVariableToFactor = GaussianMessage[]
    nextFactorToVariable = GaussianMessage[]

    for edge in graph.edges
        message = initial[edge.variableIndex]

        push!(variableToFactor, copyGaussianMessage(message))
        push!(factorToVariable, copyGaussianMessage(message))

    end

    marginal = GaussianMessage[]

    for message in initial
        push!(marginal, copyGaussianMessage(message))
    end

    return GaussianMomentInference(
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

function setMomentInitialFromUnaryFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    factorData::GaussianFactor
)
    if !factorData.initialize
        return nothing
    end

    variableIdx = variableIndex(graph.referenceIndex, only(factorData.variables))
    message = initialMessageFromUnaryFactor(factorData)

    copyMessage!(inference.initial[variableIdx], message)
    copyMessage!(inference.marginal[variableIdx], message)

    for edgeId in graph.variableEdges[variableIdx]
        copyMessage!(inference.variableToFactor[edgeId], message)
        copyMessage!(inference.factorToVariable[edgeId], message)
    end

    return nothing
end

function assertGaussianMomentInferenceMatchesGraph(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference
)
    assertInferenceCurrent(graph, inference)

    assertVariableStorageMatchesGraph(
        graph,
        "Moment inference object does not match the graph variables.",
        inference.initial,
        inference.marginal,
        inference.frozenVariables
    )

    assertCommonInferenceStorageMatchesGraph(
        graph,
        inference,
        "Moment inference object does not match the graph."
    )

    return nothing
end

function assertInferenceMatchesGraph(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference
)
    return assertGaussianMomentInferenceMatchesGraph(graph, inference)
end

function extendGaussianMomentInferenceForAddedFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    firstNewEdge::Int
)
    push!(inference.frozenFactors, false)

    for edgeId in firstNewEdge:length(graph.edges)
        edge = graph.edges[edgeId]
        message = inference.initial[edge.variableIndex]

        push!(inference.variableToFactor, copyGaussianMessage(message))
        push!(inference.factorToVariable, copyGaussianMessage(message))

        push!(inference.frozenEdges, false)
        push!(inference.dampedEdges, false)
        push!(inference.edgeDampingProb, 0.6)
        push!(inference.edgeDampingAlpha, 0.4)
    end

    return inference
end

function extendGaussianMomentInferenceForAddedVariable!(
    inference::GaussianMomentInference,
    variable::GaussianVariable;
    mean = inference.defaultMean,
    covariance = inference.defaultCovariance
)
    message = initialMessage(
        variable;
        mean = mean,
        covariance = covariance
    )

    push!(inference.initial, message)
    push!(inference.marginal, copyGaussianMessage(message))
    push!(inference.frozenVariables, false)

    return inference
end

function extendInferenceForAddedVariable!(
    inference::GaussianMomentInference,
    variable::GaussianVariable;
    mean = inference.defaultMean,
    covariance = inference.defaultCovariance
)
    return extendGaussianMomentInferenceForAddedVariable!(
        inference,
        variable;
        mean = mean,
        covariance = covariance
    )
end

function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    factorData::GaussianFactor
)
    assertGaussianMomentInferenceMatchesGraph(graph, inference)

    firstNewEdge = length(graph.edges) + 1
    added = addFactor!(graph, factorData)
    extendGaussianMomentInferenceForAddedFactor!(graph, inference, firstNewEdge)
    setMomentInitialFromUnaryFactor!(graph, inference, added)
    setInferenceGraphVersion!(graph, inference)

    return added
end

function addFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
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
    inference::GaussianMomentInference,
    factorRef::FactorRef;
    mean = nothing,
    coefficient = nothing,
    covariance = nothing,
    initialize::Union{Nothing, Bool} = nothing
)
    assertGaussianMomentInferenceMatchesGraph(graph, inference)

    updated = updateFactor!(
        graph,
        factorRef;
        mean = mean,
        coefficient = coefficient,
        covariance = covariance,
        initialize = initialize
    )

    if initialize === true
        setMomentInitialFromUnaryFactor!(graph, inference, updated)
    end

    return updated
end

function updateFactor!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference;
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

function precisionInformation(message::GaussianMessage)
    return message.precision, message.information
end


function applyFactorToVariableDamping!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    output::Vector{GaussianMessage},
    previous::Vector{GaussianMessage};
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

function gaussianFromPrecisionInformation!(
    message::GaussianMessage,
    precision::AbstractMatrix,
    information::AbstractVector
)
    copySymmetricPart!(message.precision, precision)

    if size(message.precision, 1) != size(message.precision, 2)
        error("Precision matrix must be square.")
    end

    if size(message.precision, 1) != length(information)
        error("Information vector dimension must match precision matrix dimension.")
    end

    if !isposdef(Symmetric(message.precision))
        error("Precision matrix must be positive definite.")
    end

    covarianceMatrix = inverseFromPositiveDefinite(message.precision)

    copyto!(message.covariance, covarianceMatrix)
    copyto!(message.information, information)
    mul!(message.mean, message.covariance, message.information)

    return message
end

function coefficientBlockView(
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

function edgeVariable(graph::GaussianFactorGraph, edgeId::Int)
    edge = graph.edges[edgeId]

    return graph.variables[edge.variableIndex]
end

function factorToVariableMessage!(
    output::GaussianMessage,
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    edgeId::Int
)
    edge = graph.edges[edgeId]
    factorData = graph.factors[edge.factorIndex]
    targetVariable = graph.variables[edge.variableIndex]

    Htarget = coefficientBlockView(graph, factorData, targetVariable.id)

    residual = copy(factorData.mean)
    innovationCovariance = copy(factorData.covariance)

    for otherEdgeId in graph.factorEdges[edge.factorIndex]
        if otherEdgeId == edgeId
            continue
        end

        otherVariable = edgeVariable(graph, otherEdgeId)
        otherMessage = inference.variableToFactor[otherEdgeId]
        Hother = coefficientBlockView(graph, factorData, otherVariable.id)

        residual -= Hother * otherMessage.mean
        innovationCovariance += Hother * otherMessage.covariance * Hother'
    end

    innovationCovariance = symmetricPart(innovationCovariance)

    if !isposdef(Symmetric(innovationCovariance))
        error("Innovation covariance is not positive definite for GaussianFactor $(factorData.label).")
    end

    solvedH = Symmetric(innovationCovariance) \ Htarget
    solvedResidual = Symmetric(innovationCovariance) \ residual

    precision = symmetricPart(Htarget' * solvedH)
    information = Htarget' * solvedResidual

    if !isposdef(Symmetric(precision))
        error(
            "GaussianFactor-to-GaussianVariable precision is not positive definite for message " *
            "$(factorData.label) -> $(targetVariable.label). " *
            "Moment form requires a proper Gaussian with finite covariance."
        )
    end

    return gaussianFromPrecisionInformation!(output, precision, information)
end

function updateFactorToVariableMessagesStandard!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    output::Vector{GaussianMessage}
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
    inference::GaussianMomentInference,
    output::Vector{GaussianMessage}
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

        totalResidual = copy(factorData.mean)
        totalInnovationCovariance = copy(factorData.covariance)

        meanContributions = Vector{Vector{Float64}}(undef, length(edgeIds))
        covarianceContributions = Vector{Matrix{Float64}}(undef, length(edgeIds))

        for localIndex in eachindex(edgeIds)
            edgeId = edgeIds[localIndex]
            variable = edgeVariable(graph, edgeId)
            message = inference.variableToFactor[edgeId]

            H = coefficientBlockView(graph, factorData, variable.id)

            meanContribution = H * message.mean
            covarianceContribution = H * message.covariance * H'

            meanContributions[localIndex] = meanContribution
            covarianceContributions[localIndex] = symmetricPart(covarianceContribution)

            totalResidual -= meanContribution
            totalInnovationCovariance += covarianceContribution
        end

        totalInnovationCovariance = symmetricPart(totalInnovationCovariance)

        for localIndex in eachindex(edgeIds)
            edgeId = edgeIds[localIndex]
            edge = graph.edges[edgeId]
            targetVariable = graph.variables[edge.variableIndex]

            if inference.frozenEdges[edgeId]
                copyMessage!(output[edgeId], inference.factorToVariable[edgeId])
                continue
            end

            Htarget = coefficientBlockView(graph, factorData, targetVariable.id)

            residual = totalResidual + meanContributions[localIndex]
            innovationCovariance =
                totalInnovationCovariance - covarianceContributions[localIndex]

            innovationCovariance = symmetricPart(innovationCovariance)

            if !isposdef(Symmetric(innovationCovariance))
                error(
                    "Innovation covariance is not positive definite for message " *
                    "$(factorData.label) -> $(targetVariable.label)."
                )
            end

            solvedH = Symmetric(innovationCovariance) \ Htarget
            solvedResidual = Symmetric(innovationCovariance) \ residual

            precision = symmetricPart(Htarget' * solvedH)
            information = Htarget' * solvedResidual

            if !isposdef(Symmetric(precision))
                error(
                    "GaussianFactor-to-GaussianVariable precision is not positive definite for message " *
                    "$(factorData.label) -> $(targetVariable.label). " *
                    "Moment form requires a proper Gaussian with finite covariance."
                )
            end

            gaussianFromPrecisionInformation!(
                output[edgeId],
                precision,
                information
            )
        end
    end

    return nothing
end

function updateFactorToVariableMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    output::Vector{GaussianMessage};
    broadcast::Bool = false
)
    if broadcast
        return updateFactorToVariableMessagesBroadcast!(graph, inference, output)
    end

    return updateFactorToVariableMessagesStandard!(graph, inference, output)
end

function variableToFactorMessage!(
    output::GaussianMessage,
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    edgeId::Int
)
    edge = graph.edges[edgeId]
    variable = graph.variables[edge.variableIndex]

    fill!(output.precision, 0.0)
    fill!(output.information, 0.0)
    hasIncomingMessage = false

    for otherEdgeId in graph.variableEdges[edge.variableIndex]
        if otherEdgeId == edgeId
            continue
        end

        output.precision .+= inference.factorToVariable[otherEdgeId].precision
        output.information .+= inference.factorToVariable[otherEdgeId].information
        hasIncomingMessage = true
    end

    if !hasIncomingMessage
        return setGaussianMessage!(
            output,
            inference.initial[edge.variableIndex].mean,
            inference.initial[edge.variableIndex].covariance
        )
    end

    return gaussianFromPrecisionInformation!(output, output.precision, output.information)
end

function updateVariableToFactorMessagesStandard!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    output::Vector{GaussianMessage}
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
    inference::GaussianMomentInference,
    output::Vector{GaussianMessage}
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

        edgePrecisions = Vector{Matrix{Float64}}(undef, length(edgeIds))
        edgeInformations = Vector{Vector{Float64}}(undef, length(edgeIds))

        for localIndex in eachindex(edgeIds)
            edgeId = edgeIds[localIndex]
            precision, information = precisionInformation(inference.factorToVariable[edgeId])

            edgePrecisions[localIndex] = precision
            edgeInformations[localIndex] = information

            totalPrecision += precision
            totalInformation += information
        end

        for localIndex in eachindex(edgeIds)
            edgeId = edgeIds[localIndex]

            if inference.frozenEdges[edgeId]
                copyMessage!(output[edgeId], inference.variableToFactor[edgeId])
                continue
            end

            if length(edgeIds) == 1
                setGaussianMessage!(
                    output[edgeId],
                    inference.initial[variableIndex].mean,
                    inference.initial[variableIndex].covariance
                )

                continue
            end

            outgoingPrecision =
                symmetricPart(totalPrecision - edgePrecisions[localIndex])

            outgoingInformation =
                totalInformation - edgeInformations[localIndex]

            gaussianFromPrecisionInformation!(
                output[edgeId],
                outgoingPrecision,
                outgoingInformation
            )
        end
    end

    return nothing
end

function updateVariableToFactorMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    output::Vector{GaussianMessage};
    broadcast::Bool = false
)
    if broadcast
        return updateVariableToFactorMessagesBroadcast!(graph, inference, output)
    end

    return updateVariableToFactorMessagesStandard!(graph, inference, output)
end

function marginals!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference
)
    assertGaussianMomentInferenceMatchesGraph(graph, inference)

    for variableIndex in eachindex(graph.variables)
        edgeIds = graph.variableEdges[variableIndex]
        variable = graph.variables[variableIndex]
        dimension = variable.dimension

        if isempty(edgeIds)
            setGaussianMessage!(
                inference.marginal[variableIndex],
                inference.initial[variableIndex].mean,
                inference.initial[variableIndex].covariance
            )

            continue
        end

        fill!(inference.marginal[variableIndex].precision, 0.0)
        fill!(inference.marginal[variableIndex].information, 0.0)

        for edgeId in edgeIds
            inference.marginal[variableIndex].precision .+=
                inference.factorToVariable[edgeId].precision

            inference.marginal[variableIndex].information .+=
                inference.factorToVariable[edgeId].information
        end

        gaussianFromPrecisionInformation!(
            inference.marginal[variableIndex],
            inference.marginal[variableIndex].precision,
            inference.marginal[variableIndex].information
        )
    end

    return nothing
end

function commitFactorToVariableMessages!(inference::GaussianMomentInference)
    inference.factorToVariable, inference.nextFactorToVariable =
        inference.nextFactorToVariable, inference.factorToVariable

    return nothing
end

function commitVariableToFactorMessages!(inference::GaussianMomentInference)
    inference.variableToFactor, inference.nextVariableToFactor =
        inference.nextVariableToFactor, inference.variableToFactor

    return nothing
end

function factorToVariableMessages!(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    assertGaussianMomentInferenceMatchesGraph(graph, inference)
    validateDamping(prob, alpha)

    hasDamping = damping || any(inference.dampedEdges)
    previousFactorToVariable =
        hasDamping ? copyGaussianMessages(inference.factorToVariable) : GaussianMessage[]

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
    inference::GaussianMomentInference;
    broadcast::Bool = false
)
    assertGaussianMomentInferenceMatchesGraph(graph, inference)

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
    inference::GaussianMomentInference;
    damping::Bool = false,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    schedule = nothing,
    updateFraction = nothing,
    updateCount = nothing,
    broadcast::Bool = false,
    rng = Random.default_rng()
)
    assertGaussianMomentInferenceMatchesGraph(graph, inference)
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
        hasDamping ? copyGaussianMessages(inference.factorToVariable) : GaussianMessage[]

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
    inference::GaussianMomentInference;
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
    inference::GaussianMomentInference;
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
    assertGaussianMomentInferenceMatchesGraph(graph, inference)

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
            copyGaussianMessages(inference.marginal)

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

"""
    marginalMean(
        graph::GaussianFactorGraph, inference::GaussianSumProductInference,
        variable::VariableRef
    )

Return the current marginal mean for one variable.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching moment or canonical inference object.
- `variable`: Variable ID or label.

# Returns

The marginal mean vector.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

marginalMean(graph, inference, :x1)
```
"""
function marginalMean(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    variableRef::VariableRef
)
    return marginal(graph, inference, variableRef).mean
end

function marginalMean(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    variableRef::VariableRef
)
    return marginal(graph, inference, variableRef).mean
end

"""
    marginalCovariance(
        graph::GaussianFactorGraph, inference::GaussianSumProductInference,
        variable::VariableRef
    )

Return the current marginal covariance for one variable.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching moment or canonical inference object.
- `variable`: Variable ID or label.

# Returns

The marginal covariance matrix.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

marginalCovariance(graph, inference, :x1)
```
"""
function marginalCovariance(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    variableRef::VariableRef
)
    return marginal(graph, inference, variableRef).covariance
end

function marginalCovariance(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference,
    variableRef::VariableRef
)
    return marginal(graph, inference, variableRef).covariance
end
