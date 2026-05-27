"""
    WeightedLeastSquaresResult

Weighted least-squares solution returned by [`solveWLS`](@ref).

# Fields

- `mean`: Full stacked state mean.
- `covariance`: Full stacked state covariance.
- `variableMean`: Per-variable mean blocks.
- `variableCovariance`: Per-variable covariance blocks.

# Notes

Fields include the full `mean` and `covariance`, plus per-GaussianVariable
`variableMean` and `variableCovariance` blocks.

The full state is ordered in the same order as `graph.variables`. Per-GaussianVariable
blocks make it easier to compare WLS with Gaussian belief propagation marginals.

# Example

```julia
result = solveWLS(graph)
```
"""
struct WeightedLeastSquaresResult
    mean::Vector{Float64}
    covariance::Matrix{Float64}
    variableMean::Vector{Vector{Float64}}
    variableCovariance::Vector{Matrix{Float64}}
end

function inverseSPDWLS(matrix::AbstractMatrix)
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

function totalStateDimension(graph::GaussianFactorGraph)
    dimension = 0

    for variable in graph.variables
        dimension += variable.dimension
    end

    return dimension
end

function stateRanges(graph::GaussianFactorGraph)
    starts = Int[]
    stops = Int[]

    startIndex = 1

    for variable in graph.variables
        stopIndex = startIndex + variable.dimension - 1

        push!(starts, startIndex)
        push!(stops, stopIndex)

        startIndex = stopIndex + 1
    end

    return starts, stops
end

function globalCoefficientMatrix(
    graph::GaussianFactorGraph,
    factorData::GaussianFactor,
    starts::Vector{Int},
    stops::Vector{Int}
)
    meanDimension = length(factorData.mean)
    stateDimension = totalStateDimension(graph)

    Hglobal = zeros(meanDimension, stateDimension)

    localStart = 1

    for variableRef in factorData.variables
        index = variableIndex(graph.referenceIndex, variableRef)
        dimension = graph.variables[index].dimension

        localStop = localStart + dimension - 1

        localRange = localStart:localStop
        globalRange = starts[index]:stops[index]

        Hglobal[:, globalRange] .= factorData.coefficient[:, localRange]

        localStart = localStop + 1
    end

    return Hglobal
end

"""
    solveWLS(graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}})

Solve the equivalent centralized weighted least-squares problem.

# Arguments

- `graph`: Gaussian factor graph or tree view.

# Returns

A [`WeightedLeastSquaresResult`](@ref).

# Notes

This is useful as a reference solution for Gaussian belief propagation. The
solver assembles the global normal equations from all factors and returns a
[`WeightedLeastSquaresResult`](@ref). An error is thrown if the normal matrix is not positive
definite, which usually means the model is unobservable or insufficiently
constrained.

# Example

```julia
result = solveWLS(graph)
```
"""
function solveWLS(graph::GaussianFactorGraph)
    starts, stops = stateRanges(graph)

    stateDimension = totalStateDimension(graph)

    normalPrecision = zeros(stateDimension, stateDimension)
    information = zeros(stateDimension)

    for factorData in graph.factors
        H = globalCoefficientMatrix(graph, factorData, starts, stops)
        covarianceInverse = inverseSPDWLS(factorData.covariance)

        normalPrecision += H' * covarianceInverse * H
        information += H' * covarianceInverse * factorData.mean
    end

    normalPrecision = symmetricPart(normalPrecision)

    if !isposdef(Symmetric(normalPrecision))
        error(
            "The WLS normal matrix is not positive definite. " *
            "The model may be unobservable or insufficiently constrained."
        )
    end

    covariance = inverseSPDWLS(normalPrecision)
    mean = covariance * information

    variableMean = Vector{Float64}[]
    variableCovariance = Matrix{Float64}[]

    for index in eachindex(graph.variables)
        range = starts[index]:stops[index]

        push!(variableMean, mean[range])
        push!(variableCovariance, covariance[range, range])
    end

    return WeightedLeastSquaresResult(mean, covariance, variableMean, variableCovariance)
end

function solveWLS(tree::TreeFactorGraph{GaussianFactorGraph})
    return solveWLS(tree.graph)
end

"""
    printWLS(
        graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}},
        result::WeightedLeastSquaresResult
    )

Print a weighted least-squares result.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `result`: Weighted least-squares result.

# Notes

Prints one mean and covariance block per GaussianVariable.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
result = solveWLS(graph)

printWLS(graph, result)
```
"""
function printWLS(graph::GaussianFactorGraph, result::WeightedLeastSquaresResult)
    println("WLS result")

    for (index, variable) in pairs(graph.variables)
        println("Estimate for GaussianVariable ", variable.label)
        println("  mean = ", result.variableMean[index])
        println("  covariance = ")
        show(stdout, "text/plain", result.variableCovariance[index])
        println()
        println()
    end

    return nothing
end

function printWLS(tree::TreeFactorGraph{GaussianFactorGraph}, result::WeightedLeastSquaresResult)
    return printWLS(tree.graph, result)
end

function factorMarginalMean(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    factorData::GaussianFactor
)
    values = Float64[]

    for variableRef in factorData.variables
        variableIdx = variableIndex(graph.referenceIndex, variableRef)
        append!(values, inference.marginal[variableIdx].mean)
    end

    return values
end

function factorResidual(factorData::GaussianFactor, stateMean::Vector{Float64})
    return factorData.mean - factorData.coefficient * stateMean
end

function residualRows(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference;
    normalize::Bool
)
    assertInferenceMatchesGraph(graph, inference)

    rows = Vector{NamedTuple{(:factor, :index, :value), Tuple{String, Int, Vector{Float64}}}}()

    for (index, factorData) in pairs(graph.factors)
        stateMean = factorMarginalMean(graph, inference, factorData)
        value = factorResidual(factorData, stateMean)

        if normalize
            value = value ./ sqrt.(diag(factorData.covariance))
        end

        push!(
            rows,
            (
                factor = factorData.label,
                index = index,
                value = value
            )
        )
    end

    return rows
end

"""
    residuals(
        graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}},
        inference::GaussianSumProductInference
    )

Return one residual vector per GaussianFactor using the current marginal means.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching Gaussian inference object with stored marginal means.

# Returns

A vector of named tuples with `factor`, `index`, and `value` fields.

# Notes

For each GaussianFactor this computes `GaussianFactor.mean - GaussianFactor.coefficient * x`, where
`x` is the stacked vector of current marginal means for the variables connected
to that GaussianFactor.

# Example

```julia
inference = moment(graph)
gbp!(graph, inference)

residuals(graph, inference)
```
"""
function residuals(graph::GaussianFactorGraph, inference::GaussianSumProductInference)
    return residualRows(graph, inference; normalize = false)
end

function residuals(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference
)
    return residuals(tree.graph, inference)
end

"""
    normalizedResiduals(
        graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}},
        inference::GaussianSumProductInference
    )

Return residuals scaled by each GaussianFactor component's standard deviation.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching Gaussian inference object with stored marginal means.

# Returns

A vector of named tuples with `factor`, `index`, and normalized `value` fields.

# Notes

This uses the diagonal entries of the GaussianFactor covariance matrix, so each
component is computed as `residual / sqrt(covariance_component)`. The return
format matches [`residuals`](@ref).

# Example

```julia
inference = moment(graph)
gbp!(graph, inference)

normalizedResiduals(graph, inference)
```
"""
function normalizedResiduals(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference
)
    return residualRows(graph, inference; normalize = true)
end

function normalizedResiduals(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference
)
    return normalizedResiduals(tree.graph, inference)
end

"""
    compareMeanWithWLS(
        graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}},
        inference::GaussianSumProductInference, result::WeightedLeastSquaresResult
    )

Print per-GaussianVariable mean differences between Gaussian belief propagation marginals and WLS.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching Gaussian inference object.
- `result`: Weighted least-squares result.

# Notes

Use [`maxMeanError`](@ref) when a numeric error value is needed instead of
formatted output.

# Example

```julia
inference = moment(graph)
gbp!(graph, inference)

result = solveWLS(graph)
compareMeanWithWLS(graph, inference, result)
```
"""
function compareMeanWithWLS(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    assertInferenceMatchesGraph(graph, inference)

    println("Gaussian belief propagation and WLS comparison")

    for (index, variable) in pairs(graph.variables)
        gbpMean = inference.marginal[index].mean
        wlsMean = result.variableMean[index]

        println("GaussianVariable ", variable.label)
        println("  Gaussian belief propagation mean = ", gbpMean)
        println("  WLS mean = ", wlsMean)
        println("  error    = ", gbpMean - wlsMean)
        println("  norm     = ", norm(gbpMean - wlsMean))
        println()
    end

    return nothing
end

function compareMeanWithWLS(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    return compareMeanWithWLS(tree.graph, inference, result)
end

"""
    maxMeanError(
        graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}},
        inference::GaussianSumProductInference, result::WeightedLeastSquaresResult
    )

Return the maximum norm of per-GaussianVariable mean errors against WLS.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching Gaussian inference object.
- `result`: Weighted least-squares result.

# Returns

The maximum Euclidean mean error across variables.

# Notes

The returned scalar is `maximum(norm(gbpMean - wlsMean))` across all variables.

# Example

```julia
inference = moment(graph)
gbp!(graph, inference)

result = solveWLS(graph)
maxMeanError(graph, inference, result)
```
"""
function maxMeanError(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    assertInferenceMatchesGraph(graph, inference)

    maximumError = 0.0

    for index in eachindex(graph.variables)
        gbpMean = inference.marginal[index].mean
        wlsMean = result.variableMean[index]

        maximumError = max(maximumError, norm(gbpMean - wlsMean))
    end

    return maximumError
end

function maxMeanError(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    return maxMeanError(tree.graph, inference, result)
end

"""
    compareCovarianceWithWLS(
        graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}},
        inference::GaussianSumProductInference, result::WeightedLeastSquaresResult
    )

Print per-GaussianVariable covariance differences between Gaussian belief propagation
marginals and WLS.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching Gaussian inference object.
- `result`: Weighted least-squares result.

# Notes

Use [`maxCovarianceError`](@ref) when a numeric error value is needed.

# Example

```julia
inference = moment(graph)
gbp!(graph, inference)

result = solveWLS(graph)
compareCovarianceWithWLS(graph, inference, result)
```
"""
function compareCovarianceWithWLS(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    assertInferenceMatchesGraph(graph, inference)

    println("Gaussian belief propagation and WLS covariance comparison")

    for (index, variable) in pairs(graph.variables)
        gbpCovariance = inference.marginal[index].covariance
        wlsCovariance = result.variableCovariance[index]

        println("GaussianVariable ", variable.label)
        println("  Gaussian belief propagation covariance = ")
        display(gbpCovariance)
        println("  WLS covariance = ")
        display(wlsCovariance)
        println("  error norm = ", norm(gbpCovariance - wlsCovariance))
        println()
    end

    return nothing
end

function compareCovarianceWithWLS(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    return compareCovarianceWithWLS(tree.graph, inference, result)
end

"""
    maxCovarianceError(
        graph::Union{GaussianFactorGraph, TreeFactorGraph{GaussianFactorGraph}},
        inference::GaussianSumProductInference, result::WeightedLeastSquaresResult
    )

Return the maximum norm of per-GaussianVariable covariance errors against WLS.

# Arguments

- `graph`: Gaussian factor graph or tree view.
- `inference`: Matching Gaussian inference object.
- `result`: Weighted least-squares result.

# Returns

The maximum covariance-matrix norm error across variables.

# Notes

The returned scalar is `maximum(norm(gbpCovariance - wlsCovariance))` across all
variables.

# Example

```julia
inference = moment(graph)
gbp!(graph, inference)

result = solveWLS(graph)
maxCovarianceError(graph, inference, result)
```
"""
function maxCovarianceError(
    graph::GaussianFactorGraph,
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    assertInferenceMatchesGraph(graph, inference)

    maximumError = 0.0

    for index in eachindex(graph.variables)
        gbpCovariance = inference.marginal[index].covariance
        wlsCovariance = result.variableCovariance[index]

        maximumError = max(maximumError, norm(gbpCovariance - wlsCovariance))
    end

    return maximumError
end

function maxCovarianceError(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianSumProductInference,
    result::WeightedLeastSquaresResult
)
    return maxCovarianceError(tree.graph, inference, result)
end
