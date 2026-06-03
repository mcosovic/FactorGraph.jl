function graphVariables(graph::AbstractFactorGraph)
    return graph.variables
end

function graphFactors(graph::AbstractFactorGraph)
    return graph.factors
end

function graphReferenceIndex(graph::AbstractFactorGraph)
    return graph.referenceIndex
end

"""
    printGraph(graph::AbstractFactorGraph)

Print variables, factors, and edge connectivity as tables.

# Arguments

- `graph`: Gaussian or discrete factor graph, or tree view.

# Notes

The output shows each variable with connected factor labels and edge IDs, each
factor with connected variable labels and edge IDs, and finally the edge list.
This function is intended for interactive inspection.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])

printGraph(graph)
```
"""
function printGraph(graph::GaussianFactorGraph)
    println("Factor graph")
    println("Number of variable nodes: ", length(graphVariables(graph)))
    println("Number of factor nodes: ", length(graphFactors(graph)))
    println("Number of edges: ", length(graph.edges))
    println()

    variableRows = graphVariableRows(graph)
    factorRows = graphFactorRows(graph)
    edgeRows = graphEdgeRows(graph)

    printTable(["Variable", "Dimension", "Components", "Factors", "Edges"], variableRows)

    println()
    printTable(["Factor", "Variables", "Edges"], factorRows)

    println()
    printTable(["Edge", "Variable", "Factor"], edgeRows)

    return nothing
end

function printGraph(tree::TreeFactorGraph{GaussianFactorGraph})
    return printGraph(tree.graph)
end

function printGraph(graph::DiscreteFactorGraph)
    println("Discrete factor graph")
    println("Number of variable nodes: ", length(graphVariables(graph)))
    println("Number of factor nodes: ", length(graphFactors(graph)))
    println("Number of edges: ", length(graph.edges))
    println()

    variableRows = graphVariableRows(graph)
    factorRows = graphFactorRows(graph)
    edgeRows = graphEdgeRows(graph)

    printTable(["Variable", "Cardinality", "States", "Factors", "Edges"], variableRows)

    println()
    printTable(["Factor", "Variables", "Edges"], factorRows)

    println()
    printTable(["Edge", "Variable", "Factor"], edgeRows)

    return nothing
end

function printGraph(tree::TreeFactorGraph{DiscreteFactorGraph})
    return printGraph(tree.graph)
end

function graphVariableRows(graph::GaussianFactorGraph)
    rows = Vector{Vector{String}}()

    for (variableIndexValue, variable) in pairs(graphVariables(graph))
        factorLabels = String[]

        for edgeId in graph.variableEdges[variableIndexValue]
            edge = graph.edges[edgeId]
            factorData = graphFactors(graph)[edge.factorIndex]
            push!(factorLabels, factorData.label)
        end

        push!(
            rows,
            [
                variable.label,
                string(variable.dimension),
                formattedRefs(variable.components),
                join(factorLabels, ", "),
                join(graph.variableEdges[variableIndexValue], ", ")
            ]
        )
    end

    return rows
end

function graphVariableRows(graph::DiscreteFactorGraph)
    rows = Vector{Vector{String}}()

    for (variableIndexValue, variable) in pairs(graphVariables(graph))
        factorLabels = String[]

        for edgeId in graph.variableEdges[variableIndexValue]
            edge = graph.edges[edgeId]
            factorData = graphFactors(graph)[edge.factorIndex]
            push!(factorLabels, factorData.label)
        end

        push!(
            rows,
            [
                variable.label,
                string(variable.cardinality),
                formattedRefs(variable.states),
                join(factorLabels, ", "),
                join(graph.variableEdges[variableIndexValue], ", ")
            ]
        )
    end

    return rows
end

function graphFactorRows(graph::AbstractFactorGraph)
    rows = Vector{Vector{String}}()

    for (factorIndexValue, factorData) in pairs(graphFactors(graph))
        variableLabels = String[]

        for variableRef in factorData.variables
            variableIndexValue = variableIndex(graphReferenceIndex(graph), variableRef)
            push!(variableLabels, graphVariables(graph)[variableIndexValue].label)
        end

        push!(
            rows,
            [
                factorData.label,
                join(variableLabels, ", "),
                join(graph.factorEdges[factorIndexValue], ", ")
            ]
        )
    end

    return rows
end

function graphEdgeRows(graph::AbstractFactorGraph)
    rows = Vector{Vector{String}}()

    for edge in graph.edges
        factorData = graphFactors(graph)[edge.factorIndex]
        variable = graphVariables(graph)[edge.variableIndex]

        push!(
            rows,
            [
                string(edge.id),
                variable.label,
                factorData.label
            ]
        )
    end

    return rows
end

function columnWidths(rows::Vector{Vector{String}}, headers::Vector{String})
    columnCount = length(headers)
    widths = zeros(Int, columnCount)

    for column in eachindex(headers)
        widths[column] = length(headers[column])
    end

    for row in rows
        for column in eachindex(row)
            widths[column] = max(widths[column], length(row[column]))
        end
    end

    return widths
end

function printTable(headers::Vector{String}, rows::Vector{Vector{String}})
    widths = columnWidths(rows, headers)

    printTableBorder(widths)
    printTableRow(headers, widths)
    printTableBorder(widths)

    for row in rows
        printTableRow(row, widths)
    end

    printTableBorder(widths)

    return nothing
end

function printTableBorder(widths::Vector{Int})
    print("  ")

    for width in widths
        print("|")
        print(repeat("-", width + 2))
    end

    println("|")

    return nothing
end

function printTableRow(row::Vector{String}, widths::Vector{Int})
    print("  ")

    for column in eachindex(row)
        print("| ")
        print(rpad(row[column], widths[column]))
        print(" ")
    end

    println("|")

    return nothing
end

"""
    printEdges(graph::AbstractFactorGraph)

Print one line per graph edge.

# Arguments

- `graph`: Gaussian or discrete factor graph, or tree view.

# Notes

Each line has the form `edge_id: variable <-> factor`.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])

printEdges(graph)
```
"""
function printEdges(graph::AbstractFactorGraph)
    for edge in graph.edges
        factorData = graphFactors(graph)[edge.factorIndex]
        variable = graphVariables(graph)[edge.variableIndex]

        println(
            "  ",
            edge.id,
            ": ",
            variable.label,
            " <-> ",
            factorData.label
        )
    end

    return nothing
end

function printEdges(tree::TreeFactorGraph{GaussianFactorGraph})
    return printEdges(tree.graph)
end

function printEdges(tree::TreeFactorGraph{DiscreteFactorGraph})
    return printEdges(tree.graph)
end

"""
    printMessages(
        graph::AbstractFactorGraph, inference::AbstractInference;
        variable = nothing, factor = nothing
    )

Print messages using exactly one keyword filter:

# Arguments

- `graph`: Gaussian or discrete factor graph, or tree view.
- `inference`: Matching inference object.

# Keywords

- `variable`: Variable ID or label for variable-to-factor messages.
- `factor`: Factor index or label for factor-to-variable messages.

# Notes

Use exactly one keyword filter. `variable = label` prints messages from the variable node
to connected factor nodes. `factor = label` prints messages from the factor node to
connected variable nodes.

The printed fields depend on the inference form: moment form prints mean and
covariance, canonical form prints information and precision, and min-sum form
prints costs or quadratic-message data.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)

messages!(graph, inference)
printMessages(graph, inference; variable = :x1)
```
"""
function printMessages(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference;
    variable = nothing,
    factor = nothing
)
    return printMessages(graph, inference, "moment"; variable = variable, factor = factor)
end

function printMessages(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMomentInference;
    kwargs...
)
    return printMessages(tree.graph, inference; kwargs...)
end

function printMessages(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
    variable = nothing,
    factor = nothing
)
    return printMessages(graph, inference, "canonical"; variable = variable, factor = factor)
end

function printMessages(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianCanonicalInference;
    kwargs...
)
    return printMessages(tree.graph, inference; kwargs...)
end

function printMessages(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    variable = nothing,
    factor = nothing
)
    return printMessages(graph, inference, "min-sum"; variable = variable, factor = factor)
end

function printMessages(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMinSumInference;
    kwargs...
)
    return printMessages(tree.graph, inference; kwargs...)
end

function printMessages(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
    variable = nothing,
    factor = nothing
)
    return printMessages(graph, inference, "sum-product"; variable = variable, factor = factor)
end

function printMessages(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference;
    kwargs...
)
    return printMessages(tree.graph, inference; kwargs...)
end

function printMessages(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    variable = nothing,
    factor = nothing
)
    return printMessages(graph, inference, "min-sum"; variable = variable, factor = factor)
end

function printMessages(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference;
    kwargs...
)
    return printMessages(tree.graph, inference; kwargs...)
end

function printMessages(
    graph::Union{GaussianFactorGraph, DiscreteFactorGraph},
    inference,
    formName::String;
    variable = nothing,
    factor = nothing
)
    assertInferenceMatchesGraph(graph, inference)

    if (variable === nothing && factor === nothing) ||
            (variable !== nothing && factor !== nothing)
        error("Call printMessages with exactly one keyword: variable = label or factor = label.")
    end

    if variable !== nothing
        variableIdx = variableIndex(graphReferenceIndex(graph), variable)
        variableData = graphVariables(graph)[variableIdx]

        println(
            "Messages from variable node \"",
            variableData.label,
            "\" to factor nodes (",
            formName,
            " form):"
        )
        println()

        for edgeId in graph.variableEdges[variableIdx]
            edge = graph.edges[edgeId]
            factorData = graphFactors(graph)[edge.factorIndex]
            message = inference.variableToFactor[edgeId]

            println("Message to factor node \"", factorData.label, "\":")
            printMessageData(message, formName)
            println()
        end

        return nothing
    end

    factorIdx = factorIndex(graph, factor)
    factorData = graphFactors(graph)[factorIdx]

    println(
        "Messages from factor node \"",
        factorData.label,
        "\" to variable nodes (",
        formName,
        " form):"
    )
    println()

    for edgeId in graph.factorEdges[factorIdx]
        edge = graph.edges[edgeId]
        variableData = graphVariables(graph)[edge.variableIndex]
        message = inference.factorToVariable[edgeId]

        println("Message to variable node \"", variableData.label, "\":")
        printMessageData(message, formName)
        println()
    end

    return nothing
end

"""
    printMarginal(
        graph::AbstractFactorGraph, inference::AbstractSumProductInference;
        variable = nothing
    )

Print variable marginals. If `variable` is omitted, all variable marginals are printed.

# Arguments

- `graph`: Gaussian or discrete factor graph, or tree view.
- `inference`: Matching inference object.

# Keywords

- `variable`: Optional variable ID or label.

# Notes

Moment-form marginals print mean and covariance. Canonical-form marginals also
print information and precision.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = moment(graph)
gbp!(graph, inference)

printMarginal(graph, inference; variable = :x1)
```
"""
function printMarginal(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference;
    variable = nothing
)
    return printMarginal(graph, inference, "moment"; variable = variable)
end

function printMarginal(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMomentInference;
    kwargs...
)
    return printMarginal(tree.graph, inference; kwargs...)
end

function printMarginal(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
    variable = nothing
)
    return printMarginal(graph, inference, "canonical"; variable = variable)
end

function printMarginal(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianCanonicalInference;
    kwargs...
)
    return printMarginal(tree.graph, inference; kwargs...)
end

function printMarginal(
    graph::DiscreteFactorGraph,
    inference::DiscreteSumProductInference;
    variable = nothing
)
    return printMarginal(graph, inference, "sum-product"; variable = variable)
end

function printMarginal(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteSumProductInference;
    kwargs...
)
    return printMarginal(tree.graph, inference; kwargs...)
end

function printMarginal(
    graph::Union{GaussianFactorGraph, DiscreteFactorGraph},
    inference,
    formName::String;
    variable = nothing
)
    assertInferenceMatchesGraph(graph, inference)

    if variable === nothing
        println("Variable marginals (", formName, " form):")
        println()

        for variableIndexValue in eachindex(graphVariables(graph))
            printVariableMarginal(graph, inference, variableIndexValue)
        end

        return nothing
    end

    variableIndexValue = variableIndex(graphReferenceIndex(graph), variable)
    variableData = graphVariables(graph)[variableIndexValue]

    println(
        "Marginal for variable node \"",
        variableData.label,
        "\" (",
        formName,
        " form):"
    )
    printMarginalData(inference.marginal[variableIndexValue])

    return nothing
end

function printVariableMarginal(
    graph::Union{GaussianFactorGraph, DiscreteFactorGraph},
    inference,
    variableIndexValue::Int
)
    variableData = graphVariables(graph)[variableIndexValue]

    println("Marginal for variable node \"", variableData.label, "\":")
    printMarginalData(inference.marginal[variableIndexValue])
    println()

    return nothing
end

"""
    printEstimate(
        graph::AbstractFactorGraph, inference::AbstractMinSumInference;
        variable = nothing
    )

Print min-sum MAP estimates. If `variable` is omitted, all variable estimates are printed.

# Arguments

- `graph`: Gaussian or discrete factor graph, or tree view.
- `inference`: Matching min-sum inference object.

# Keywords

- `variable`: Optional variable ID or label.

# Example

```julia
x1 = GaussianVariable(:x1, 1)
f1 = GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "f1")

graph = factorGraph([x1], [f1])
inference = minsum(graph)
gbp!(graph, inference)

printEstimate(graph, inference; variable = :x1)
```
"""
function printEstimate(
    graph::GaussianFactorGraph,
    inference::GaussianMinSumInference;
    variable = nothing
)
    assertGaussianMinSumInferenceMatchesGraph(graph, inference)

    return printEstimateImpl(graph, inference; variable = variable)
end

function printEstimate(
    tree::TreeFactorGraph{GaussianFactorGraph},
    inference::GaussianMinSumInference;
    kwargs...
)
    return printEstimate(tree.graph, inference; kwargs...)
end

function printEstimate(
    graph::DiscreteFactorGraph,
    inference::DiscreteMinSumInference;
    variable = nothing
)
    assertDiscreteMinSumInferenceMatchesGraph(graph, inference)

    return printEstimateImpl(graph, inference; variable = variable)
end

function printEstimate(
    tree::TreeFactorGraph{DiscreteFactorGraph},
    inference::DiscreteMinSumInference;
    kwargs...
)
    return printEstimate(tree.graph, inference; kwargs...)
end

function printEstimateImpl(
    graph::AbstractFactorGraph,
    inference::AbstractMinSumInference;
    variable = nothing
)
    if variable === nothing
        println("Variable MAP estimates (min-sum):")
        println()

        for variableIndexValue in eachindex(graphVariables(graph))
            printVariableEstimate(graph, inference, variableIndexValue)
        end

        return nothing
    end

    variableIndexValue = variableIndex(graphReferenceIndex(graph), variable)
    variableData = graphVariables(graph)[variableIndexValue]

    println("MAP estimate for variable node \"", variableData.label, "\" (min-sum):")
    printEstimateData(inference.estimate[variableIndexValue])

    return nothing
end

function printVariableEstimate(
    graph::AbstractFactorGraph,
    inference::AbstractMinSumInference,
    variableIndexValue::Int
)
    variableData = graphVariables(graph)[variableIndexValue]

    println("MAP estimate for variable node \"", variableData.label, "\":")
    printEstimateData(inference.estimate[variableIndexValue])
    println()

    return nothing
end

function printMessageData(message::GaussianMessage)
    println("  mean = ", message.mean)
    printMatrix("covariance", message.covariance)

    return nothing
end

function printMessageData(message::CanonicalMessage)
    println("  information = ", message.information)
    printMatrix("precision", message.precision)

    return nothing
end

function printMessageData(message::QuadraticMessage)
    printMatrix("J", message.J)
    println("  h = ", message.h)
    println("  c = ", message.c)

    return nothing
end

function printMessageData(message::Vector{Float64})
    println("  probability = ", message)

    return nothing
end

function printMessageData(message, formName::String)
    return printMessageData(message)
end

function printMessageData(message::Vector{Float64}, formName::String)
    label = formName == "min-sum" ? "cost" : "probability"

    println("  ", label, " = ", message)

    return nothing
end

function printMarginalData(marginal::GaussianMessage)
    println("  mean = ", marginal.mean)
    printMatrix("covariance", marginal.covariance)

    return nothing
end

function printMarginalData(marginal::CanonicalMarginal)
    println("  mean = ", marginal.mean)
    printMatrix("covariance", marginal.covariance)
    println("  information = ", marginal.information)
    printMatrix("precision", marginal.precision)

    return nothing
end

function printMarginalData(marginal::Vector{Float64})
    println("  probability = ", marginal)

    return nothing
end

function printEstimateData(estimate::Vector{Float64})
    println("  estimate = ", estimate)

    return nothing
end

function printEstimateData(estimate::StateRef)
    println("  estimate = ", estimate)

    return nothing
end

function printMatrix(label::String, matrix::AbstractMatrix)
    prefix = "  " * label * " = "
    padding = repeat(" ", length(prefix))

    for (rowNumber, row) in enumerate(axes(matrix, 1))
        print(rowNumber == 1 ? prefix : padding)
        printMatrixRow(matrix, row)
    end

    return nothing
end

function printMatrixRow(matrix::AbstractMatrix, row)
    print("[")

    for column in axes(matrix, 2)
        print(matrix[row, column])

        if column != last(axes(matrix, 2))
            print("  ")
        end
    end

    println("]")

    return nothing
end
