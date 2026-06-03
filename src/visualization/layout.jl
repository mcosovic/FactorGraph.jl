function nodeRadii(
    variables::AbstractVector{<:Union{GaussianVariable, DiscreteVariable}},
    defaultRadius::Int,
    labelPlacement::Symbol,
    showLabels::Bool,
    fontSize::Int
)
    if labelPlacement == :outside || !showLabels
        return fill(max(12, min(defaultRadius, 16)), length(variables))
    end

    return [
        max(defaultRadius, ceil(Int, approximateTextWidth(variable.label, fontSize) / 2 + 12))
        for variable in variables
    ]
end

function factorNodeWidths(
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    defaultSize::Int,
    labelPlacement::Symbol,
    showLabels::Bool,
    fontSize::Int
)
    if labelPlacement == :outside || !showLabels
        return fill(max(12, min(defaultSize, 18)), length(factors))
    end

    return [
        max(defaultSize, ceil(Int, approximateTextWidth(factor.label, fontSize) + 18))
        for factor in factors
    ]
end

function factorNodeHeights(
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    defaultSize::Int,
    labelPlacement::Symbol,
    showLabels::Bool
)
    if labelPlacement == :outside || !showLabels
        return fill(max(12, min(defaultSize, 18)), length(factors))
    end

    return fill(max(defaultSize, 34), length(factors))
end

function verticalFactorNodeWidths(
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    defaultSize::Int,
    labelPlacement::Symbol,
    showLabels::Bool,
    fontSize::Int,
    rotateLabels::AbstractVector{Bool}
)
    if labelPlacement == :outside || !showLabels
        return fill(max(12, min(defaultSize, 18)), length(factors))
    end

    return [
        rotateLabels[index] ?
            max(defaultSize, 34) :
            max(defaultSize, ceil(Int, approximateTextWidth(factor.label, fontSize) + 18))
        for (index, factor) in pairs(factors)
    ]
end

function verticalFactorNodeHeights(
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    defaultSize::Int,
    labelPlacement::Symbol,
    showLabels::Bool,
    fontSize::Int,
    rotateLabels::AbstractVector{Bool}
)
    if labelPlacement == :outside || !showLabels
        return fill(max(12, min(defaultSize, 18)), length(factors))
    end

    return [
        rotateLabels[index] ?
            max(defaultSize, ceil(Int, approximateTextWidth(factor.label, fontSize) + 18)) :
            max(defaultSize, 34)
        for (index, factor) in pairs(factors)
    ]
end

function approximateTextWidth(label::AbstractString, fontSize::Int)
    return 0.62 * fontSize * length(label)
end

function labelReserve(
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    factorIndices::AbstractVector{Int},
    fontSize::Int
)
    if isempty(factorIndices)
        return 0.0
    end

    return maximum(
        approximateTextWidth(factors[index].label, fontSize) for index in factorIndices
    ) + 16
end

function minimumNodeGap(
    nodeHeights::AbstractVector,
    nodeIndices::AbstractVector{Int},
    nodeGap::Int
)
    if isempty(nodeIndices)
        return 0
    end

    return maximum(nodeHeights[index] for index in nodeIndices) + nodeGap
end

function verticalFactorLabelRotations(
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    rawFactorX::AbstractVector,
    unaryFactors::AbstractVector{Int},
    multiFactors::AbstractVector{Int},
    defaultSize::Int,
    labelPlacement::Symbol,
    showLabels::Bool,
    fontSize::Int,
    nodeGap::Int
)
    rotateLabels = fill(false, length(factors))

    if !showLabels
        return rotateLabels
    end

    if hasOverlappingVerticalLabels(
        factors,
        rawFactorX,
        unaryFactors,
        defaultSize,
        labelPlacement,
        fontSize,
        nodeGap
    )
        rotateLabels[unaryFactors] .= true
    end

    if hasOverlappingVerticalLabels(
        factors,
        rawFactorX,
        multiFactors,
        defaultSize,
        labelPlacement,
        fontSize,
        nodeGap
    )
        rotateLabels[multiFactors] .= true
    end

    return rotateLabels
end

function hasOverlappingVerticalLabels(
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    rawFactorX::AbstractVector,
    factorIndices::AbstractVector{Int},
    defaultSize::Int,
    labelPlacement::Symbol,
    fontSize::Int,
    nodeGap::Int
)
    if length(factorIndices) <= 1
        return false
    end

    sortedIndices = sort(collect(factorIndices); by = index -> rawFactorX[index])
    previousIndex = first(sortedIndices)

    for index in sortedIndices[2:end]
        previousWidth = horizontalVerticalFactorLabelWidth(
            factors[previousIndex],
            defaultSize,
            labelPlacement,
            fontSize
        )
        currentWidth = horizontalVerticalFactorLabelWidth(
            factors[index],
            defaultSize,
            labelPlacement,
            fontSize
        )
        minimumDistance = previousWidth / 2 + currentWidth / 2 + nodeGap

        if rawFactorX[index] - rawFactorX[previousIndex] < minimumDistance
            return true
        end

        previousIndex = index
    end

    return false
end

function horizontalVerticalFactorLabelWidth(
    factor::Union{GaussianFactor, DiscreteFactor},
    defaultSize::Int,
    labelPlacement::Symbol,
    fontSize::Int
)
    nodeWidth = labelPlacement == :outside ? max(12, min(defaultSize, 18)) : max(defaultSize, 34)
    labelWidth = approximateTextWidth(factor.label, fontSize)

    return labelPlacement == :outside ? max(nodeWidth, labelWidth) : max(nodeWidth, labelWidth + 18)
end

function centeredColumns(
    width::Real,
    margin::Int,
    leftLabelReserve::Real,
    rightLabelReserve::Real,
    maxFactorWidth::Real,
    columnSpacing
)
    unaryFactorSpacing = graphFigureSpacingAt(columnSpacing, 1)
    multiFactorSpacing = graphFigureSpacingAt(columnSpacing, 2)
    leftNodeHalf = maxFactorWidth / 2
    rightNodeHalf = maxFactorWidth / 2
    graphWidth = leftLabelReserve + leftNodeHalf + unaryFactorSpacing +
        multiFactorSpacing + rightNodeHalf + rightLabelReserve
    graphLeft = (width - graphWidth) / 2

    unaryFactorX = graphLeft + leftLabelReserve + leftNodeHalf
    variableX = unaryFactorX + unaryFactorSpacing
    multiFactorX = variableX + multiFactorSpacing

    return variableX, unaryFactorX, multiFactorX
end

function centeredRows(
    height::Real,
    margin::Int,
    topLabelReserve::Real,
    bottomLabelReserve::Real,
    maxFactorHeight::Real,
    rowSpacing
)
    unaryFactorSpacing = graphFigureSpacingAt(rowSpacing, 1)
    multiFactorSpacing = graphFigureSpacingAt(rowSpacing, 2)
    topNodeHalf = maxFactorHeight / 2
    bottomNodeHalf = maxFactorHeight / 2
    graphHeight = topLabelReserve + topNodeHalf + unaryFactorSpacing +
        multiFactorSpacing + bottomNodeHalf + bottomLabelReserve
    graphTop = (height - graphHeight) / 2

    unaryFactorY = graphTop + topLabelReserve + topNodeHalf
    variableY = unaryFactorY + unaryFactorSpacing
    multiFactorY = variableY + multiFactorSpacing

    return variableY, unaryFactorY, multiFactorY
end

function spreadCloseRows(values::AbstractVector, minimumGap::Real)
    if isempty(values)
        return Float64[]
    end

    indexed = collect(enumerate(Float64.(values)))
    sort!(indexed; by = pair -> pair[2])

    spread = similar(Float64.(values))
    previous = -Inf

    for (index, value) in indexed
        y = max(value, previous + minimumGap)
        spread[index] = y
        previous = y
    end

    centerShift = average(spread) - average(Float64.(values))

    return spread .- centerShift
end

function average(values)
    total = 0.0
    count = 0

    for value in values
        total += value
        count += 1
    end

    if count == 0
        error("Cannot compute average of empty values.")
    end

    return total / count
end

function graphFigureSpacingAt(spacing::Real, index::Int)
    return spacing
end

function graphFigureSpacingAt(spacing::Union{Tuple, AbstractVector}, index::Int)
    return spacing[min(index, length(spacing))]
end

function graphFigureSpacingSum(spacing, gapCount::Int)
    total = 0.0

    for index in 1:gapCount
        total += graphFigureSpacingAt(spacing, index)
    end

    return total
end

function graphFigureSpacingMaximum(spacing::Real)
    return spacing
end

function graphFigureSpacingMaximum(spacing::Union{Tuple, AbstractVector})
    return maximum(spacing)
end

function nodeRows(count::Int, rowCount::Int, margin::Int, nodeSpacing)
    if count == 0
        return Float64[]
    end

    startOffset =
        (graphFigureSpacingSum(nodeSpacing, rowCount - 1) -
            graphFigureSpacingSum(nodeSpacing, count - 1)) / 2
    rows = Float64[]
    position = margin + startOffset

    for index in 1:count
        push!(rows, position)
        position += graphFigureSpacingAt(nodeSpacing, index)
    end

    return rows
end

function graphFigurePositions(count::Int, canvasSize::Real, spacing)
    if count == 0
        return Float64[]
    end

    totalSpacing = graphFigureSpacingSum(spacing, count - 1)
    positions = Float64[]
    position = (canvasSize - totalSpacing) / 2

    for index in 1:count
        push!(positions, position)
        position += graphFigureSpacingAt(spacing, index)
    end

    return positions
end

function graphFigureCenteredRows(count::Int, rowCount::Int, canvasSize::Real, spacing)
    if count == 0
        return Float64[]
    end

    totalRows = graphFigureSpacingSum(spacing, count - 1)
    totalSlots = graphFigureSpacingSum(spacing, rowCount - 1)
    position = (canvasSize - totalSlots) / 2 + (totalSlots - totalRows) / 2
    rows = Float64[]
    for index in 1:count
        push!(rows, position)
        position += graphFigureSpacingAt(spacing, index)
    end
    return rows
end

function fittedCanvasSizeAndShift(
    extentMin::Real,
    extentMax::Real,
    canvasSize::Int,
    canvasPadding::Int
)
    extent = extentMax - extentMin
    fittedSize = ceil(Int, extent + 2 * canvasPadding)
    finalSize = max(canvasSize, fittedSize)
    shift = (finalSize - extent) / 2 - extentMin

    return finalSize, shift
end

function fitVerticalCanvas!(
    variableY::AbstractVector,
    factorY::Union{AbstractVector, AbstractDict{Int}},
    variableRadii::AbstractVector,
    factorHeights::AbstractVector,
    canvasHeight::Int,
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    outsideLabelGap::Int,
    fontSize::Int,
    canvasPadding::Int
    ;
    variableIndices = eachindex(variableY),
    factorIndices = eachindex(factorHeights)
)
    yMin = Inf
    yMax = -Inf
    variableBottomPadding = labelPlacement == :outside && showVariableLabels ?
        outsideLabelGap + fontSize : 0
    factorLabelPadding = labelPlacement == :outside && showFactorLabels ?
        outsideLabelGap + fontSize : 0

    for index in variableIndices
        yMin = min(yMin, variableY[index] - variableRadii[index])
        yMax = max(yMax, variableY[index] + variableRadii[index] + variableBottomPadding)
    end

    for index in factorIndices
        y = factorY[index]
        halfHeight = factorHeights[index] / 2
        yMin = min(yMin, y - halfHeight - factorLabelPadding)
        yMax = max(yMax, y + halfHeight + factorLabelPadding)
    end

    if !isfinite(yMin)
        return canvasHeight
    end

    finalHeight, shift = fittedCanvasSizeAndShift(yMin, yMax, canvasHeight, canvasPadding)

    if !iszero(shift)
        for index in variableIndices
            variableY[index] += shift
        end

        for index in factorIndices
            factorY[index] += shift
        end
    end

    return finalHeight
end

function fitGeneralHorizontalCanvas!(
    variableX::Real,
    factorX::AbstractDict{Int},
    variableRadii::AbstractVector,
    factorWidths::AbstractVector,
    canvasWidth::Int,
    variables::AbstractVector{<:Union{GaussianVariable, DiscreteVariable}},
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    outsideLabelGap::Int,
    fontSize::Int,
    canvasPadding::Int
    ;
    variableIndices = eachindex(variableRadii),
    factorIndices = eachindex(factorWidths)
)
    xMin = minimum(variableX - variableRadii[index] for index in variableIndices)
    xMax = maximum(variableX + variableRadii[index] for index in variableIndices)

    if labelPlacement == :outside && showVariableLabels
        variableLabelWidth = maximum(
            approximateTextWidth(variables[index].label, fontSize) for index in variableIndices
        )
        xMin = min(xMin, variableX - variableLabelWidth / 2)
        xMax = max(xMax, variableX + variableLabelWidth / 2)
    end

    for index in factorIndices
        labelWidth = labelPlacement == :outside && showFactorLabels ?
            outsideLabelGap + approximateTextWidth(factors[index].label, fontSize) : 0.0
        xMin = min(xMin, factorX[index] - factorWidths[index] / 2 - labelWidth)
        xMax = max(xMax, factorX[index] + factorWidths[index] / 2 + labelWidth)
    end

    finalWidth, shift = fittedCanvasSizeAndShift(xMin, xMax, canvasWidth, canvasPadding)

    if !iszero(shift)
        variableX += shift

        for index in factorIndices
            factorX[index] += shift
        end
    end

    return finalWidth, variableX
end

function fitVerticalGraphHorizontalCanvas!(
    variableX::AbstractVector,
    factorX::AbstractDict{Int},
    variableRadii::AbstractVector,
    factorWidths::AbstractVector,
    canvasWidth::Int,
    variables::AbstractVector{<:Union{GaussianVariable, DiscreteVariable}},
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    rotateFactorLabels::AbstractVector{Bool},
    outsideLabelGap::Int,
    fontSize::Int,
    canvasPadding::Int
    ;
    variableIndices = eachindex(variableX),
    factorIndices = eachindex(factorWidths)
)
    xMin = Inf
    xMax = -Inf

    for index in variableIndices
        xMin = min(xMin, variableX[index] - variableRadii[index])
        xMax = max(xMax, variableX[index] + variableRadii[index])

        if labelPlacement == :outside && showVariableLabels
            labelPadding = outsideLabelGap + approximateTextWidth(variables[index].label, fontSize)
            xMax = max(xMax, variableX[index] + variableRadii[index] + labelPadding)
        end
    end

    for index in factorIndices
        xMin = min(xMin, factorX[index] - factorWidths[index] / 2)
        xMax = max(xMax, factorX[index] + factorWidths[index] / 2)

        if labelPlacement == :outside && showFactorLabels && !rotateFactorLabels[index]
            labelHalfWidth = approximateTextWidth(factors[index].label, fontSize) / 2
            xMin = min(xMin, factorX[index] - labelHalfWidth)
            xMax = max(xMax, factorX[index] + labelHalfWidth)
        end
    end

    if !isfinite(xMin)
        return canvasWidth
    end

    finalWidth, shift = fittedCanvasSizeAndShift(xMin, xMax, canvasWidth, canvasPadding)

    if !iszero(shift)
        for index in variableIndices
            variableX[index] += shift
        end

        for index in factorIndices
            factorX[index] += shift
        end
    end

    return finalWidth
end

function fitVerticalGraphVerticalCanvas!(
    variableY::Real,
    factorY::AbstractDict{Int},
    variableRadii::AbstractVector,
    factorHeights::AbstractVector,
    canvasHeight::Int,
    variables::AbstractVector{<:Union{GaussianVariable, DiscreteVariable}},
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    rotateFactorLabels::AbstractVector{Bool},
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    outsideLabelGap::Int,
    fontSize::Int,
    canvasPadding::Int
    ;
    variableIndices = eachindex(variableRadii),
    factorIndices = eachindex(factorHeights)
)
    yMin = minimum(variableY - variableRadii[index] for index in variableIndices)
    yMax = maximum(variableY + variableRadii[index] for index in variableIndices)

    if labelPlacement == :outside && showVariableLabels
        labelY = variableY + maximum(
            variableRadii[index] for index in variableIndices; init = 0
        ) + outsideLabelGap + fontSize / 2
        yMax = max(yMax, labelY + fontSize / 2)
    end

    for index in factorIndices
        y = factorY[index]
        halfHeight = factorHeights[index] / 2
        yMin = min(yMin, y - halfHeight)
        yMax = max(yMax, y + halfHeight)

        if labelPlacement == :outside && showFactorLabels
            if rotateFactorLabels[index]
                labelHeight = approximateTextWidth(factors[index].label, fontSize)
                labelY = y < variableY ?
                    y - halfHeight - outsideLabelGap :
                    y + halfHeight + outsideLabelGap
                yMin = min(yMin, y < variableY ? labelY - labelHeight : labelY)
                yMax = max(yMax, y < variableY ? labelY : labelY + labelHeight)
            else
                labelHalfHeight = fontSize / 2
                labelY = y < variableY ?
                    y - halfHeight - outsideLabelGap - labelHalfHeight :
                    y + halfHeight + outsideLabelGap + labelHalfHeight
                yMin = min(yMin, labelY - labelHalfHeight)
                yMax = max(yMax, labelY + labelHalfHeight)
            end
        end
    end

    finalHeight, shift = fittedCanvasSizeAndShift(yMin, yMax, canvasHeight, canvasPadding)

    if !iszero(shift)
        variableY += shift

        for index in factorIndices
            factorY[index] += shift
        end
    end

    return finalHeight, variableY
end

function fitHorizontalCanvas(
    variableX::AbstractVector,
    factorX::AbstractVector,
    variableRadii::AbstractVector,
    factorWidths::AbstractVector,
    canvasWidth::Int,
    variables::AbstractVector{<:Union{GaussianVariable, DiscreteVariable}},
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    outsideLabelGap::Int,
    fontSize::Int,
    canvasPadding::Int
    ;
    variableIndices = eachindex(variableX),
    factorIndices = eachindex(factorX)
)
    xMin = Inf
    xMax = -Inf

    for index in variableIndices
        xMin = min(xMin, variableX[index] - variableRadii[index])
        xMax = max(xMax, variableX[index] + variableRadii[index])

        if labelPlacement == :outside && showVariableLabels
            labelHalfWidth = approximateTextWidth(variables[index].label, fontSize) / 2
            xMin = min(xMin, variableX[index] - labelHalfWidth)
            xMax = max(xMax, variableX[index] + labelHalfWidth)
        end
    end

    for index in factorIndices
        labelPadding = labelPlacement == :outside && showFactorLabels ?
            outsideLabelGap + approximateTextWidth(factors[index].label, fontSize) : 0.0
        xMin = min(xMin, factorX[index] - factorWidths[index] / 2)
        xMax = max(xMax, factorX[index] + factorWidths[index] / 2 + labelPadding)
    end

    if !isfinite(xMin)
        return canvasWidth
    end

    finalWidth, shift = fittedCanvasSizeAndShift(xMin, xMax, canvasWidth, canvasPadding)

    if !iszero(shift)
        for index in variableIndices
            variableX[index] += shift
        end

        for index in factorIndices
            factorX[index] += shift
        end
    end

    return finalWidth
end

function fitSideLabelHorizontalCanvas!(
    variableX::AbstractVector,
    factorX::AbstractVector,
    variableRadii::AbstractVector,
    factorWidths::AbstractVector,
    canvasWidth::Int,
    variables::AbstractVector{<:Union{GaussianVariable, DiscreteVariable}},
    factors::AbstractVector{<:Union{GaussianFactor, DiscreteFactor}},
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    outsideLabelGap::Int,
    fontSize::Int,
    canvasPadding::Int
    ;
    variableIndices = eachindex(variableX),
    factorIndices = eachindex(factorX)
)
    xMin = Inf
    xMax = -Inf

    for index in variableIndices
        xMin = min(xMin, variableX[index] - variableRadii[index])
        xMax = max(xMax, variableX[index] + variableRadii[index])

        if labelPlacement == :outside && showVariableLabels
            labelPadding = outsideLabelGap + approximateTextWidth(variables[index].label, fontSize)
            xMax = max(xMax, variableX[index] + variableRadii[index] + labelPadding)
        end
    end

    for index in factorIndices
        xMin = min(xMin, factorX[index] - factorWidths[index] / 2)
        xMax = max(xMax, factorX[index] + factorWidths[index] / 2)

        if labelPlacement == :outside && showFactorLabels
            labelPadding = outsideLabelGap + approximateTextWidth(factors[index].label, fontSize)
            xMax = max(xMax, factorX[index] + factorWidths[index] / 2 + labelPadding)
        end
    end

    if !isfinite(xMin)
        return canvasWidth
    end

    finalWidth, shift = fittedCanvasSizeAndShift(xMin, xMax, canvasWidth, canvasPadding)

    if !iszero(shift)
        for index in variableIndices
            variableX[index] += shift
        end

        for index in factorIndices
            factorX[index] += shift
        end
    end

    return finalWidth
end
