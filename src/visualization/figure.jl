"""
    graphFigure(graph::AbstractFactorGraph; kwargs...)

Return a SVG visualization of a factor graph.

# Arguments

- `graph`: Factor graph or tree view to draw.

# Keywords

Canvas options control the size, padding, and displayed scale of the SVG:
- `canvas = NamedTuple()`: Override canvas options. Supported keys are:
  - `width`: Minimum SVG canvas width. If `nothing`, fit the graph tightly with
    `padding`.
  - `height`: Minimum SVG canvas height. If `nothing`, fit the graph tightly
    with `padding`.
  - `padding`: Minimum padding around the drawn graph.
  - `zoom`: Scale the displayed SVG size while preserving the same view box.

Layout options control graph orientation, spacing, and edge geometry:
- `layout = NamedTuple()`: Override layout options. Supported keys are:
  - `orientation`: Use `:horizontal` or `:vertical` layout.
  - `rowSpacing`: Vertical spacing between node rows. May be a scalar or a
    tuple/vector of per-gap spacings; the last per-gap value is reused.
  - `columnSpacing`: Horizontal spacing between node columns. May be a scalar
    or a tuple/vector of per-gap spacings; the last per-gap value is reused.
  - `curvedEdges`: Draw edges as cubic curves instead of straight lines.

Node options control the size of variable and factor nodes:
- `node = NamedTuple()`: Override node options. Supported keys are:
  - `variableRadius`: Default variable node radius.
  - `factorSize`: Default factor node size.

Label options control label placement, visibility, and font size:
- `label = NamedTuple()`: Override label options. Supported keys are:
  - `placement`: Place labels `:outside` nodes or `:inside` nodes.
  - `outsideGap`: Gap between a node and its outside label.
  - `showVariables`: Draw variable labels.
  - `showFactors`: Draw factor labels.
  - `showTooltips`: Add SVG hover tooltips for variables, factors, and edges.
  - `showEdgeIds`: Draw edge id labels next to edges.
  - `tooltipDetail`: Tooltip detail level: `:summary` or `:full`.
  - `fontSize`: Label font size in SVG user units.

View options control which part of a factor graph or tree view is drawn:
- `view = NamedTuple()`: Override view options. Supported keys are:
  - `variables`: Variable ids or labels used as the focus.
  - `factors`: Factor ids or labels used as the focus.
  - `hops`: Nonnegative number of bipartite graph hops from focused variables
    and factors, or `:all` to expand through the connected component.

Style options control default graph colors and stroke widths:
- `style = NamedTuple()`: Override default SVG styles. Supported keys are:
  - `backgroundFill`: SVG background fill color.
  - `variableFill`: Variable node fill color.
  - `variableStroke`: Variable node stroke color.
  - `variableStrokeWidth`: Variable node stroke width.
  - `factorFill`: Factor node fill color.
  - `factorStroke`: Factor node stroke color.
  - `factorStrokeWidth`: Factor node stroke width.
  - `edgeStroke`: Edge stroke color.
  - `edgeStrokeWidth`: Edge stroke width.
  - `edgeOpacity`: Edge stroke opacity.
  - `labelFill`: Label text fill color.

Highlight options control highlighted variables, factors, and edges:
- `highlight = NamedTuple[]`: Highlight entries. Supported entry selectors are:
  - `(variable = x1, ...)`: Highlight a variable node.
  - `(factor = f1, ...)`: Highlight a factor node.
  - `(edge = 1, ...)`: Highlight an edge by edge id or edge object.
  - `(variable = x1, factor = f1, ...)`: Highlight the edge connecting a
    variable and factor.

  Supported style keys are:
  - `stroke`: Node or edge stroke color.
  - `fill`: Node fill color. Applies only to variable and factor entries.
  - `strokeWidth`: Stroke width for the selected variable, factor, or edge.
  - `incidentEdges`: Whether a variable or factor entry highlights incident
    edges. Defaults to `true`.
  - `edgeStroke`: Incident edge stroke color for variable and factor entries.
  - `edgeStrokeWidth`: Incident edge stroke width for variable and factor
    entries. Defaults to `strokeWidth`.

# Returns

An SVG document as a `String`.

# Notes

For a general factor graph, horizontal layout draws unary factors on the left,
variables in the middle, and multi-variable factors on the right. Vertical
layout draws unary factors above variables and multi-variable factors below
variables.

For a `TreeFactorGraph`, the layout uses graph depth. In horizontal tree layout,
`columnSpacing` separates depth levels and `rowSpacing` separates nodes within a
level. In vertical tree layout, `rowSpacing` separates depth levels and
`columnSpacing` separates nodes within a level.

# Example

```julia
variables = [
    GaussianVariable(:x1, 1; label = "x1"),
    GaussianVariable(:x2, 1; label = "x2")
]
factors = [
    GaussianFactor(:x1, 0.0, 1.0, 1.0; label = "prior"),
    GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 1.0; label = "link")
]
graph = factorGraph(variables, factors)

svg = graphFigure(
    graph;
    highlight = [(variable = :x1, stroke = "#16a34a", fill = "#dcfce7", strokeWidth = 4)]
)
```
"""
function graphFigure(
    graph::AbstractFactorGraph;
    canvas::NamedTuple = NamedTuple(),
    layout::NamedTuple = NamedTuple(),
    node::NamedTuple = NamedTuple(),
    label::NamedTuple = NamedTuple(),
    view::NamedTuple = NamedTuple(),
    style::NamedTuple = NamedTuple(),
    highlight::AbstractVector = NamedTuple[],
)
    canvasOptions = graphFigureCanvasOptions(canvas)
    layoutOptions = graphFigureLayoutOptions(layout)
    nodeOptions = graphFigureNodeOptions(node)
    labelOptions = graphFigureLabelOptions(label)
    viewOptions = graphFigureViewOptions(graph, view)
    graphStyle = graphFigureStyle(style)

    variableCount = length(viewOptions.variables)
    unaryFactors = [index for index in viewOptions.factors if isUnaryFactor(graph, index)]
    multiFactors = [index for index in viewOptions.factors if !isUnaryFactor(graph, index)]
    margin = 48
    nodeGap = 10

    if layoutOptions.orientation == :vertical
        return verticalGraphFigure(
            graph,
            variableCount,
            unaryFactors,
            multiFactors;
            width = canvasOptions.width,
            height = canvasOptions.height,
            canvasPadding = canvasOptions.padding,
            zoom = canvasOptions.zoom,
            margin = margin,
            rowSpacing = layoutOptions.rowSpacing,
            columnSpacing = layoutOptions.columnSpacing,
            variableRadius = nodeOptions.variableRadius,
            factorSize = nodeOptions.factorSize,
            labelPlacement = labelOptions.placement,
            outsideLabelGap = labelOptions.outsideGap,
            nodeGap = nodeGap,
            curvedEdges = layoutOptions.curvedEdges,
            style = graphStyle,
            highlight = highlight,
            showVariableLabels = labelOptions.showVariables,
            showFactorLabels = labelOptions.showFactors,
            showTooltips = labelOptions.showTooltips,
            showEdgeIds = labelOptions.showEdgeIds,
            tooltipDetail = labelOptions.tooltipDetail,
            fontSize = labelOptions.fontSize,
            view = viewOptions
        )
    end

    rowCount = max(variableCount, length(unaryFactors), length(multiFactors), 1)
    defaultLayoutHeight =
        2 * margin + graphFigureSpacingSum(layoutOptions.rowSpacing, rowCount - 1)
    layoutHeight = isnothing(canvasOptions.height) ?
        defaultLayoutHeight : max(canvasOptions.height, defaultLayoutHeight)
    minCanvasHeight = isnothing(canvasOptions.height) ? 0 : canvasOptions.height

    layoutNodeSpacing = layoutOptions.rowSpacing

    variableRadii = nodeRadii(
        graph.variables,
        nodeOptions.variableRadius,
        labelOptions.placement,
        labelOptions.showVariables,
        labelOptions.fontSize
    )
    factorWidths = factorNodeWidths(
        graph.factors,
        nodeOptions.factorSize,
        labelOptions.placement,
        labelOptions.showFactors,
        labelOptions.fontSize
    )
    factorHeights = factorNodeHeights(
        graph.factors,
        nodeOptions.factorSize,
        labelOptions.placement,
        labelOptions.showFactors
    )
    leftLabelReserve = labelOptions.placement == :outside && labelOptions.showFactors ?
        labelReserve(graph.factors, unaryFactors, labelOptions.fontSize) : 0.0
    rightLabelReserve = labelOptions.placement == :outside && labelOptions.showFactors ?
        labelReserve(graph.factors, multiFactors, labelOptions.fontSize) : 0.0

    maxFactorWidth = maximum(factorWidths; init = nodeOptions.factorSize)
    defaultLayoutWidth = 2 * margin + leftLabelReserve + maxFactorWidth +
        graphFigureSpacingSum(layoutOptions.columnSpacing, 2) + rightLabelReserve
    layoutWidth = isnothing(canvasOptions.width) ?
        defaultLayoutWidth : max(canvasOptions.width, defaultLayoutWidth)
    minCanvasWidth = isnothing(canvasOptions.width) ? 0 : canvasOptions.width
    variableX, unaryFactorX, multiFactorX = centeredColumns(
        layoutWidth,
        margin,
        leftLabelReserve,
        rightLabelReserve,
        maxFactorWidth,
        layoutOptions.columnSpacing
    )
    variableY = fill(0.0, length(graph.variables))
    orderedVariableY =
        graphFigureCenteredRows(variableCount, rowCount, layoutHeight, layoutNodeSpacing)
    for (position, variableIndex) in pairs(viewOptions.variables)
        variableY[variableIndex] = orderedVariableY[position]
    end
    factorY = Dict{Int, Float64}()
    factorX = Dict{Int, Float64}()

    rawFactorY = fill(0.0, length(graph.factors))
    for factorIndex in viewOptions.factors
        rawFactorY[factorIndex] = average(
            variableY[edge.variableIndex] for edge in viewOptions.edges
            if edge.factorIndex == factorIndex
        )
    end
    unaryFactorY = spreadCloseRows(
        rawFactorY[unaryFactors],
        minimumNodeGap(factorHeights, unaryFactors, nodeGap)
    )
    multiFactorY = spreadCloseRows(
        rawFactorY[multiFactors],
        minimumNodeGap(factorHeights, multiFactors, nodeGap)
    )

    for (position, factorIndex) in pairs(unaryFactors)
        factorY[factorIndex] = unaryFactorY[position]
        factorX[factorIndex] = unaryFactorX
    end

    for (position, factorIndex) in pairs(multiFactors)
        factorY[factorIndex] = multiFactorY[position]
        factorX[factorIndex] = multiFactorX
    end

    canvasWidth, variableX = fitGeneralHorizontalCanvas!(
        variableX,
        factorX,
        variableRadii,
        factorWidths,
        minCanvasWidth,
        graph.variables,
        graph.factors,
        labelOptions.placement,
        labelOptions.showVariables,
        labelOptions.showFactors,
        labelOptions.outsideGap,
        labelOptions.fontSize,
        canvasOptions.padding;
        variableIndices = viewOptions.variables,
        factorIndices = viewOptions.factors
    )
    canvasHeight = fitVerticalCanvas!(
        variableY,
        factorY,
        variableRadii,
        factorHeights,
        minCanvasHeight,
        labelOptions.placement,
        labelOptions.showVariables,
        labelOptions.showFactors,
        labelOptions.outsideGap,
        labelOptions.fontSize,
        canvasOptions.padding;
        variableIndices = viewOptions.variables,
        factorIndices = viewOptions.factors
    )

    highlightState = graphFigureHighlightState(
        graph,
        highlight,
        graphStyle
    )
    buffer = startGraphFigureBuffer(
        canvasWidth,
        canvasHeight,
        canvasOptions.zoom,
        labelOptions.fontSize,
        graphStyle
    )

    occupiedEdgeLabels = NamedTuple[]
    edgeLabelBackgrounds = String[]
    edgeLabelTexts = String[]
    edgeGeometries = NamedTuple[]

    for edge in viewOptions.edges
        currentFactorX = factorX[edge.factorIndex]

        if currentFactorX < variableX
            x1 = currentFactorX + factorWidths[edge.factorIndex] / 2
            y1 = factorY[edge.factorIndex]
            x2 = variableX - variableRadii[edge.variableIndex]
            y2 = variableY[edge.variableIndex]
        else
            x1 = variableX + variableRadii[edge.variableIndex]
            y1 = variableY[edge.variableIndex]
            x2 = currentFactorX - factorWidths[edge.factorIndex] / 2
            y2 = factorY[edge.factorIndex]
        end

        push!(
            edgeGeometries,
            edgeGeometry(edge.id, x1, y1, x2, y2, layoutOptions.curvedEdges, :horizontal)
        )
    end

    for edge in viewOptions.edges
        currentFactorX = factorX[edge.factorIndex]

        if currentFactorX < variableX
            x1 = currentFactorX + factorWidths[edge.factorIndex] / 2
            y1 = factorY[edge.factorIndex]
            x2 = variableX - variableRadii[edge.variableIndex]
            y2 = variableY[edge.variableIndex]
        else
            x1 = variableX + variableRadii[edge.variableIndex]
            y1 = variableY[edge.variableIndex]
            x2 = currentFactorX - factorWidths[edge.factorIndex] / 2
            y2 = factorY[edge.factorIndex]
        end

        edgeClass = edgeClassName(
            edge,
            highlightState.edges
        )
        tooltip = labelOptions.showTooltips ?
            graphFigureEdgeTooltip(graph, edge, labelOptions.tooltipDetail) : nothing
        println(
            buffer,
            svgEdge(
                x1,
                y1,
                x2,
                y2,
                layoutOptions.curvedEdges,
                edgeClass,
                get(highlightState.edgeStyles, edge.id, nothing),
                tooltip
            )
        )
        if labelOptions.showEdgeIds
            labelParts = svgEdgeIdLabelParts(
                edge.id,
                x1,
                y1,
                x2,
                y2,
                occupiedEdgeLabels,
                layoutOptions.curvedEdges,
                :horizontal,
                edgeGeometries,
                tooltip
            )
            push!(edgeLabelBackgrounds, labelParts.background)
            push!(edgeLabelTexts, labelParts.text)
        end
    end

    drawEdgeIdLabels!(buffer, edgeLabelBackgrounds, edgeLabelTexts)

    for index in viewOptions.variables
        variable = graph.variables[index]
        label = escapeXML(variable.label)
        y = variableY[index]
        className = variableClassName(
            index,
            highlightState.variables,
            viewOptions.focusVariableSet
        )
        tooltip = labelOptions.showTooltips ?
            graphFigureVariableTooltip(graph, index, labelOptions.tooltipDetail) : nothing
        println(
            buffer,
            svgCircle(
                variableX,
                y,
                variableRadii[index],
                className,
                get(highlightState.variableStyles, index, nothing),
                tooltip
            )
        )
        if labelOptions.showVariables
            drawVariableLabel!(
                buffer,
                variableX,
                y,
                variableRadii[index],
                label,
                labelOptions.placement,
                labelOptions.outsideGap,
                labelOptions.fontSize,
                tooltip
            )
        end
    end

    for index in viewOptions.factors
        factor = graph.factors[index]
        label = escapeXML(factor.label)
        width = factorWidths[index]
        height = factorHeights[index]
        y = factorY[index]
        x = factorX[index] - width / 2
        className = factorClassName(index, highlightState.factors, viewOptions.focusFactorSet)
        tooltip = labelOptions.showTooltips ?
            graphFigureFactorTooltip(graph, index, labelOptions.tooltipDetail) : nothing
        println(
            buffer,
            svgRect(
                x,
                y - height / 2,
                width,
                height,
                labelOptions.placement,
                className,
                get(highlightState.factorStyles, index, nothing),
                tooltip
            )
        )
        if labelOptions.showFactors
            drawFactorLabel!(
                buffer,
                factorX[index],
                y,
                width,
                variableX,
                label,
                labelOptions.placement,
                labelOptions.outsideGap,
                labelOptions.fontSize,
                tooltip
            )
        end
    end

    println(buffer, "</svg>")

    return String(take!(buffer))
end

function verticalGraphFigure(
    graph::AbstractFactorGraph,
    variableCount::Int,
    unaryFactors::AbstractVector{Int},
    multiFactors::AbstractVector{Int};
    width::Union{Nothing, Int},
    height::Union{Nothing, Int},
    canvasPadding::Int,
    zoom::Real,
    margin::Int,
    rowSpacing,
    columnSpacing,
    variableRadius::Int,
    factorSize::Int,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    nodeGap::Int,
    curvedEdges::Bool,
    style::NamedTuple,
    highlight::AbstractVector,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    showTooltips::Bool,
    showEdgeIds::Bool,
    tooltipDetail::Symbol,
    fontSize::Int,
    view::NamedTuple
)
    columnCount = max(variableCount, length(unaryFactors), length(multiFactors), 1)
    variableRadii = nodeRadii(
        graph.variables,
        variableRadius,
        labelPlacement,
        showVariableLabels,
        fontSize
    )
    columnSpacing = verticalGraphColumnSpacing(
        graph,
        view.variables,
        variableRadii,
        columnSpacing,
        outsideLabelGap,
        labelPlacement,
        showVariableLabels,
        fontSize
    )
    defaultLayoutWidth = 2 * margin + graphFigureSpacingSum(columnSpacing, columnCount - 1)
    layoutWidth = isnothing(width) ? defaultLayoutWidth : max(width, defaultLayoutWidth)
    minCanvasWidth = isnothing(width) ? 0 : width
    layoutColumnSpacing = columnSpacing

    rawVariableX = fill(0.0, length(graph.variables))
    orderedVariableX =
        graphFigureCenteredRows(variableCount, columnCount, layoutWidth, layoutColumnSpacing)
    for (position, variableIndex) in pairs(view.variables)
        rawVariableX[variableIndex] = orderedVariableX[position]
    end
    rawFactorX = fill(0.0, length(graph.factors))
    for factorIndex in view.factors
        rawFactorX[factorIndex] = average(
            rawVariableX[edge.variableIndex] for edge in view.edges
            if edge.factorIndex == factorIndex
        )
    end
    rotateFactorLabels = verticalFactorLabelRotations(
        graph.factors,
        rawFactorX,
        unaryFactors,
        multiFactors,
        factorSize,
        labelPlacement,
        showFactorLabels,
        fontSize,
        nodeGap
    )
    factorWidths = verticalFactorNodeWidths(
        graph.factors,
        factorSize,
        labelPlacement,
        showFactorLabels,
        fontSize,
        rotateFactorLabels
    )
    factorHeights = verticalFactorNodeHeights(
        graph.factors,
        factorSize,
        labelPlacement,
        showFactorLabels,
        fontSize,
        rotateFactorLabels
    )
    topLabelReserve = labelPlacement == :outside && showFactorLabels ?
        labelReserve(graph.factors, unaryFactors, fontSize) : 0.0
    bottomLabelReserve = labelPlacement == :outside && showFactorLabels ?
        labelReserve(graph.factors, multiFactors, fontSize) : 0.0

    maxFactorHeight = maximum(factorHeights; init = factorSize)
    defaultLayoutHeight = 2 * margin + topLabelReserve + maxFactorHeight +
        graphFigureSpacingSum(rowSpacing, 2) + bottomLabelReserve
    layoutHeight = isnothing(height) ? defaultLayoutHeight : max(height, defaultLayoutHeight)
    minCanvasHeight = isnothing(height) ? 0 : height
    variableY, unaryFactorY, multiFactorY = centeredRows(
        layoutHeight,
        margin,
        topLabelReserve,
        bottomLabelReserve,
        maxFactorHeight,
        rowSpacing
    )

    variableX = rawVariableX
    factorX = Dict{Int, Float64}()
    factorY = Dict{Int, Float64}()
    unaryFactorX = spreadCloseRows(
        rawFactorX[unaryFactors],
        minimumNodeGap(factorWidths, unaryFactors, nodeGap)
    )
    multiFactorX = spreadCloseRows(
        rawFactorX[multiFactors],
        minimumNodeGap(factorWidths, multiFactors, nodeGap)
    )

    for (position, factorIndex) in pairs(unaryFactors)
        factorX[factorIndex] = unaryFactorX[position]
        factorY[factorIndex] = unaryFactorY
    end

    for (position, factorIndex) in pairs(multiFactors)
        factorX[factorIndex] = multiFactorX[position]
        factorY[factorIndex] = multiFactorY
    end

    canvasWidth = fitVerticalGraphHorizontalCanvas!(
        variableX,
        factorX,
        variableRadii,
        factorWidths,
        minCanvasWidth,
        graph.variables,
        graph.factors,
        labelPlacement,
        showVariableLabels,
        showFactorLabels,
        rotateFactorLabels,
        outsideLabelGap,
        fontSize,
        canvasPadding;
        variableIndices = view.variables,
        factorIndices = view.factors
    )
    canvasHeight, variableY = fitVerticalGraphVerticalCanvas!(
        variableY,
        factorY,
        variableRadii,
        factorHeights,
        minCanvasHeight,
        graph.variables,
        graph.factors,
        rotateFactorLabels,
        labelPlacement,
        showVariableLabels,
        showFactorLabels,
        outsideLabelGap,
        fontSize,
        canvasPadding;
        variableIndices = view.variables,
        factorIndices = view.factors
    )

    highlightState = graphFigureHighlightState(
        graph,
        highlight,
        style
    )
    buffer = startGraphFigureBuffer(
        canvasWidth,
        canvasHeight,
        zoom,
        fontSize,
        style
    )

    occupiedEdgeLabels = NamedTuple[]
    edgeLabelBackgrounds = String[]
    edgeLabelTexts = String[]
    edgeGeometries = NamedTuple[]

    for edge in view.edges
        if factorY[edge.factorIndex] < variableY
            x1 = factorX[edge.factorIndex]
            y1 = factorY[edge.factorIndex] + factorHeights[edge.factorIndex] / 2
            x2 = variableX[edge.variableIndex]
            y2 = variableY - variableRadii[edge.variableIndex]
        else
            x1 = variableX[edge.variableIndex]
            y1 = variableY + variableRadii[edge.variableIndex]
            x2 = factorX[edge.factorIndex]
            y2 = factorY[edge.factorIndex] - factorHeights[edge.factorIndex] / 2
        end

        push!(
            edgeGeometries,
            edgeGeometry(edge.id, x1, y1, x2, y2, curvedEdges, :vertical)
        )
    end

    for edge in view.edges
        if factorY[edge.factorIndex] < variableY
            x1 = factorX[edge.factorIndex]
            y1 = factorY[edge.factorIndex] + factorHeights[edge.factorIndex] / 2
            x2 = variableX[edge.variableIndex]
            y2 = variableY - variableRadii[edge.variableIndex]
        else
            x1 = variableX[edge.variableIndex]
            y1 = variableY + variableRadii[edge.variableIndex]
            x2 = factorX[edge.factorIndex]
            y2 = factorY[edge.factorIndex] - factorHeights[edge.factorIndex] / 2
        end

        edgeClass = edgeClassName(
            edge,
            highlightState.edges
        )
        tooltip = showTooltips ? graphFigureEdgeTooltip(graph, edge, tooltipDetail) : nothing
        println(
            buffer,
            svgEdge(
                x1,
                y1,
                x2,
                y2,
                curvedEdges,
                edgeClass,
                get(highlightState.edgeStyles, edge.id, nothing),
                tooltip,
                :vertical
            )
        )
        if showEdgeIds
            labelParts = svgEdgeIdLabelParts(
                edge.id,
                x1,
                y1,
                x2,
                y2,
                occupiedEdgeLabels,
                curvedEdges,
                :vertical,
                edgeGeometries,
                tooltip
            )
            push!(edgeLabelBackgrounds, labelParts.background)
            push!(edgeLabelTexts, labelParts.text)
        end
    end

    drawEdgeIdLabels!(buffer, edgeLabelBackgrounds, edgeLabelTexts)

    for index in view.variables
        variable = graph.variables[index]
        label = escapeXML(variable.label)
        className = variableClassName(index, highlightState.variables, view.focusVariableSet)
        tooltip = showTooltips ? graphFigureVariableTooltip(graph, index, tooltipDetail) : nothing
        println(
            buffer,
            svgCircle(
                variableX[index],
                variableY,
                variableRadii[index],
                className,
                get(highlightState.variableStyles, index, nothing),
                tooltip
            )
        )

        if showVariableLabels
            drawSideLabel!(
                buffer,
                variableX[index] + variableRadii[index] + outsideLabelGap,
                variableY,
                label,
                labelPlacement,
                :right,
                tooltip
            )
        end
    end

    for index in view.factors
        factor = graph.factors[index]
        label = escapeXML(factor.label)
        className = factorClassName(index, highlightState.factors, view.focusFactorSet)
        x = factorX[index] - factorWidths[index] / 2
        y = factorY[index] - factorHeights[index] / 2
        println(
            buffer,
            svgRect(
                x,
                y,
                factorWidths[index],
                factorHeights[index],
                labelPlacement,
                className,
                get(highlightState.factorStyles, index, nothing),
                showTooltips ? graphFigureFactorTooltip(graph, index, tooltipDetail) : nothing
            )
        )

        if showFactorLabels
            drawAdaptiveVerticalFactorLabel!(
                buffer,
                factorX[index],
                factorY[index],
                factorHeights[index],
                variableY,
                label,
                rotateFactorLabels[index],
                labelPlacement,
                outsideLabelGap,
                fontSize,
                showTooltips ? graphFigureFactorTooltip(graph, index, tooltipDetail) : nothing
            )
        end
    end

    println(buffer, "</svg>")

    return String(take!(buffer))
end

function graphFigure(
    tree::TreeFactorGraph;
    canvas::NamedTuple = NamedTuple(),
    layout::NamedTuple = NamedTuple(),
    node::NamedTuple = NamedTuple(),
    label::NamedTuple = NamedTuple(),
    view::NamedTuple = NamedTuple(),
    style::NamedTuple = NamedTuple(),
    highlight::AbstractVector = NamedTuple[]
)
    margin = 48
    canvasOptions = graphFigureCanvasOptions(canvas)
    layoutOptions = graphFigureLayoutOptions(layout)
    nodeOptions = graphFigureNodeOptions(node)
    labelOptions = graphFigureLabelOptions(label)
    viewOptions = graphFigureViewOptions(tree.graph, view)
    graphStyle = graphFigureStyle(style)

    return treeGraphFigure(
        tree;
        width = canvasOptions.width,
        height = canvasOptions.height,
        canvasPadding = canvasOptions.padding,
        zoom = canvasOptions.zoom,
        margin = margin,
        rowSpacing = layoutOptions.rowSpacing,
        columnSpacing = layoutOptions.columnSpacing,
        variableRadius = nodeOptions.variableRadius,
        factorSize = nodeOptions.factorSize,
        labelPlacement = labelOptions.placement,
        outsideLabelGap = labelOptions.outsideGap,
        curvedEdges = layoutOptions.curvedEdges,
        orientation = layoutOptions.orientation,
        style = graphStyle,
        highlight = highlight,
        showVariableLabels = labelOptions.showVariables,
        showFactorLabels = labelOptions.showFactors,
        showTooltips = labelOptions.showTooltips,
        showEdgeIds = labelOptions.showEdgeIds,
        tooltipDetail = labelOptions.tooltipDetail,
        fontSize = labelOptions.fontSize,
        view = viewOptions
    )
end

"""
    saveGraphFigure(path::AbstractString, graph::AbstractFactorGraph; kwargs...)

Write a dependency-free SVG visualization of a factor graph to `path`.

# Arguments

- `path`: Output SVG file path.
- `graph`: Factor graph or tree view to draw.

# Keywords

All keyword arguments are forwarded to [`graphFigure`](@ref).

# Returns

The written `path`.

# Notes

The output is an SVG document. Parent directories must already exist.

# Example

```julia
saveGraphFigure(
    "graph.svg",
    graph;
    layout = (orientation = :vertical,),
    label = (placement = :inside,)
)
```
"""
function saveGraphFigure(path::AbstractString, graph::AbstractFactorGraph; kwargs...)
    open(path, "w") do io
        write(io, graphFigure(graph; kwargs...))
    end

    return path
end

function saveGraphFigure(path::AbstractString, tree::TreeFactorGraph; kwargs...)
    open(path, "w") do io
        write(io, graphFigure(tree; kwargs...))
    end

    return path
end

function treeGraphFigure(
    tree::TreeFactorGraph;
    width::Union{Nothing, Int},
    height::Union{Nothing, Int},
    canvasPadding::Int,
    zoom::Real,
    margin::Int,
    rowSpacing,
    columnSpacing,
    variableRadius::Int,
    factorSize::Int,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    curvedEdges::Bool,
    orientation::Symbol,
    style::NamedTuple,
    highlight::AbstractVector,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    showTooltips::Bool,
    showEdgeIds::Bool,
    tooltipDetail::Symbol,
    fontSize::Int,
    view::NamedTuple
)
    if !(orientation in (:horizontal, :vertical))
        error("orientation must be :horizontal or :vertical.")
    end

    if zoom <= 0
        error("zoom must be positive.")
    end

    if canvasPadding < 0
        error("canvasPadding must be non-negative.")
    end

    if orientation == :vertical
        return verticalTreeGraphFigure(
            tree;
            width = width,
            height = height,
            canvasPadding = canvasPadding,
            zoom = zoom,
            margin = margin,
            rowSpacing = rowSpacing,
            columnSpacing = columnSpacing,
            variableRadius = variableRadius,
            factorSize = factorSize,
            labelPlacement = labelPlacement,
            outsideLabelGap = outsideLabelGap,
            curvedEdges = curvedEdges,
            style = style,
            highlight = highlight,
            showVariableLabels = showVariableLabels,
            showFactorLabels = showFactorLabels,
            showTooltips = showTooltips,
            showEdgeIds = showEdgeIds,
            tooltipDetail = tooltipDetail,
            fontSize = fontSize,
            view = view
        )
    end

    graph = tree.graph

    if !(labelPlacement in (:inside, :outside))
        error("labelPlacement must be :inside or :outside.")
    end

    variableDepths, factorDepths = treeNodeDepths(tree)
    maxDepth = maximum(vcat(variableDepths[view.variables], factorDepths[view.factors]))
    depthCounts = [
        count(index -> variableDepths[index] == depth, view.variables) +
            count(index -> factorDepths[index] == depth, view.factors)
        for depth in 1:maxDepth
    ]
    maxDepthCount = maximum(depthCounts; init = 1)
    defaultLayoutHeight = 2 * margin + graphFigureSpacingSum(rowSpacing, maxDepthCount - 1)
    layoutHeight = isnothing(height) ? defaultLayoutHeight : max(height, defaultLayoutHeight)
    minCanvasHeight = isnothing(height) ? 0 : height

    layoutNodeSpacing = rowSpacing
    defaultLayoutWidth = 2 * margin + graphFigureSpacingSum(columnSpacing, maxDepth - 1)
    layoutWidth = isnothing(width) ? defaultLayoutWidth : max(width, defaultLayoutWidth)
    minCanvasWidth = isnothing(width) ? 0 : width
    levelX = graphFigurePositions(maxDepth, layoutWidth, columnSpacing)

    variableRadii = nodeRadii(
        graph.variables,
        variableRadius,
        labelPlacement,
        showVariableLabels,
        fontSize
    )
    factorWidths = factorNodeWidths(
        graph.factors,
        factorSize,
        labelPlacement,
        showFactorLabels,
        fontSize
    )
    factorHeights = factorNodeHeights(
        graph.factors,
        factorSize,
        labelPlacement,
        showFactorLabels
    )
    variableX = zeros(Float64, length(graph.variables))
    factorX = zeros(Float64, length(graph.factors))
    for index in view.variables
        variableX[index] = levelX[variableDepths[index]]
    end
    for index in view.factors
        factorX[index] = levelX[factorDepths[index]]
    end
    variableY = zeros(Float64, length(graph.variables))
    factorY = zeros(Float64, length(graph.factors))

    for depth in 1:maxDepth
        nodes = Tuple{Symbol, Int}[]

        for index in view.variables
            if variableDepths[index] == depth
                push!(nodes, (:variable, index))
            end
        end

        for index in view.factors
            if factorDepths[index] == depth
                push!(nodes, (:factor, index))
            end
        end

        rows = graphFigureCenteredRows(
            length(nodes),
            maxDepthCount,
            layoutHeight,
            layoutNodeSpacing
        )

        for (position, (kind, index)) in pairs(nodes)
            if kind == :variable
                variableY[index] = rows[position]
            else
                factorY[index] = rows[position]
            end
        end
    end

    canvasHeight = fitVerticalCanvas!(
        variableY,
        factorY,
        variableRadii,
        factorHeights,
        minCanvasHeight,
        labelPlacement,
        showVariableLabels,
        showFactorLabels,
        outsideLabelGap,
        fontSize,
        canvasPadding;
        variableIndices = view.variables,
        factorIndices = view.factors
    )
    canvasWidth = fitSideLabelHorizontalCanvas!(
        variableX,
        factorX,
        variableRadii,
        factorWidths,
        minCanvasWidth,
        graph.variables,
        graph.factors,
        labelPlacement,
        showVariableLabels,
        showFactorLabels,
        outsideLabelGap,
        fontSize,
        canvasPadding;
        variableIndices = view.variables,
        factorIndices = view.factors
    )

    highlightState = graphFigureHighlightState(
        graph,
        highlight,
        style
    )
    buffer = startGraphFigureBuffer(
        canvasWidth,
        canvasHeight,
        zoom,
        fontSize,
        style
    )

    occupiedEdgeLabels = NamedTuple[]
    edgeLabelBackgrounds = String[]
    edgeLabelTexts = String[]
    edgeGeometries = NamedTuple[]

    for edge in view.edges
        factorLeft = factorX[edge.factorIndex] - factorWidths[edge.factorIndex] / 2
        factorRight = factorX[edge.factorIndex] + factorWidths[edge.factorIndex] / 2

        if variableX[edge.variableIndex] <= factorX[edge.factorIndex]
            x1 = variableX[edge.variableIndex] + variableRadii[edge.variableIndex]
            y1 = variableY[edge.variableIndex]
            x2 = factorLeft
            y2 = factorY[edge.factorIndex]
        else
            x1 = factorRight
            y1 = factorY[edge.factorIndex]
            x2 = variableX[edge.variableIndex] - variableRadii[edge.variableIndex]
            y2 = variableY[edge.variableIndex]
        end

        push!(
            edgeGeometries,
            edgeGeometry(edge.id, x1, y1, x2, y2, curvedEdges, :horizontal)
        )
    end

    for edge in view.edges
        factorLeft = factorX[edge.factorIndex] - factorWidths[edge.factorIndex] / 2
        factorRight = factorX[edge.factorIndex] + factorWidths[edge.factorIndex] / 2

        if variableX[edge.variableIndex] <= factorX[edge.factorIndex]
            x1 = variableX[edge.variableIndex] + variableRadii[edge.variableIndex]
            y1 = variableY[edge.variableIndex]
            x2 = factorLeft
            y2 = factorY[edge.factorIndex]
        else
            x1 = factorRight
            y1 = factorY[edge.factorIndex]
            x2 = variableX[edge.variableIndex] - variableRadii[edge.variableIndex]
            y2 = variableY[edge.variableIndex]
        end

        edgeClass = edgeClassName(
            edge,
            highlightState.edges
        )
        tooltip = showTooltips ? graphFigureEdgeTooltip(graph, edge, tooltipDetail) : nothing
        println(
            buffer,
            svgEdge(
                x1,
                y1,
                x2,
                y2,
                curvedEdges,
                edgeClass,
                get(highlightState.edgeStyles, edge.id, nothing),
                tooltip
            )
        )
        if showEdgeIds
            labelParts = svgEdgeIdLabelParts(
                edge.id,
                x1,
                y1,
                x2,
                y2,
                occupiedEdgeLabels,
                curvedEdges,
                :horizontal,
                edgeGeometries,
                tooltip
            )
            push!(edgeLabelBackgrounds, labelParts.background)
            push!(edgeLabelTexts, labelParts.text)
        end
    end

    drawEdgeIdLabels!(buffer, edgeLabelBackgrounds, edgeLabelTexts)

    for index in view.variables
        variable = graph.variables[index]
        label = escapeXML(variable.label)
        className = variableClassName(index, highlightState.variables, view.focusVariableSet)
        tooltip = showTooltips ? graphFigureVariableTooltip(graph, index, tooltipDetail) : nothing
        println(
            buffer,
            svgCircle(
                variableX[index],
                variableY[index],
                variableRadii[index],
                className,
                get(highlightState.variableStyles, index, nothing),
                tooltip
            )
        )

        if showVariableLabels
            drawVariableLabel!(
                buffer,
                variableX[index],
                variableY[index],
                variableRadii[index],
                label,
                labelPlacement,
                outsideLabelGap,
                fontSize,
                tooltip
            )
        end
    end

    for index in view.factors
        factor = graph.factors[index]
        label = escapeXML(factor.label)
        className = factorClassName(index, highlightState.factors, view.focusFactorSet)
        x = factorX[index] - factorWidths[index] / 2
        y = factorY[index] - factorHeights[index] / 2
        println(
            buffer,
            svgRect(
                x,
                y,
                factorWidths[index],
                factorHeights[index],
                labelPlacement,
                className,
                get(highlightState.factorStyles, index, nothing),
                showTooltips ? graphFigureFactorTooltip(graph, index, tooltipDetail) : nothing
            )
        )

        if showFactorLabels
            drawTreeFactorLabel!(
                buffer,
                factorX[index],
                factorY[index],
                factorHeights[index],
                label,
                labelPlacement,
                outsideLabelGap,
                fontSize,
                showTooltips ? graphFigureFactorTooltip(graph, index, tooltipDetail) : nothing
            )
        end
    end

    println(buffer, "</svg>")

    return String(take!(buffer))
end

function verticalTreeGraphFigure(
    tree::TreeFactorGraph;
    width::Union{Nothing, Int},
    height::Union{Nothing, Int},
    canvasPadding::Int,
    zoom::Real,
    margin::Int,
    rowSpacing,
    columnSpacing,
    variableRadius::Int,
    factorSize::Int,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    curvedEdges::Bool,
    style::NamedTuple,
    highlight::AbstractVector,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    showTooltips::Bool,
    showEdgeIds::Bool,
    tooltipDetail::Symbol,
    fontSize::Int,
    view::NamedTuple
)
    if zoom <= 0
        error("zoom must be positive.")
    end

    if canvasPadding < 0
        error("canvasPadding must be non-negative.")
    end

    graph = tree.graph
    variableDepths, factorDepths = treeNodeDepths(tree)
    maxDepth = maximum(vcat(variableDepths[view.variables], factorDepths[view.factors]))
    depthCounts = [
        count(index -> variableDepths[index] == depth, view.variables) +
            count(index -> factorDepths[index] == depth, view.factors)
        for depth in 1:maxDepth
    ]
    maxDepthCount = maximum(depthCounts; init = 1)
    minCanvasWidth = isnothing(width) ? 0 : width
    defaultLayoutHeight = 2 * margin + graphFigureSpacingSum(rowSpacing, maxDepth - 1)
    layoutHeight = isnothing(height) ? defaultLayoutHeight : max(height, defaultLayoutHeight)
    minCanvasHeight = isnothing(height) ? 0 : height

    levelY = graphFigurePositions(maxDepth, layoutHeight, rowSpacing)

    variableRadii = nodeRadii(
        graph.variables,
        variableRadius,
        labelPlacement,
        showVariableLabels,
        fontSize
    )
    factorWidths = factorNodeWidths(
        graph.factors,
        factorSize,
        labelPlacement,
        showFactorLabels,
        fontSize
    )
    factorHeights = factorNodeHeights(
        graph.factors,
        factorSize,
        labelPlacement,
        showFactorLabels
    )
    variableX = zeros(Float64, length(graph.variables))
    factorX = zeros(Float64, length(graph.factors))
    variableY = zeros(Float64, length(graph.variables))
    factorY = zeros(Float64, length(graph.factors))
    for index in view.variables
        variableY[index] = levelY[variableDepths[index]]
    end
    for index in view.factors
        factorY[index] = levelY[factorDepths[index]]
    end

    for depth in 1:maxDepth
        nodes = Tuple{Symbol, Int}[]

        for index in view.variables
            if variableDepths[index] == depth
                push!(nodes, (:variable, index))
            end
        end

        for index in view.factors
            if factorDepths[index] == depth
                push!(nodes, (:factor, index))
            end
        end

        levelSpacing = treeLevelSpacing(
            graph,
            nodes,
            factorWidths,
            variableRadii,
            columnSpacing,
            outsideLabelGap,
            labelPlacement,
            showVariableLabels,
            showFactorLabels,
            fontSize
        )
        columns = nodeRows(length(nodes), maxDepthCount, margin, levelSpacing)

        for (position, (kind, index)) in pairs(nodes)
            if kind == :variable
                variableX[index] = columns[position]
            else
                factorX[index] = columns[position]
            end
        end
    end

    canvasWidth = fitHorizontalCanvas(
        variableX,
        factorX,
        variableRadii,
        factorWidths,
        minCanvasWidth,
        graph.variables,
        graph.factors,
        labelPlacement,
        showVariableLabels,
        showFactorLabels,
        outsideLabelGap,
        fontSize,
        canvasPadding;
        variableIndices = view.variables,
        factorIndices = view.factors
    )
    canvasHeight = fitVerticalCanvas!(
        variableY,
        factorY,
        variableRadii,
        factorHeights,
        minCanvasHeight,
        labelPlacement,
        showVariableLabels,
        showFactorLabels,
        outsideLabelGap,
        fontSize,
        canvasPadding;
        variableIndices = view.variables,
        factorIndices = view.factors
    )

    highlightState = graphFigureHighlightState(
        graph,
        highlight,
        style
    )
    buffer = startGraphFigureBuffer(
        canvasWidth,
        canvasHeight,
        zoom,
        fontSize,
        style
    )

    occupiedEdgeLabels = NamedTuple[]
    edgeLabelBackgrounds = String[]
    edgeLabelTexts = String[]
    edgeGeometries = NamedTuple[]

    for edge in view.edges
        if variableY[edge.variableIndex] <= factorY[edge.factorIndex]
            x1 = variableX[edge.variableIndex]
            y1 = variableY[edge.variableIndex] + variableRadii[edge.variableIndex]
            x2 = factorX[edge.factorIndex]
            y2 = factorY[edge.factorIndex] - factorHeights[edge.factorIndex] / 2
        else
            x1 = factorX[edge.factorIndex]
            y1 = factorY[edge.factorIndex] + factorHeights[edge.factorIndex] / 2
            x2 = variableX[edge.variableIndex]
            y2 = variableY[edge.variableIndex] - variableRadii[edge.variableIndex]
        end

        push!(
            edgeGeometries,
            edgeGeometry(edge.id, x1, y1, x2, y2, curvedEdges, :vertical)
        )
    end

    for edge in view.edges
        if variableY[edge.variableIndex] <= factorY[edge.factorIndex]
            x1 = variableX[edge.variableIndex]
            y1 = variableY[edge.variableIndex] + variableRadii[edge.variableIndex]
            x2 = factorX[edge.factorIndex]
            y2 = factorY[edge.factorIndex] - factorHeights[edge.factorIndex] / 2
        else
            x1 = factorX[edge.factorIndex]
            y1 = factorY[edge.factorIndex] + factorHeights[edge.factorIndex] / 2
            x2 = variableX[edge.variableIndex]
            y2 = variableY[edge.variableIndex] - variableRadii[edge.variableIndex]
        end

        edgeClass = edgeClassName(
            edge,
            highlightState.edges
        )
        tooltip = showTooltips ? graphFigureEdgeTooltip(graph, edge, tooltipDetail) : nothing
        println(
            buffer,
            svgEdge(
                x1,
                y1,
                x2,
                y2,
                curvedEdges,
                edgeClass,
                get(highlightState.edgeStyles, edge.id, nothing),
                tooltip,
                :vertical
            )
        )
        if showEdgeIds
            labelParts = svgEdgeIdLabelParts(
                edge.id,
                x1,
                y1,
                x2,
                y2,
                occupiedEdgeLabels,
                curvedEdges,
                :vertical,
                edgeGeometries,
                tooltip
            )
            push!(edgeLabelBackgrounds, labelParts.background)
            push!(edgeLabelTexts, labelParts.text)
        end
    end

    drawEdgeIdLabels!(buffer, edgeLabelBackgrounds, edgeLabelTexts)

    for index in view.variables
        variable = graph.variables[index]
        label = escapeXML(variable.label)
        className = variableClassName(index, highlightState.variables, view.focusVariableSet)
        tooltip = showTooltips ? graphFigureVariableTooltip(graph, index, tooltipDetail) : nothing
        println(
            buffer,
            svgCircle(
                variableX[index],
                variableY[index],
                variableRadii[index],
                className,
                get(highlightState.variableStyles, index, nothing),
                tooltip
            )
        )

        if showVariableLabels
            drawSideLabel!(
                buffer,
                variableX[index] + variableRadii[index] + outsideLabelGap,
                variableY[index],
                label,
                labelPlacement,
                :right,
                tooltip
            )
        end
    end

    for index in view.factors
        factor = graph.factors[index]
        label = escapeXML(factor.label)
        className = factorClassName(index, highlightState.factors, view.focusFactorSet)
        x = factorX[index] - factorWidths[index] / 2
        y = factorY[index] - factorHeights[index] / 2
        println(
            buffer,
            svgRect(
                x,
                y,
                factorWidths[index],
                factorHeights[index],
                labelPlacement,
                className,
                get(highlightState.factorStyles, index, nothing),
                showTooltips ? graphFigureFactorTooltip(graph, index, tooltipDetail) : nothing
            )
        )

        if showFactorLabels
            drawSideLabel!(
                buffer,
                factorX[index] + factorWidths[index] / 2 + outsideLabelGap,
                factorY[index],
                label,
                labelPlacement,
                :right,
                showTooltips ? graphFigureFactorTooltip(graph, index, tooltipDetail) : nothing
            )
        end
    end

    println(buffer, "</svg>")

    return String(take!(buffer))
end

function treeNodeDepths(tree::TreeFactorGraph)
    graph = tree.graph
    variableDepths = zeros(Int, length(graph.variables))
    factorDepths = zeros(Int, length(graph.factors))
    variableDepths[tree.rootVariableIndex] = 1

    for edgeId in tree.backwardOrder
        edge = graph.edges[edgeId]

        if tree.factorParentEdge[edge.factorIndex] == edgeId
            factorDepths[edge.factorIndex] = variableDepths[edge.variableIndex] + 1
        elseif tree.variableParentEdge[edge.variableIndex] == edgeId
            variableDepths[edge.variableIndex] = factorDepths[edge.factorIndex] + 1
        end
    end

    return variableDepths, factorDepths
end

function verticalGraphColumnSpacing(
    graph::AbstractFactorGraph,
    variableIndices::AbstractVector{Int},
    variableRadii::AbstractVector,
    columnSpacing,
    outsideLabelGap::Int,
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    fontSize::Int
)
    baseColumnSpacing = graphFigureSpacingMaximum(columnSpacing)

    if labelPlacement == :inside || !showVariableLabels
        return columnSpacing
    end

    maxWidth = 0.0

    for index in variableIndices
        labelWidth = approximateTextWidth(graph.variables[index].label, fontSize)
        maxWidth = max(maxWidth, 2 * variableRadii[index] + labelWidth)
    end

    requiredSpacing = ceil(Int, maxWidth + 2 * outsideLabelGap)

    return requiredSpacing <= baseColumnSpacing ? columnSpacing : requiredSpacing
end

function treeLevelSpacing(
    graph::AbstractFactorGraph,
    nodes::AbstractVector{Tuple{Symbol, Int}},
    factorWidths::AbstractVector,
    variableRadii::AbstractVector,
    columnSpacing,
    outsideLabelGap::Int,
    labelPlacement::Symbol,
    showVariableLabels::Bool,
    showFactorLabels::Bool,
    fontSize::Int
)
    baseColumnSpacing = graphFigureSpacingMaximum(columnSpacing)

    if labelPlacement == :inside
        return columnSpacing
    end

    maxWidth = 0.0

    for (kind, index) in nodes
        if kind == :variable
            labelWidth = showVariableLabels ?
                outsideLabelGap + approximateTextWidth(graph.variables[index].label, fontSize) : 0.0
            maxWidth = max(maxWidth, 2 * variableRadii[index] + labelWidth)
        else
            labelWidth = showFactorLabels ?
                outsideLabelGap + approximateTextWidth(graph.factors[index].label, fontSize) : 0.0
            maxWidth = max(maxWidth, factorWidths[index] + labelWidth)
        end
    end

    requiredSpacing = ceil(Int, maxWidth + outsideLabelGap)

    return requiredSpacing <= baseColumnSpacing ? columnSpacing : requiredSpacing
end


function isUnaryFactor(graph::AbstractFactorGraph, factorIndex::Int)
    return length(graph.factorEdges[factorIndex]) == 1
end

function edgeClassName(
    edge::Edge,
    highlightedEdges::Set{Int}
)
    highlighted = edge.id in highlightedEdges

    return highlighted ? "edge edge-highlight" : "edge"
end

function variableClassName(
    variableIndex::Int,
    highlightedVariables::Set{Int},
    focusVariables::Set{Int}
)
    if variableIndex in highlightedVariables
        return "variable variable-highlight"
    elseif variableIndex in focusVariables
        return "variable"
    end

    return "variable variable-context"
end

function factorClassName(
    factorIndexValue::Int,
    highlightedFactors::Set{Int},
    focusFactors::Set{Int}
)
    if factorIndexValue in highlightedFactors
        return "factor factor-highlight"
    elseif factorIndexValue in focusFactors
        return "factor"
    end

    return "factor factor-context"
end

function graphFigureViewOptions(graph::AbstractFactorGraph, view::NamedTuple)
    defaultView = (
        variables = nothing,
        factors = nothing,
        hops = 1
    )

    for key in propertynames(view)
        if !(key in propertynames(defaultView))
            error("Unknown graph figure view key $key.")
        end
    end

    options = merge(defaultView, view)

    if !((options.hops isa Integer && options.hops >= 0) || options.hops == :all)
        error("view.hops must be a nonnegative integer or :all.")
    end

    if isnothing(options.variables) && isnothing(options.factors)
        variableIndices = collect(eachindex(graph.variables))
        factorIndices = collect(eachindex(graph.factors))
        visibleEdges = graph.edges

        return (
            variables = variableIndices,
            focusVariableSet = Set(variableIndices),
            variableSet = Set(variableIndices),
            factors = factorIndices,
            focusFactorSet = Set(factorIndices),
            factorSet = Set(factorIndices),
            edges = visibleEdges,
            hops = options.hops
        )
    end

    variableIndices = isnothing(options.variables) ? Int[] : [
        graphFigureVariableIndex(graph, variable) for variable in options.variables
    ]
    factorIndices = isnothing(options.factors) ? Int[] : [
        graphFigureFactorIndex(graph, factor) for factor in options.factors
    ]

    if isempty(variableIndices) && isempty(factorIndices)
        error("view.variables or view.factors must contain at least one node.")
    end

    if length(unique(variableIndices)) != length(variableIndices)
        error("view.variables must not contain duplicate variables.")
    end

    if length(unique(factorIndices)) != length(factorIndices)
        error("view.factors must not contain duplicate factors.")
    end

    startVariableSet = Set(variableIndices)
    startFactorSet = Set(factorIndices)
    variableSet = Set(variableIndices)
    factorSet = Set(factorIndices)
    variableDepths = Dict(index => 0 for index in variableIndices)
    factorDepths = Dict(index => 0 for index in factorIndices)
    frontierVariables = Set(variableIndices)
    frontierFactors = Set(factorIndices)
    hopLimit = options.hops == :all ?
        length(graph.variables) + length(graph.factors) :
        options.hops

    for hop in 1:hopLimit
        nextFactors = Set{Int}()
        nextVariables = Set{Int}()

        for edge in graph.edges
            if edge.variableIndex in frontierVariables && !(edge.factorIndex in factorSet)
                push!(nextFactors, edge.factorIndex)
                factorDepths[edge.factorIndex] = hop
            end

            if edge.factorIndex in frontierFactors && !(edge.variableIndex in variableSet)
                push!(nextVariables, edge.variableIndex)
                variableDepths[edge.variableIndex] = hop
            end
        end

        union!(factorSet, nextFactors)
        union!(variableSet, nextVariables)
        frontierVariables = nextVariables
        frontierFactors = nextFactors

        if isempty(frontierVariables) && isempty(frontierFactors)
            break
        end
    end

    expandedVariables = [
        index for index in eachindex(graph.variables)
        if index in variableSet && !(index in startVariableSet)
    ]
    sort!(expandedVariables; by = index -> (variableDepths[index], index))
    variableIndices = vcat(variableIndices, expandedVariables)
    expandedFactors = [
        index for index in eachindex(graph.factors)
        if index in factorSet && !(index in startFactorSet)
    ]
    sort!(expandedFactors; by = index -> (factorDepths[index], index))
    factorIndices = vcat(factorIndices, expandedFactors)
    visibleEdges = [
        edge for edge in graph.edges
        if edge.variableIndex in variableSet && edge.factorIndex in factorSet
    ]
    focusVariableSet = copy(startVariableSet)
    focusFactorSet = copy(startFactorSet)

    for edge in visibleEdges
        if edge.variableIndex in startVariableSet
            push!(focusFactorSet, edge.factorIndex)
        end

        if edge.factorIndex in startFactorSet
            push!(focusVariableSet, edge.variableIndex)
        end
    end

    return (
        variables = variableIndices,
        focusVariableSet = focusVariableSet,
        variableSet = variableSet,
        factors = factorIndices,
        focusFactorSet = focusFactorSet,
        factorSet = factorSet,
        edges = visibleEdges,
        hops = options.hops
    )
end

function graphFigureVariableTooltip(
    graph::AbstractFactorGraph,
    variableIndex::Int,
    detail::Symbol
)
    variable = graph.variables[variableIndex]
    identityLines = [
        "type: $(graphFigureModelKind(variable))",
        "index: $variableIndex",
        "id: $(graphFigureValueText(variable.id))",
        graphFigureVariableSizeTooltip(variable),
        "degree: $(length(graph.variableEdges[variableIndex]))"
    ]
    lines = graphFigureTooltipLines(
        "Variable $(variable.label)",
        identityLines,
        graphFigureVariableDataTooltip(variable, detail),
        detail
    )

    return escapeXML(join(lines, "\n"))
end

function graphFigureFactorTooltip(
    graph::AbstractFactorGraph,
    factorIndex::Int,
    detail::Symbol
)
    factor = graph.factors[factorIndex]
    variableLabels = [
        graph.variables[edge.variableIndex].label for edge in graph.edges
        if edge.factorIndex == factorIndex
    ]
    identityLines = [
        "type: $(graphFigureModelKind(factor))",
        "index: $factorIndex",
        "variables: $(join(variableLabels, ", "))",
        "degree: $(length(graph.factorEdges[factorIndex]))"
    ]
    lines = graphFigureTooltipLines(
        "Factor $(factor.label)",
        identityLines,
        graphFigureFactorDataTooltip(graph, factorIndex, detail),
        detail
    )

    return escapeXML(join(lines, "\n"))
end

graphFigureModelKind(::Union{GaussianVariable, GaussianFactor}) = "Gaussian"
graphFigureModelKind(::Union{DiscreteVariable, DiscreteFactor}) = "Discrete"

function graphFigureEdgeTooltip(graph::AbstractFactorGraph, edge::Edge, detail::Symbol)
    variable = graph.variables[edge.variableIndex]
    factor = graph.factors[edge.factorIndex]
    identityLines = [
        "variable: $(variable.label)",
        "factor: $(factor.label)",
        "variable index: $(edge.variableIndex)",
        "factor index: $(edge.factorIndex)"
    ]
    lines = graphFigureTooltipLines(
        "Edge $(edge.id)",
        identityLines,
        Pair{String, Vector{String}}[],
        detail
    )

    return escapeXML(join(lines, "\n"))
end

function graphFigureVariableSizeTooltip(variable::GaussianVariable)
    return "dimension: $(variable.dimension)"
end

function graphFigureVariableSizeTooltip(variable::DiscreteVariable)
    return "cardinality: $(variable.cardinality)"
end

function graphFigureVariableDataTooltip(variable::GaussianVariable, detail::Symbol)
    if detail == :full
        return graphFigureFullVariableDataTooltip(variable)
    end

    return [
        "Initialization" => [
            "mean size: $(graphFigureOptionalSizeText(variable.mean))",
            "covariance size: $(graphFigureOptionalSizeText(variable.covariance))"
        ],
        "Components" => [
            "component count: $(length(variable.components))"
        ]
    ]
end

function graphFigureVariableDataTooltip(variable::DiscreteVariable, detail::Symbol)
    if detail == :full
        return [
            "Initialization" => [
                graphFigureOptionalFieldText("probability", variable.probability)
            ],
            "States" => [
                "states: $(graphFigureReferenceVectorText(variable.states))"
            ]
        ]
    end

    return [
        "Initialization" => [
            "probability size: $(graphFigureOptionalSizeText(variable.probability))"
        ],
        "States" => [
            "state count: $(length(variable.states))"
        ]
    ]
end

function graphFigureFullVariableDataTooltip(variable::GaussianVariable)
    if variable.dimension == 1
        return [
            "Initialization" => [
                graphFigureOptionalFieldText("mean", variable.mean),
                "variance: $(graphFigureScalarVarianceText(variable.covariance))"
            ],
            "Components" => [
                "component: $(graphFigureValueText(only(variable.components)))"
            ]
        ]
    end

    return [
        "Initialization" => [
            graphFigureOptionalFieldText("mean", variable.mean),
            graphFigureOptionalFieldText("covariance", variable.covariance)
        ],
        "Components" => [
            "components: $(graphFigureReferenceVectorText(variable.components))"
        ]
    ]
end

function graphFigureFactorDataTooltip(graph::AbstractFactorGraph, factorIndex::Int, detail::Symbol)
    return graphFigureFactorDataTooltip(
        graph,
        graph.factors[factorIndex],
        detail,
        length(graph.factorEdges[factorIndex]) == 1
    )
end

function graphFigureFactorDataTooltip(
    graph::AbstractFactorGraph,
    factor::GaussianFactor,
    detail::Symbol,
    isUnary::Bool
)
    if detail == :full
        sections = [
            "Parameters" => [
                graphFigureFieldText("mean", factor.mean),
                graphFigureFieldText("coefficient", factor.coefficient),
                graphFigureFieldText("covariance", factor.covariance)
            ]
        ]

        if isUnary
            push!(sections, "Options" => ["initialize: $(factor.initialize)"])
        end

        return sections
    end

    sections = [
        "Parameters" => [
            "mean size: $(graphFigureSizeText(factor.mean))",
            "coefficient size: $(graphFigureSizeText(factor.coefficient))",
            "covariance size: $(graphFigureSizeText(factor.covariance))"
        ]
    ]

    if isUnary
        push!(sections, "Options" => ["initialize: $(factor.initialize)"])
    end

    return sections
end

function graphFigureFactorDataTooltip(
    graph::AbstractFactorGraph,
    factor::DiscreteFactor,
    detail::Symbol,
    isUnary::Bool
)
    if detail == :full
        sections = [
            "Table" => graphFigureDiscreteTableTooltipLines(graph, factor)
        ]

        if isUnary
            push!(sections, "Options" => ["initialize: $(factor.initialize)"])
        end

        return sections
    end

    sections = [
        "Table" => ["size: $(graphFigureSizeText(factor.table))"]
    ]

    if isUnary
        push!(sections, "Options" => ["initialize: $(factor.initialize)"])
    end

    return sections
end

function graphFigureTooltipLines(
    title::AbstractString,
    identityLines::AbstractVector{String},
    sections,
    detail::Symbol
)
    lines = [title, "", "Identity"]
    append!(lines, graphFigureIndentedTooltipLines(identityLines))

    for section in sections
        if !isempty(section.second)
            push!(lines, "")
            push!(lines, section.first)
            append!(lines, graphFigureIndentedTooltipLines(section.second))
        end
    end

    return lines
end

function graphFigureIndentedTooltipLines(sectionLines::AbstractVector{String})
    lines = String[]

    for line in sectionLines
        for part in split(line, "\n"; keepempty = true)
            push!(lines, isempty(part) ? "" : "  $part")
        end
    end

    return lines
end

function graphFigureOptionalFieldText(name::AbstractString, value)
    return isnothing(value) ? "$name: not set" : graphFigureFieldText(name, value)
end

function graphFigureFieldText(name::AbstractString, value)
    if value isa AbstractArray
        return "$name:\n$(graphFigureArrayText(value))"
    end

    return "$name: $(graphFigureValueText(value))"
end

function graphFigureDiscreteTableTooltipLines(graph::AbstractFactorGraph, factor::DiscreteFactor)
    variableLabels = [
        graph.variables[variableIndex(graph.referenceIndex, variableRef)].label
        for variableRef in factor.variables
    ]
    lines = ["size: $(graphFigureSizeText(factor.table))"]

    if ndims(factor.table) == 3 && length(factor.variables) >= 3
        sliceVariable = graph.variables[variableIndex(graph.referenceIndex, factor.variables[3])]
        rowLabel = variableLabels[1]
        columnLabel = variableLabels[2]

        push!(lines, "axes: $(join(variableLabels, " x "))")
        append!(
            lines,
            split(
                graphFigureThreeDimensionalArrayText(
                    factor.table,
                    sliceVariable,
                    rowLabel,
                    columnLabel;
                    includeSize = false
                ),
                "\n"
            )
        )

        return lines
    end

    axisLines = graphFigureDiscreteTableAxisLines(factor.table, variableLabels)
    append!(lines, axisLines)
    append!(lines, split(graphFigureArrayText(factor.table), "\n"))

    return lines
end

function graphFigureDiscreteTableAxisLines(array::AbstractArray, variableLabels::AbstractVector)
    if ndims(array) == 1 && length(variableLabels) >= 1
        return ["axis: $(variableLabels[1])"]
    elseif ndims(array) == 2 && length(variableLabels) >= 2
        return ["rows: $(variableLabels[1]), columns: $(variableLabels[2])"]
    end

    return String[]
end

function graphFigureScalarVarianceText(matrix::Nothing)
    return "not set"
end

function graphFigureScalarVarianceText(matrix::AbstractMatrix)
    return graphFigureValueText(matrix[1, 1])
end

function graphFigureOptionalSizeText(value)
    return isnothing(value) ? "not set" : graphFigureSizeText(value)
end

function graphFigureSizeText(array::AbstractArray)
    return join(size(array), " x ")
end

function graphFigureReferenceVectorText(values::AbstractVector)
    return "[$(join((graphFigureValueText(value) for value in values), ", "))]"
end

function graphFigureValueText(value)
    return sprint(show, value)
end

function graphFigureArrayText(array::AbstractVector)
    if isempty(array)
        return "  []"
    end

    return "  [$(join((graphFigureValueText(value) for value in array), ", "))]"
end

function graphFigureArrayText(array::AbstractMatrix)
    return graphFigureMatrixText(array, 2)
end

function graphFigureArrayText(array::AbstractArray)
    if ndims(array) == 0
        return "  $(graphFigureValueText(only(array)))"
    elseif ndims(array) == 1
        return graphFigureArrayText(vec(array))
    elseif ndims(array) == 2
        return graphFigureMatrixText(array, 2)
    elseif ndims(array) == 3
        return graphFigureThreeDimensionalArrayText(array)
    end

    return "  size: $(graphFigureSizeText(array))\n  $(graphFigureValueText(array))"
end

function graphFigureThreeDimensionalArrayText(array::AbstractArray)
    lines = ["  size: $(graphFigureSizeText(array))"]

    for index in axes(array, 3)
        push!(lines, "")
        push!(lines, "  slice 3 = $index")
        append!(lines, split(graphFigureMatrixText(@view(array[:, :, index]), 4), "\n"))
    end

    return join(lines, "\n")
end

function graphFigureThreeDimensionalArrayText(
    array::AbstractArray,
    sliceVariable::DiscreteVariable,
    rowLabel::AbstractString,
    columnLabel::AbstractString;
    includeSize::Bool = true
)
    lines = includeSize ? ["  size: $(graphFigureSizeText(array))"] : String[]
    labelIndent = includeSize ? "  " : ""
    matrixIndent = includeSize ? 4 : 2

    for index in axes(array, 3)
        push!(lines, "")
        push!(
            lines,
            "$(labelIndent)slice $(sliceVariable.label) = " *
                graphFigureValueText(sliceVariable.states[index])
        )
        push!(lines, "$(labelIndent)rows: $rowLabel, columns: $columnLabel")
        append!(lines, split(graphFigureMatrixText(@view(array[:, :, index]), matrixIndent), "\n"))
    end

    return join(lines, "\n")
end

function graphFigureMatrixText(matrix::AbstractMatrix, indent::Int)
    spaces = repeat(" ", indent)

    if isempty(matrix)
        return "$(spaces)[]"
    end

    cellText = [
        graphFigureValueText(matrix[row, column])
        for row in axes(matrix, 1), column in axes(matrix, 2)
    ]
    widths = [
        maximum(length(cellText[row, column]) for row in axes(cellText, 1))
        for column in axes(cellText, 2)
    ]
    lines = String[]

    for row in axes(cellText, 1)
        cells = [
            lpad(cellText[row, column], widths[column])
            for column in axes(cellText, 2)
        ]
        leftBracket = row == first(axes(cellText, 1)) ? "[" : " "
        rightBracket = row == last(axes(cellText, 1)) ? "]" : ""
        push!(lines, "$spaces$leftBracket$(join(cells, "  "))$rightBracket")
    end

    return join(lines, "\n")
end

function graphFigureHighlightState(
    graph::AbstractFactorGraph,
    highlight::AbstractVector,
    style::NamedTuple
)
    variableStyles = Dict{Int, NamedTuple}()
    factorStyles = Dict{Int, NamedTuple}()
    edgeStyles = Dict{Int, NamedTuple}()

    applyNamedTupleHighlights!(
        variableStyles,
        factorStyles,
        edgeStyles,
        graph,
        highlight,
        style
    )

    return (
        variableStyles = variableStyles,
        factorStyles = factorStyles,
        edgeStyles = edgeStyles,
        variables = Set(keys(variableStyles)),
        factors = Set(keys(factorStyles)),
        edges = Set(keys(edgeStyles))
    )
end

function applyNamedTupleHighlights!(
    variableStyles::AbstractDict{Int},
    factorStyles::AbstractDict{Int},
    edgeStyles::AbstractDict{Int},
    graph::AbstractFactorGraph,
    highlights::AbstractVector,
    graphStyle::NamedTuple
)
    for highlight in highlights
        if !(highlight isa NamedTuple)
            error("Graph figure highlight entries must be named tuples.")
        end

        hasVariable = hasproperty(highlight, :variable)
        hasFactor = hasproperty(highlight, :factor)
        hasEdge = hasproperty(highlight, :edge)

        if hasEdge && (hasVariable || hasFactor)
            error("A highlight entry cannot combine edge with variable or factor.")
        elseif hasEdge
            edgeIndexValue = graphFigureEdgeIndex(graph, highlight.edge)
            edgeStyles[edgeIndexValue] = graphFigureEdgeHighlightStyle(
                highlight,
                graphStyle.edgeStroke,
                graphStyle.edgeStrokeWidth
            )
        elseif hasVariable && hasFactor
            edgeIndexValue = graphFigureEdgeIndex(graph, highlight.variable, highlight.factor)
            edgeStyles[edgeIndexValue] = graphFigureEdgeHighlightStyle(
                highlight,
                graphStyle.edgeStroke,
                graphStyle.edgeStrokeWidth
            )
        elseif hasVariable
            variableIndexValue = graphFigureVariableIndex(graph, highlight.variable)
            nodeStyle = graphFigureNodeHighlightStyle(
                highlight,
                graphStyle.variableStroke,
                graphStyle.variableFill,
                graphStyle.variableStrokeWidth
            )
            variableStyles[variableIndexValue] = nodeStyle

            if graphFigureHighlightOption(highlight, :incidentEdges, true)
                addVariableIncidentEdgeStyles!(
                    edgeStyles,
                    graph,
                    variableIndexValue,
                    graphFigureIncidentEdgeHighlightStyle(highlight, nodeStyle)
                )
            end
        elseif hasFactor
            factorIndexValue = graphFigureFactorIndex(graph, highlight.factor)
            nodeStyle = graphFigureNodeHighlightStyle(
                highlight,
                graphStyle.factorStroke,
                graphStyle.factorFill,
                graphStyle.factorStrokeWidth
            )
            factorStyles[factorIndexValue] = nodeStyle

            if graphFigureHighlightOption(highlight, :incidentEdges, true)
                addFactorIncidentEdgeStyles!(
                    edgeStyles,
                    graph,
                    factorIndexValue,
                    graphFigureIncidentEdgeHighlightStyle(highlight, nodeStyle)
                )
            end
        else
            error("A highlight entry must contain variable, factor, edge, or variable and factor.")
        end
    end

    return variableStyles, factorStyles, edgeStyles
end

function graphFigureHighlightOption(highlight::NamedTuple, key::Symbol, default)
    return hasproperty(highlight, key) ? getproperty(highlight, key) : default
end

function graphFigureNodeHighlightStyle(
    highlight::NamedTuple,
    defaultStroke,
    defaultFill,
    defaultStrokeWidth
)
    return (
        stroke = something(
            graphFigureHighlightOption(highlight, :stroke, defaultStroke),
            defaultStroke
        ),
        fill = something(graphFigureHighlightOption(highlight, :fill, defaultFill), defaultFill),
        strokeWidth = something(
            graphFigureHighlightOption(highlight, :strokeWidth, defaultStrokeWidth),
            defaultStrokeWidth
        )
    )
end

function graphFigureEdgeHighlightStyle(
    highlight::NamedTuple,
    defaultStroke,
    defaultStrokeWidth
)
    return (
        stroke = something(
            graphFigureHighlightOption(highlight, :stroke, defaultStroke),
            defaultStroke
        ),
        fill = nothing,
        strokeWidth = something(
            graphFigureHighlightOption(highlight, :strokeWidth, defaultStrokeWidth),
            defaultStrokeWidth
        )
    )
end

function graphFigureIncidentEdgeHighlightStyle(highlight::NamedTuple, nodeStyle::NamedTuple)
    return (
        stroke = something(
            graphFigureHighlightOption(highlight, :edgeStroke, nodeStyle.stroke),
            nodeStyle.stroke
        ),
        fill = nothing,
        strokeWidth = something(
            graphFigureHighlightOption(highlight, :edgeStrokeWidth, nodeStyle.strokeWidth),
            nodeStyle.strokeWidth
        )
    )
end

function addVariableIncidentEdgeStyles!(
    edgeStyles::AbstractDict{Int},
    graph::AbstractFactorGraph,
    variableIndexValue::Int,
    style::NamedTuple
)
    for edge in graph.edges
        if edge.variableIndex == variableIndexValue && !haskey(edgeStyles, edge.id)
            edgeStyles[edge.id] = style
        end
    end

    return edgeStyles
end

function addFactorIncidentEdgeStyles!(
    edgeStyles::AbstractDict{Int},
    graph::AbstractFactorGraph,
    factorIndexValue::Int,
    style::NamedTuple
)
    for edge in graph.edges
        if edge.factorIndex == factorIndexValue && !haskey(edgeStyles, edge.id)
            edgeStyles[edge.id] = style
        end
    end

    return edgeStyles
end

function graphFigureVariableIndex(graph::AbstractFactorGraph, variable)
    return variableIndex(graph, variableRef(variable))
end

function graphFigureFactorIndex(graph::AbstractFactorGraph, factor)
    if factor isa FactorRef
        return factorIndex(graph, factor)
    end

    error("Factor references must be Int, Symbol, String, GaussianFactor, or DiscreteFactor.")
end

function graphFigureFactorIndex(
    graph::AbstractFactorGraph,
    factor::Union{GaussianFactor, DiscreteFactor}
)
    if 1 <= factor.id <= length(graph.factors) && graph.factors[factor.id] == factor
        return factor.id
    end

    if !isempty(factor.label)
        return factorIndex(graph, factor.label)
    end

    matches = Int[]

    for (index, graphFactor) in pairs(graph.factors)
        if sameGraphFigureFactor(graphFactor, factor)
            push!(matches, index)
        end
    end

    if length(matches) == 1
        return only(matches)
    elseif isempty(matches)
        error("Factor node is not defined in the graph.")
    end

    error("Factor node matches multiple graph factors; use a factor index or label.")
end

function sameGraphFigureFactor(graphFactor::GaussianFactor, factor::GaussianFactor)
    return graphFactor.variables == factor.variables &&
        graphFactor.mean == factor.mean &&
        graphFactor.coefficient == factor.coefficient &&
        graphFactor.covariance == factor.covariance &&
        graphFactor.initialize == factor.initialize &&
        (isempty(factor.label) || graphFactor.label == factor.label)
end

function sameGraphFigureFactor(graphFactor::DiscreteFactor, factor::DiscreteFactor)
    return graphFactor.variables == factor.variables &&
        graphFactor.table == factor.table &&
        graphFactor.initialize == factor.initialize &&
        (isempty(factor.label) || graphFactor.label == factor.label)
end

sameGraphFigureFactor(_, _) = false

function graphFigureEdgeIndex(graph::AbstractFactorGraph, variable, factor)
    variableIndexValue = graphFigureVariableIndex(graph, variable)
    factorIndexValue = graphFigureFactorIndex(graph, factor)

    for edge in graph.edges
        if edge.variableIndex == variableIndexValue && edge.factorIndex == factorIndexValue
            return edge.id
        end
    end

    error("Variable and factor are not connected by an edge.")
end

function graphFigureEdgeIndex(graph::AbstractFactorGraph, edge::Int)
    if edge < 1 || edge > length(graph.edges)
        error("Edge reference $edge is not defined.")
    end

    return edge
end

function graphFigureEdgeIndex(graph::AbstractFactorGraph, edge::Edge)
    if 1 <= edge.id <= length(graph.edges) && graph.edges[edge.id] == edge
        return edge.id
    end

    return graphFigureEdgeIndex(graph, edge.variableIndex, edge.factorIndex)
end

function graphFigureEdgeIndex(graph::AbstractFactorGraph, selector::Tuple)
    return graphFigureEdgeIndex(graph, selector[1], selector[2])
end
