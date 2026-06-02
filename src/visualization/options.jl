function graphFigureCanvasOptions(canvas::NamedTuple)
    defaultCanvas = (
        width = nothing,
        height = nothing,
        padding = 16,
        zoom = 1.0
    )

    for key in propertynames(canvas)
        if !(key in propertynames(defaultCanvas))
            error("Unknown graph figure canvas key $key.")
        end
    end

    options = merge(defaultCanvas, canvas)

    if options.zoom <= 0
        error("canvas.zoom must be positive.")
    end

    if options.padding < 0
        error("canvas.padding must be non-negative.")
    end

    return options
end

function graphFigureLayoutOptions(layout::NamedTuple)
    defaultLayout = (
        orientation = :horizontal,
        rowSpacing = 80,
        columnSpacing = 160,
        unaryFactorOffset = 90,
        multiFactorOffset = 210,
        curvedEdges = true
    )

    for key in propertynames(layout)
        if !(key in propertynames(defaultLayout))
            error("Unknown graph figure layout key $key.")
        end
    end

    options = merge(defaultLayout, layout)

    if !(options.orientation in (:horizontal, :vertical))
        error("layout.orientation must be :horizontal or :vertical.")
    end

    if options.unaryFactorOffset <= 0 || options.multiFactorOffset <= 0
        error("layout unary and multi-factor offsets must be positive.")
    end

    if options.rowSpacing <= 0 || options.columnSpacing <= 0
        error("layout rowSpacing and columnSpacing must be positive.")
    end

    return options
end

function graphFigureNodeOptions(node::NamedTuple)
    defaultNode = (
        variableRadius = 24,
        factorSize = 22
    )

    for key in propertynames(node)
        if !(key in propertynames(defaultNode))
            error("Unknown graph figure node key $key.")
        end
    end

    options = merge(defaultNode, node)

    if options.variableRadius <= 0 || options.factorSize <= 0
        error("node variableRadius and factorSize must be positive.")
    end

    return options
end

function graphFigureLabelOptions(label::NamedTuple)
    defaultLabel = (
        placement = :outside,
        outsideGap = 6,
        showVariables = true,
        showFactors = true,
        showTooltips = true,
        showEdgeIds = false,
        tooltipDetail = :summary,
        fontSize = 14
    )

    for key in propertynames(label)
        if !(key in propertynames(defaultLabel))
            error("Unknown graph figure label key $key.")
        end
    end

    options = merge(defaultLabel, label)

    if !(options.placement in (:inside, :outside))
        error("label.placement must be :inside or :outside.")
    end

    if options.outsideGap < 0
        error("label.outsideGap must be non-negative.")
    end

    if options.fontSize <= 0
        error("label.fontSize must be positive.")
    end

    if !(options.tooltipDetail in (:summary, :full))
        error("label.tooltipDetail must be :summary or :full.")
    end

    return options
end


function graphFigureStyle(style::NamedTuple)
    defaultStyle = (
        backgroundFill = "#ffffff",
        variableFill = "#e0f2fe",
        variableStroke = "#0369a1",
        variableStrokeWidth = 1.8,
        factorFill = "#dc2626",
        factorStroke = "#991b1b",
        factorStrokeWidth = 1.8,
        edgeStroke = "#64748b",
        edgeStrokeWidth = 1.6,
        edgeOpacity = 1.0,
        labelFill = "#111827"
    )

    for key in propertynames(style)
        if !(key in propertynames(defaultStyle))
            error("Unknown graph figure style key $key.")
        end
    end

    graphStyle = merge(defaultStyle, style)

    if graphStyle.variableStrokeWidth <= 0 ||
            graphStyle.factorStrokeWidth <= 0 ||
            graphStyle.edgeStrokeWidth <= 0
        error("Graph figure style stroke widths must be positive.")
    end

    return graphStyle
end
