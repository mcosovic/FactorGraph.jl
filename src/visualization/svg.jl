function escapeXML(text::AbstractString)
    escaped = replace(text, '&' => "&amp;")
    escaped = replace(escaped, '<' => "&lt;")
    escaped = replace(escaped, '>' => "&gt;")
    escaped = replace(escaped, '"' => "&quot;")
    escaped = replace(escaped, '\'' => "&apos;")

    return escaped
end

function svgHeader(width::Int, height::Int, zoom::Real = 1.0)
    displayWidth = ceil(Int, width * zoom)
    displayHeight = ceil(Int, height * zoom)

    return string(
        """<svg xmlns="http://www.w3.org/2000/svg" """,
        """width="$displayWidth" height="$displayHeight" viewBox="0 0 $width $height">
  <title>Factor graph</title>"""
    )
end

function svgStyle(
    fontSize::Int,
    style::NamedTuple
)
    return """
  <style>
    .background { fill: $(style.backgroundFill); }
    .edge {
      fill: none;
      stroke: $(style.edgeStroke);
      stroke-width: $(style.edgeStrokeWidth);
      stroke-opacity: $(style.edgeOpacity);
    }
    .edge-hitbox {
      fill: none;
      stroke: #000000;
      stroke-width: 14;
      stroke-opacity: 0;
      pointer-events: stroke;
    }
    .variable {
      fill: $(style.variableFill);
      stroke: $(style.variableStroke);
      stroke-width: $(style.variableStrokeWidth);
    }
    .factor {
      fill: $(style.factorFill);
      stroke: $(style.factorStroke);
      stroke-width: $(style.factorStrokeWidth);
    }
    .edge-highlight { stroke-opacity: 1.0; }
    .edge-label {
      fill: $(style.edgeStroke);
      font-size: $(ceil(Int, 0.85 * fontSize))px;
      paint-order: stroke;
      stroke: $(style.backgroundFill);
      stroke-linejoin: round;
      stroke-width: 3px;
      text-anchor: middle;
      dominant-baseline: middle;
    }
    .label {
      fill: $(style.labelFill);
      font-family: Helvetica, Arial, sans-serif;
      font-size: $(fontSize)px;
      text-anchor: middle;
      dominant-baseline: middle;
    }
    .label-start { text-anchor: start; }
    .label-end { text-anchor: end; }
  </style>"""
end

function startGraphFigureBuffer(
    canvasWidth::Int,
    canvasHeight::Int,
    zoom::Real,
    fontSize::Int,
    style::NamedTuple
)
    buffer = IOBuffer()

    println(buffer, svgHeader(canvasWidth, canvasHeight, zoom))
    println(
        buffer,
        svgStyle(
            fontSize,
            style
        )
    )
    println(
        buffer,
        """  <rect class="background" x="0" y="0" width="$canvasWidth" height="$canvasHeight"/>"""
    )

    return buffer
end

function svgLine(
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    className::AbstractString,
    style::Union{Nothing, NamedTuple} = nothing,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    attributes = string(
        "class=\"$className\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\"",
        svgStyleAttribute(style)
    )

    return svgShape("line", attributes, tooltip)
end

function svgPath(
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    className::AbstractString,
    style::Union{Nothing, NamedTuple} = nothing,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    dx = abs(x2 - x1)
    direction = x2 >= x1 ? 1 : -1
    c1x = x1 + direction * 0.45 * dx
    c2x = x2 - direction * 0.45 * dx

    attributes = string(
        "class=\"$className\" d=\"M $x1 $y1 C $c1x $y1, $c2x $y2, $x2 $y2\"",
        svgStyleAttribute(style)
    )

    return svgShape("path", attributes, tooltip)
end

function svgEdge(
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    curvedEdges::Bool,
    className::AbstractString,
    style::Union{Nothing, NamedTuple} = nothing,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    visibleEdge = curvedEdges ?
        svgPath(x1, y1, x2, y2, className, style) :
        svgLine(x1, y1, x2, y2, className, style)

    if isnothing(tooltip)
        return visibleEdge
    end

    hitbox = curvedEdges ?
        svgPath(x1, y1, x2, y2, "edge-hitbox", nothing, tooltip) :
        svgLine(x1, y1, x2, y2, "edge-hitbox", nothing, tooltip)

    return string(visibleEdge, "\n", hitbox)
end

function svgEdgeIdLabel(edgeId::Int, x1::Real, y1::Real, x2::Real, y2::Real)
    x = (x1 + x2) / 2
    y = (y1 + y2) / 2

    return """  <text class="edge-label" x="$x" y="$y">e$edgeId</text>"""
end

function svgStyleAttribute(style::Union{Nothing, NamedTuple})
    if isnothing(style)
        return ""
    end

    parts = String[]

    if !isnothing(style.stroke)
        push!(parts, "stroke: $(style.stroke)")
    end

    if !isnothing(style.fill)
        push!(parts, "fill: $(style.fill)")
    end

    if !isnothing(style.strokeWidth)
        push!(parts, "stroke-width: $(style.strokeWidth)")
    end

    return isempty(parts) ? "" : " style=\"$(join(parts, "; "))\""
end

function svgShape(
    name::AbstractString,
    attributes::AbstractString,
    tooltip::Union{Nothing, AbstractString}
)
    if isnothing(tooltip)
        return """  <$name $attributes/>"""
    end

    return """
  <$name $attributes>
    <title>$tooltip</title>
  </$name>"""
end

function svgCircle(
    x::Real,
    y::Real,
    radius::Int,
    className::AbstractString,
    style::Union{Nothing, NamedTuple} = nothing,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    attributes = string(
        "class=\"$className\" cx=\"$x\" cy=\"$y\" r=\"$radius\"",
        svgStyleAttribute(style)
    )

    return svgShape("circle", attributes, tooltip)
end

function svgRect(
    x::Real,
    y::Real,
    width::Int,
    height::Int,
    labelPlacement::Symbol,
    className::AbstractString,
    style::Union{Nothing, NamedTuple} = nothing,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    radius = labelPlacement == :inside ? 4 : 3

    attributes = string(
        "class=\"$className\" x=\"$x\" y=\"$y\" width=\"$width\" height=\"$height\" rx=\"$radius\"",
        svgStyleAttribute(style)
    )

    return svgShape("rect", attributes, tooltip)
end

function svgText(x::Real, y::Real, label::AbstractString)
    return """  <text class="label" x="$x" y="$y">$label</text>"""
end

function svgTextRotated(x::Real, y::Real, label::AbstractString)
    return """  <text class="label" x="$x" y="$y" transform="rotate(-90 $x $y)">$label</text>"""
end

function drawVariableLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    radius::Int,
    label::AbstractString,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    fontSize::Int
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label))
    else
        println(buffer, svgText(x, y + radius + outsideLabelGap + fontSize / 2, label))
    end

    return nothing
end

function drawFactorLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    width::Int,
    variableX::Real,
    label::AbstractString,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    fontSize::Int
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label))
    elseif x < variableX
        labelX = x - width / 2 - outsideLabelGap
        println(buffer, svgTextEnd(labelX, y, label))
    else
        labelX = x + width / 2 + outsideLabelGap
        println(buffer, svgTextStart(labelX, y, label))
    end

    return nothing
end

function drawRotatedFactorLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    height::Int,
    variableY::Real,
    label::AbstractString,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    fontSize::Int
)
    if labelPlacement == :inside
        println(buffer, svgTextRotated(x, y, label))
    elseif y < variableY
        labelY = y - height / 2 - outsideLabelGap - fontSize / 2
        println(buffer, svgTextRotated(x, labelY, label))
    else
        labelY = y + height / 2 + outsideLabelGap + fontSize / 2
        println(buffer, svgTextRotated(x, labelY, label))
    end

    return nothing
end

function drawAdaptiveVerticalFactorLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    height::Int,
    variableY::Real,
    label::AbstractString,
    rotateLabel::Bool,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    fontSize::Int
)
    if rotateLabel
        drawRotatedFactorLabel!(
            buffer,
            x,
            y,
            height,
            variableY,
            label,
            labelPlacement,
            outsideLabelGap,
            fontSize
        )

        return nothing
    end

    if labelPlacement == :inside
        println(buffer, svgText(x, y, label))
    elseif y < variableY
        labelY = y - height / 2 - outsideLabelGap - fontSize / 2
        println(buffer, svgText(x, labelY, label))
    else
        labelY = y + height / 2 + outsideLabelGap + fontSize / 2
        println(buffer, svgText(x, labelY, label))
    end

    return nothing
end

function drawTreeFactorLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    height::Int,
    label::AbstractString,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    fontSize::Int
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label))
    else
        labelY = y - height / 2 - outsideLabelGap - fontSize / 2
        println(buffer, svgText(x, labelY, label))
    end

    return nothing
end

function drawSideLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    label::AbstractString,
    labelPlacement::Symbol,
    side::Symbol
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label))
    elseif side == :right
        println(buffer, svgTextStart(x, y, label))
    elseif side == :left
        println(buffer, svgTextEnd(x, y, label))
    else
        error("Label side must be :left or :right.")
    end

    return nothing
end

function svgTextStart(x::Real, y::Real, label::AbstractString)
    return """  <text class="label label-start" x="$x" y="$y">$label</text>"""
end

function svgTextEnd(x::Real, y::Real, label::AbstractString)
    return """  <text class="label label-end" x="$x" y="$y">$label</text>"""
end
