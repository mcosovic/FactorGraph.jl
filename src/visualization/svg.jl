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
    .variable-context { opacity: 0.38; }
    .factor {
      fill: $(style.factorFill);
      stroke: $(style.factorStroke);
      stroke-width: $(style.factorStrokeWidth);
    }
    .factor-context { opacity: 0.42; }
    .edge-highlight { stroke-opacity: 1.0; }
    .edge-label {
      fill: $(style.edgeStroke);
      font-size: $(ceil(Int, 0.85 * fontSize))px;
      text-anchor: middle;
      dominant-baseline: middle;
    }
    .edge-label-background {
      fill: $(style.backgroundFill);
      stroke: $(style.backgroundFill);
      stroke-linejoin: round;
      stroke-width: 4px;
      font-size: $(ceil(Int, 0.85 * fontSize))px;
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
    tooltip::Union{Nothing, AbstractString} = nothing,
    orientation::Symbol = :horizontal
)
    if orientation == :horizontal
        dx = abs(x2 - x1)
        direction = x2 >= x1 ? 1 : -1
        c1x = x1 + direction * 0.45 * dx
        c2x = x2 - direction * 0.45 * dx
        path = "M $x1 $y1 C $c1x $y1, $c2x $y2, $x2 $y2"
    elseif orientation == :vertical
        dy = abs(y2 - y1)
        direction = y2 >= y1 ? 1 : -1
        c1y = y1 + direction * 0.45 * dy
        c2y = y2 - direction * 0.45 * dy
        path = "M $x1 $y1 C $x1 $c1y, $x2 $c2y, $x2 $y2"
    else
        error("orientation must be :horizontal or :vertical.")
    end

    attributes = string(
        "class=\"$className\" d=\"$path\"",
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
    tooltip::Union{Nothing, AbstractString} = nothing,
    orientation::Symbol = :horizontal
)
    visibleEdge = curvedEdges ?
        svgPath(x1, y1, x2, y2, className, style, nothing, orientation) :
        svgLine(x1, y1, x2, y2, className, style)

    if isnothing(tooltip)
        return visibleEdge
    end

    hitbox = curvedEdges ?
        svgPath(x1, y1, x2, y2, "edge-hitbox", nothing, tooltip, orientation) :
        svgLine(x1, y1, x2, y2, "edge-hitbox", nothing, tooltip)

    return string(visibleEdge, "\n", hitbox)
end

function svgEdgeIdLabelParts(
    edgeId::Int,
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    occupiedLabels::AbstractVector,
    curvedEdges::Bool,
    orientation::Symbol,
    edgeGeometries::AbstractVector,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    x, y, width, height = edgeIdLabelPosition!(
        occupiedLabels,
        edgeId,
        x1,
        y1,
        x2,
        y2,
        curvedEdges,
        orientation,
        edgeGeometries
    )

    return svgEdgeIdLabelParts(edgeId, x, y, width, height, tooltip)
end

function svgEdgeIdLabelParts(
    edgeId::Int,
    x::Real,
    y::Real,
    width::Int,
    height::Int,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    label = "e$edgeId"
    text = isnothing(tooltip) ?
        """  <text class="edge-label" x="$x" y="$y">$label</text>""" :
        """  <text class="edge-label" x="$x" y="$y"><title>$tooltip</title>$label</text>"""

    return (
        background = """  <text class="edge-label-background" x="$x" y="$y">$label</text>""",
        text = text
    )
end

function drawEdgeIdLabels!(
    buffer::IO,
    backgrounds::AbstractVector{<:AbstractString},
    texts::AbstractVector{<:AbstractString}
)
    for background in backgrounds
        println(buffer, background)
    end

    for text in texts
        println(buffer, text)
    end

    return nothing
end

function edgeIdLabelPosition!(
    occupiedLabels::AbstractVector,
    edgeId::Int,
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    curvedEdges::Bool,
    orientation::Symbol,
    edgeGeometries::AbstractVector
)
    label = "e$edgeId"
    fontSize = 12
    width = ceil(Int, approximateTextWidth(label, fontSize) + 5)
    height = fontSize + 4
    baseX, baseY = edgePoint(x1, y1, x2, y2, curvedEdges, orientation, 0.5)
    fallback = nothing
    fallbackScore = (typemax(Int), Inf)

    for t in edgeIdLabelCandidatePositions()
        x, y = edgePoint(x1, y1, x2, y2, curvedEdges, orientation, t)
        overlapsLabel = edgeIdLabelOverlaps(occupiedLabels, x, y, width, height)

        if overlapsLabel
            continue
        end

        if isnothing(fallback)
            fallback = (x = x, y = y, width = width, height = height)
        end

        intersectionCount = edgeIdLabelIntersectionCount(
            edgeGeometries,
            edgeId,
            x,
            y,
            width,
            height
        )

        score = (intersectionCount, abs(t - 0.5))

        if score < fallbackScore
            fallback = (x = x, y = y, width = width, height = height)
            fallbackScore = score
        end

        if intersectionCount == 0
            push!(occupiedLabels, (x = x, y = y, width = width, height = height))
            return x, y, width, height
        end
    end

    if !isnothing(fallback)
        push!(occupiedLabels, fallback)
        return fallback.x, fallback.y, fallback.width, fallback.height
    end

    push!(occupiedLabels, (x = baseX, y = baseY, width = width, height = height))
    return baseX, baseY, width, height
end

function edgeGeometry(
    edgeId::Int,
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    curvedEdges::Bool,
    orientation::Symbol
)
    return (
        id = edgeId,
        x1 = x1,
        y1 = y1,
        x2 = x2,
        y2 = y2,
        curvedEdges = curvedEdges,
        orientation = orientation
    )
end

function edgeIdLabelCandidatePositions()
    positions = [0.5]

    for offset in 0.02:0.02:0.40
        push!(positions, 0.5 - offset)
        push!(positions, 0.5 + offset)
    end

    return positions
end

function edgePoint(
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    curvedEdges::Bool,
    orientation::Symbol,
    t::Real
)
    if !curvedEdges
        return x1 + (x2 - x1) * t, y1 + (y2 - y1) * t
    end

    if orientation == :horizontal
        dx = abs(x2 - x1)
        direction = x2 >= x1 ? 1 : -1
        c1x = x1 + direction * 0.45 * dx
        c1y = y1
        c2x = x2 - direction * 0.45 * dx
        c2y = y2
    elseif orientation == :vertical
        dy = abs(y2 - y1)
        direction = y2 >= y1 ? 1 : -1
        c1x = x1
        c1y = y1 + direction * 0.45 * dy
        c2x = x2
        c2y = y2 - direction * 0.45 * dy
    else
        error("orientation must be :horizontal or :vertical.")
    end

    mt = 1 - t
    x = mt^3 * x1 + 3 * mt^2 * t * c1x + 3 * mt * t^2 * c2x + t^3 * x2
    y = mt^3 * y1 + 3 * mt^2 * t * c1y + 3 * mt * t^2 * c2y + t^3 * y2

    return x, y
end

function edgeIdLabelOverlaps(
    occupiedLabels::AbstractVector,
    x::Real,
    y::Real,
    width::Int,
    height::Int
)
    for label in occupiedLabels
        if abs(x - label.x) < (width + label.width) / 2 + 3 &&
                abs(y - label.y) < (height + label.height) / 2 + 3
            return true
        end
    end

    return false
end

function edgeIdLabelIntersectionCount(
    edgeGeometries::AbstractVector,
    edgeId::Int,
    x::Real,
    y::Real,
    width::Int,
    height::Int
)
    halfWidth = width / 2 + 5
    halfHeight = height / 2 + 5
    count = 0

    for edge in edgeGeometries
        if edge.id == edgeId
            continue
        end

        previousPoint = nothing

        for t in 0.04:0.02:0.96
            edgeX, edgeY = edgePoint(
                edge.x1,
                edge.y1,
                edge.x2,
                edge.y2,
                edge.curvedEdges,
                edge.orientation,
                t
            )

            if abs(edgeX - x) <= halfWidth && abs(edgeY - y) <= halfHeight
                count += 1
            end

            if !isnothing(previousPoint) &&
                    segmentIntersectsRectangle(
                        previousPoint.x,
                        previousPoint.y,
                        edgeX,
                        edgeY,
                        x - halfWidth,
                        y - halfHeight,
                        x + halfWidth,
                        y + halfHeight
                    )
                count += 1
            end

            previousPoint = (x = edgeX, y = edgeY)
        end
    end

    return count
end

function segmentIntersectsRectangle(
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real,
    left::Real,
    top::Real,
    right::Real,
    bottom::Real
)
    if pointInRectangle(x1, y1, left, top, right, bottom) ||
            pointInRectangle(x2, y2, left, top, right, bottom)
        return true
    end

    return segmentsIntersect(x1, y1, x2, y2, left, top, right, top) ||
        segmentsIntersect(x1, y1, x2, y2, right, top, right, bottom) ||
        segmentsIntersect(x1, y1, x2, y2, right, bottom, left, bottom) ||
        segmentsIntersect(x1, y1, x2, y2, left, bottom, left, top)
end

function pointInRectangle(
    x::Real,
    y::Real,
    left::Real,
    top::Real,
    right::Real,
    bottom::Real
)
    return left <= x <= right && top <= y <= bottom
end

function segmentsIntersect(
    ax::Real,
    ay::Real,
    bx::Real,
    by::Real,
    cx::Real,
    cy::Real,
    dx::Real,
    dy::Real
)
    denominator = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx)

    if denominator == 0
        return false
    end

    u = ((cx - ax) * (dy - cy) - (cy - ay) * (dx - cx)) / denominator
    v = ((cx - ax) * (by - ay) - (cy - ay) * (bx - ax)) / denominator

    return 0 <= u <= 1 && 0 <= v <= 1
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

function svgTextElement(
    className::AbstractString,
    x::Real,
    y::Real,
    label::AbstractString,
    transform::AbstractString = "";
    tooltip::Union{Nothing, AbstractString} = nothing
)
    attributes = "class=\"$className\" x=\"$x\" y=\"$y\""
    if !isempty(transform)
        attributes = "$attributes transform=\"$transform\""
    end

    if isnothing(tooltip)
        return """  <text $attributes>$label</text>"""
    end

    return """  <text $attributes><title>$tooltip</title>$label</text>"""
end

function svgText(
    x::Real,
    y::Real,
    label::AbstractString,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    return svgTextElement("label", x, y, label; tooltip = tooltip)
end

function svgTextRotated(
    x::Real,
    y::Real,
    label::AbstractString,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    return svgTextElement("label", x, y, label, "rotate(-90 $x $y)"; tooltip = tooltip)
end

function svgTextRotatedStart(
    x::Real,
    y::Real,
    label::AbstractString,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    return svgTextElement("label label-start", x, y, label, "rotate(-90 $x $y)"; tooltip = tooltip)
end

function svgTextRotatedEnd(
    x::Real,
    y::Real,
    label::AbstractString,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    return svgTextElement("label label-end", x, y, label, "rotate(-90 $x $y)"; tooltip = tooltip)
end

function drawVariableLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    radius::Int,
    label::AbstractString,
    labelPlacement::Symbol,
    outsideLabelGap::Int,
    fontSize::Int,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label, tooltip))
    else
        println(buffer, svgText(x, y + radius + outsideLabelGap + fontSize / 2, label, tooltip))
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
    fontSize::Int,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label, tooltip))
    elseif x < variableX
        labelX = x - width / 2 - outsideLabelGap
        println(buffer, svgTextEnd(labelX, y, label, tooltip))
    else
        labelX = x + width / 2 + outsideLabelGap
        println(buffer, svgTextStart(labelX, y, label, tooltip))
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
    fontSize::Int,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    if labelPlacement == :inside
        println(buffer, svgTextRotated(x, y, label, tooltip))
    elseif y < variableY
        labelY = y - height / 2 - outsideLabelGap
        println(buffer, svgTextRotatedStart(x, labelY, label, tooltip))
    else
        labelY = y + height / 2 + outsideLabelGap
        println(buffer, svgTextRotatedEnd(x, labelY, label, tooltip))
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
    fontSize::Int,
    tooltip::Union{Nothing, AbstractString} = nothing
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
            fontSize,
            tooltip
        )

        return nothing
    end

    if labelPlacement == :inside
        println(buffer, svgText(x, y, label, tooltip))
    elseif y < variableY
        labelY = y - height / 2 - outsideLabelGap - fontSize / 2
        println(buffer, svgText(x, labelY, label, tooltip))
    else
        labelY = y + height / 2 + outsideLabelGap + fontSize / 2
        println(buffer, svgText(x, labelY, label, tooltip))
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
    fontSize::Int,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label, tooltip))
    else
        labelY = y - height / 2 - outsideLabelGap - fontSize / 2
        println(buffer, svgText(x, labelY, label, tooltip))
    end

    return nothing
end

function drawSideLabel!(
    buffer::IO,
    x::Real,
    y::Real,
    label::AbstractString,
    labelPlacement::Symbol,
    side::Symbol,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    if labelPlacement == :inside
        println(buffer, svgText(x, y, label, tooltip))
    elseif side == :right
        println(buffer, svgTextStart(x, y, label, tooltip))
    elseif side == :left
        println(buffer, svgTextEnd(x, y, label, tooltip))
    else
        error("Label side must be :left or :right.")
    end

    return nothing
end

function svgTextStart(
    x::Real,
    y::Real,
    label::AbstractString,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    return svgTextElement("label label-start", x, y, label; tooltip = tooltip)
end

function svgTextEnd(
    x::Real,
    y::Real,
    label::AbstractString,
    tooltip::Union{Nothing, AbstractString} = nothing
)
    return svgTextElement("label label-end", x, y, label; tooltip = tooltip)
end
