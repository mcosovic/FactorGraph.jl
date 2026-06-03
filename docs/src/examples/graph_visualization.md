# Graph Visualization

This example shows common ways to customize [`graphFigure`](@ref) and
[`saveGraphFigure`](@ref). The same visualization API works for Gaussian
factor graphs, discrete factor graphs, and tree views.

The examples below intentionally use a small graph so the effect of each option
is easy to see.

## Base Graph

Start with a small Gaussian factor graph:

```@example graph_visualization
using FactorGraph

variables = [
    GaussianVariable(:x1, 1; label = "x1"),
    GaussianVariable(:x2, 1; label = "x2"),
    GaussianVariable(:x3, 1; label = "x3")
]

factors = [
    GaussianFactor(:x1, 0.0, 1.0, 0.5; label = "f1", initialize = true),
    GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.3; label = "f2"),
    GaussianFactor(:x2, :x3, 0.0, [1.0 -1.0], 0.3; label = "f3"),
    GaussianFactor(:x3, 2.0, 1.0, 0.8; label = "f4")
]

graph = factorGraph(variables, factors)

nothing # hide
```

The graph has three variable nodes and four factor nodes. The factor labels are
kept short (`f1`, `f2`, ...) so label placement and edge id labels remain easy
to inspect in the generated SVGs.

---

## Default Figure

The default figure draws variable labels, factor labels, SVG hover tooltips,
and curved edges. This is usually enough for a quick model check:

```@example graph_visualization
saveGraphFigure("../../gv_default.svg", graph)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../gv_default.svg"
    type="image/svg+xml"
    aria-label="Default graph figure"
    style="width: 42%; height: auto;">
    <a href="../../gv_default.svg">Default graph figure</a>
  </object>
</div>
```

Hover over a node or edge in the rendered SVG to inspect its tooltip. Tooltips
show identity metadata by default and can be expanded with `tooltipDetail =
:full`.

---

## Layout And Labels

Layout options control orientation, spacing, and edge geometry. Label options
control node labels, edge ids, tooltip detail, and font size. This example
switches to a vertical layout, draws straight edges, and enables edge id labels:

```@example graph_visualization
saveGraphFigure(
    "../../gv_vertical.svg",
    graph;
    layout = (
        orientation = :vertical,
        rowSpacing = (90, 150),
        columnSpacing = 120,
        curvedEdges = false
    ),
    label = (
        placement = :outside,
        outsideGap = 8,
        showEdgeIds = true,
        tooltipDetail = :full,
        fontSize = 13
    )
)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../gv_vertical.svg"
    type="image/svg+xml"
    aria-label="Vertical graph figure with edge ids"
    style="width: 34%; height: auto;">
    <a href="../../gv_vertical.svg">Vertical graph figure with edge ids</a>
  </object>
</div>
```

Spacing options may be scalars for uniform gaps, or tuples/vectors for per-gap
spacing. In the default horizontal layout, `columnSpacing = (90, 210)` separates
unary factors from variables, then variables from multi-variable factors. In
vertical layout, the same idea applies through `rowSpacing`.

---

## Focused Views

Use `view.variables` and `view.factors` to choose focus nodes. The `hops`
keyword expands through the bipartite graph by edge count. With `hops = 2`,
the figure includes `x1`, its neighboring factors, and the variables connected
to those factors:

```@example graph_visualization
saveGraphFigure(
    "../../gv_focus.svg",
    graph;
    view = (variables = [:x1], hops = 2),
    label = (showEdgeIds = true,)
)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../gv_focus.svg"
    type="image/svg+xml"
    aria-label="Focused graph figure"
    style="width: 36%; height: auto;">
    <a href="../../gv_focus.svg">Focused graph figure</a>
  </object>
</div>
```

`hops = 0` draws only the selected seed nodes, while `hops = :all` expands
through the connected component. Focused nodes are drawn normally and expanded
context nodes are drawn with the context style.

---

## Style And Highlights

Style options set the default colors and stroke widths. Highlight entries can
select variables, factors, edges, or a variable-factor edge pair. Highlights
are useful when a figure is meant to explain one local relationship in a larger
model:

```@example graph_visualization
saveGraphFigure(
    "../../gv_highlight.svg",
    graph;
    style = (
        backgroundFill = "#f8fafc",
        variableFill = "#e0f2fe",
        variableStroke = "#0369a1",
        variableStrokeWidth = 2.2,
        factorFill = "#fee2e2",
        factorStroke = "#b91c1c",
        factorStrokeWidth = 2.2,
        edgeStroke = "#64748b",
        edgeStrokeWidth = 1.8,
        edgeOpacity = 0.9,
        labelFill = "#0f172a"
    ),
    highlight = [
        (variable = :x2, stroke = "#16a34a", fill = "#dcfce7", strokeWidth = 4),
        (factor = "f2", stroke = "#f59e0b", fill = "#fef3c7", strokeWidth = 4),
        (variable = :x3, factor = "f4", stroke = "#7c3aed", strokeWidth = 4)
    ]
)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../gv_highlight.svg"
    type="image/svg+xml"
    aria-label="Highlighted graph figure"
    style="width: 42%; height: auto;">
    <a href="../../gv_highlight.svg">Highlighted graph figure</a>
  </object>
</div>
```

Variable and factor highlights can also highlight incident edges. Edge-specific
highlights can be selected either by edge id or by a `(variable, factor)` pair.

---

## Tree Views

Tree figures use graph depth for placement. The same label, style, highlight,
and view options are available. Here the graph is converted to a tree view and
drawn horizontally from root variable `x1`:

```@example graph_visualization
tree = treeFactorGraph(graph; root = :x1)

saveGraphFigure(
    "../../gv_tree.svg",
    tree;
    layout = (orientation = :horizontal, rowSpacing = 70, columnSpacing = (95, 120, 145)),
    view = (variables = [:x1], hops = :all),
    label = (showEdgeIds = true, tooltipDetail = :full),
    highlight = [(variable = :x1, stroke = "#16a34a", fill = "#dcfce7", strokeWidth = 4)]
)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="../../gv_tree.svg"
    type="image/svg+xml"
    aria-label="Tree graph figure"
    style="width: 90%; height: auto;">
    <a href="../../gv_tree.svg">Tree graph figure</a>
  </object>
</div>
```

In horizontal tree layout, `columnSpacing` separates depth levels and
`rowSpacing` separates nodes within the same level. The interpretation is
reversed for vertical tree orientation. Tuple or vector spacing sets level gaps
one by one; the last value is reused when the graph has more gaps.

---

## Full Option Sketch

This compact call shows every option group in one place. It is meant as a
reference pattern rather than a recommended visual style:

```julia
graphFigure(
    graph;
    canvas = (width = 700, height = nothing, padding = 24, zoom = 1.0),
    layout = (
        orientation = :horizontal,
        rowSpacing = 90,
        columnSpacing = (90, 210),
        curvedEdges = true
    ),
    node = (variableRadius = 24, factorSize = 22),
    label = (
        placement = :outside,
        outsideGap = 6,
        showVariables = true,
        showFactors = true,
        showTooltips = true,
        showEdgeIds = false,
        tooltipDetail = :summary,
        fontSize = 14
    ),
    view = (variables = [:x1], factors = nothing, hops = :all),
    style = (
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
    ),
    highlight = [
        (variable = :x1, stroke = "#16a34a", fill = "#dcfce7", strokeWidth = 4),
        (factor = "f2", stroke = "#f59e0b", fill = "#fef3c7", strokeWidth = 4),
        (edge = 1, stroke = "#7c3aed", strokeWidth = 4)
    ]
)
```
