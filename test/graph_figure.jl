include("setup.jl")

function graphFigureTestGraph()
    variables = [
        DiscreteVariable(:x1, 2; label = "x1"),
        DiscreteVariable(:x2, 2; label = "x2"),
        DiscreteVariable(:x3, 2; label = "x3")
    ]

    factors = [
        DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1"),
        DiscreteFactor(:x2, [0.5, 0.5]; label = "prior_x2"),
        DiscreteFactor(:x1, :x2, [0.9 0.1; 0.2 0.8]; label = "link_1_2"),
        DiscreteFactor(:x2, :x3, [0.8 0.2; 0.3 0.7]; label = "link_2_3")
    ]

    return factorGraph(variables, factors)
end

function graphFigureLineLengths(svg::AbstractString)
    pattern = Regex(
        "<line class=\"edge(?: edge-highlight)?\" " *
        "x1=\"([^\"]+)\" y1=\"([^\"]+)\" x2=\"([^\"]+)\" y2=\"([^\"]+)\""
    )

    return [
        hypot(
            parse(Float64, match.captures[3]) - parse(Float64, match.captures[1]),
            parse(Float64, match.captures[4]) - parse(Float64, match.captures[2])
        )
        for match in eachmatch(pattern, svg)
    ]
end

function graphFigureGaussianTooltipGraph()
    variables = [
        GaussianVariable(:x1, 1; label = "x1", mean = [0.0], covariance = [1.0;;]),
        GaussianVariable(
            :x2,
            2;
            label = "x2",
            components = [:a, :b],
            mean = [2.0, -1.0],
            covariance = [4.0 0.5; 0.5 2.0]
        )
    ]
    factors = [
        GaussianFactor(:x1, 0.0, 1.0, 0.5; label = "prior_x1", initialize = true),
        GaussianFactor(:x1, :x2, [0.0], [1.0 -1.0 0.5], [1.2;;]; label = "link")
    ]

    return factorGraph(variables, factors)
end

@testset verbose = true "Graph figures" begin
    graph = graphFigureTestGraph()

    @testset "SVG string output" begin
        svg = graphFigure(
            graph;
            label = (placement = :inside,),
            layout = (curvedEdges = false,)
        )

        @test startswith(svg, "<svg")
        @test occursin("<circle", svg)
        @test occursin("<rect", svg)
        @test occursin("<line", svg)
        @test occursin("Factor graph", svg)

        lengths = graphFigureLineLengths(svg)
        @test lengths[1] < lengths[2]
    end

    @testset "Hides factor labels with compact nodes" begin
        svg = graphFigure(
            graph;
            label = (placement = :inside, showFactors = false)
        )

        @test !occursin(">prior_x1<", svg)
        @test occursin("width=\"18\" height=\"18\"", svg)
    end

    @testset "Vertical orientation" begin
        svg = graphFigure(
            graph;
            layout = (orientation = :vertical, curvedEdges = false),
            label = (placement = :inside,)
        )

        @test startswith(svg, "<svg")
        @test occursin(">x1<", svg)
        @test occursin(">link_1_2<", svg)

        lengths = graphFigureLineLengths(svg)
        @test lengths[1] < lengths[2]
    end

    @testset "Per-item highlight styles" begin
        svg = graphFigure(
            graph;
            highlight = [
                (variable = :x2, stroke = "#16a34a", fill = "#dcfce7", strokeWidth = 5),
                (factor = "link_1_2", stroke = "#2563eb", strokeWidth = 4),
                (variable = :x1, factor = "link_1_2", stroke = "#7c3aed", strokeWidth = 6)
            ]
        )

        @test occursin("stroke: #16a34a", svg)
        @test occursin("fill: #dcfce7", svg)
        @test occursin("stroke: #2563eb", svg)
        @test occursin("stroke: #7c3aed", svg)
        @test occursin("stroke-width: 6", svg)
        @test count("stroke-width: 5", svg) > 1
    end

    @testset "Custom graph style" begin
        svg = graphFigure(
            graph;
            style = (
                backgroundFill = "#f8fafc",
                variableFill = "#eef2ff",
                variableStroke = "#4338ca",
                factorFill = "#fecaca",
                edgeStroke = "#475569",
                labelFill = "#0f172a"
            )
        )

        @test occursin(".background { fill: #f8fafc; }", svg)
        @test occursin("fill: #eef2ff", svg)
        @test occursin("stroke: #4338ca", svg)
        @test occursin("fill: #fecaca", svg)
        @test occursin("stroke: #475569", svg)
        @test occursin("fill: #0f172a", svg)
    end

    @testset "Tooltips" begin
        svg = graphFigure(graph)

        @test occursin("edge-hitbox", svg)
        @test occursin("<title>Edge 1", svg)
        @test occursin(r"<text class=\"label[^\"]*\"[^>]*>\s*<title>Variable x1", svg)
        @test occursin(r"<text class=\"label[^\"]*\"[^>]*>\s*<title>Factor prior", svg)
        @test occursin("Identity", svg)
        @test occursin("type: Discrete", svg)
        @test !occursin("type: DiscreteVariable", svg)
        @test !occursin("type: DiscreteFactor", svg)
        @test occursin("Table\n  size: 2 x 2", svg)
        @test !occursin("table: [0.9 0.1; 0.2 0.8]", svg)
        @test !occursin("Factor link_1_2" * "\n\nIdentity" *
            "\n  type: Discrete" *
            "\n  index: 3" *
            "\n  id: 3" *
            "\n  variables: x1, x2" *
            "\n  degree: 2" *
            "\n\nTable" *
            "\n  size: 2 x 2" *
            "\n\nOptions", svg)

        fullSvg = graphFigure(graph; label = (tooltipDetail = :full,))
        @test occursin("Table\n  size: 2\n  axis: x1", fullSvg)
        @test occursin(
            "Table\n  size: 2 x 2\n  rows: x1, columns: x2" *
            "\n    [0.9  0.1\n     0.2  0.8]",
            fullSvg
        )
        @test !occursin("Factor link_1_2" * "\n\nIdentity" *
            "\n  type: Discrete" *
            "\n  index: 3" *
            "\n  id: 3" *
            "\n  variables: x1, x2" *
            "\n  degree: 2" *
            "\n\nTable" *
            "\n  size: 2 x 2" *
            "\n  rows: x1, columns: x2" *
            "\n    [0.9  0.1\n     0.2  0.8]" *
            "\n\nOptions", fullSvg)

        table3d = cat(
            [0.9 0.1 0.2; 0.1 0.9 0.8],
            [0.7 0.4 0.3; 0.3 0.6 0.7];
            dims = 3
        )
        graph3d = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
                DiscreteVariable(:x2, 3; label = "x2", states = [:low, :mid, :high]),
                DiscreteVariable(:x3, 2; label = "x3", states = [:closed, :open])
            ],
            [DiscreteFactor(:x1, :x2, :x3, table3d; label = "table3d")]
        )
        full3dSvg = graphFigure(graph3d; label = (tooltipDetail = :full,))
        @test occursin("Table\n  size: 2 x 3 x 2\n  axes: x1 x x2 x x3", full3dSvg)
        @test occursin("  rows: x1, columns: x2", full3dSvg)
        @test occursin("  slice x3 = :closed", full3dSvg)
        @test occursin("  slice x3 = :open", full3dSvg)
        @test !occursin("slice 3 =", full3dSvg)

        noTooltipSvg = graphFigure(graph; label = (showTooltips = false,))
        @test !occursin("class=\"edge-hitbox\"", noTooltipSvg)
        @test !occursin("<title>Edge", noTooltipSvg)
        @test !occursin("<title>Variable", noTooltipSvg)
        @test !occursin("<title>Factor prior", noTooltipSvg)

        gaussianSvg = graphFigure(graphFigureGaussianTooltipGraph())
        @test occursin("type: Gaussian", gaussianSvg)
        @test !occursin("type: GaussianVariable", gaussianSvg)
        @test !occursin("type: GaussianFactor", gaussianSvg)
        @test occursin("mean size: 2", gaussianSvg)
        @test occursin("covariance size: 2 x 2", gaussianSvg)
        @test occursin("component count: 2", gaussianSvg)

        gaussianFullSvg = graphFigure(
            graphFigureGaussianTooltipGraph();
            label = (tooltipDetail = :full,)
        )
        @test occursin("covariance:\n    [4.0  0.5\n     0.5  2.0]", gaussianFullSvg)
        @test occursin("component: 1", gaussianFullSvg)
        @test occursin("components: [:a, :b]", gaussianFullSvg)
    end

    @testset "Edge ids" begin
        svg = graphFigure(graph; label = (showEdgeIds = true,))

        @test occursin("class=\"edge-label\"", svg)
        @test occursin(r"<text class=\"edge-label\"[^>]*><title>Edge 1", svg)
        @test occursin(r"<text class=\"edge-label\"[^>]*><title>Edge 6", svg)
        @test occursin(">e1</text>", svg)
        @test occursin(">e6</text>", svg)

        hiddenSvg = graphFigure(graph)
        @test !occursin("class=\"edge-label\"", hiddenSvg)

        noTooltipSvg = graphFigure(graph; label = (showEdgeIds = true, showTooltips = false))
        @test occursin("class=\"edge-label\"", noTooltipSvg)
        @test !occursin(r"<text class=\"edge-label\"[^>]*><title>Edge", noTooltipSvg)
    end

    @testset "Focused view" begin
        svg = graphFigure(graph; view = (variables = [:x1, :x2],))

        @test occursin(">x1<", svg)
        @test occursin(">x2<", svg)
        @test !occursin(">x3<", svg)
        @test occursin(">link_1_2<", svg)
        @test occursin(">link_2_3<", svg)

        seedSvg = graphFigure(graph; view = (variables = [:x1], hops = 0))
        @test occursin(">x1<", seedSvg)
        @test !occursin(">link_1_2<", seedSvg)

        twoHopSvg = graphFigure(graph; view = (variables = [:x1, :x2], hops = 2))
        @test occursin(">x3<", twoHopSvg)

        allHopSvg = graphFigure(graph; view = (variables = [:x1], hops = :all))
        @test occursin(">x3<", allHopSvg)

        @test_throws ErrorException graphFigure(
            graph;
            view = (variables = [:x1, :x2], hops = -1)
        )

        verticalSvg = graphFigure(
            graph;
            layout = (orientation = :vertical,),
            view = (variables = [:x1, :x2],)
        )
        @test !occursin(">x3<", verticalSvg)
    end

    @testset "Tree graph output" begin
        tree = treeFactorGraph(graph; root = :x1)
        svg = graphFigure(
            tree;
            layout = (orientation = :vertical,),
            label = (placement = :outside, showEdgeIds = true)
        )

        @test startswith(svg, "<svg")
        @test occursin("label-start", svg)
        @test occursin("edge-hitbox", svg)
        @test occursin("<title>Edge 1", svg)
        @test occursin("Identity", svg)
        @test occursin("Table\n  size: 2 x 2", svg)
        @test occursin("class=\"edge-label\"", svg)
        @test occursin(r"<text class=\"edge-label\"[^>]*><title>Edge 1", svg)
        @test occursin(">e1</text>", svg)

        noTooltipSvg = graphFigure(tree; label = (showTooltips = false,))
        @test !occursin("<title>Edge", noTooltipSvg)
        @test !occursin("class=\"edge-hitbox\"", noTooltipSvg)

        viewSvg = graphFigure(tree; view = (variables = [:x1], hops = 2))
        @test occursin("variable-context", viewSvg)
        @test occursin("<title>Variable x1", viewSvg)
        @test occursin(r"<text class=\"label[^\"]*\"[^>]*>\s*<title>Variable x1", viewSvg)
    end

    @testset "File output" begin
        path = tempname() * ".svg"

        try
            returnedPath = saveGraphFigure(path, graph; label = (placement = :inside,))

            @test returnedPath == path
            @test isfile(path)
            @test startswith(read(path, String), "<svg")
        finally
            rm(path; force = true)
        end
    end

    @testset "Validation" begin
        tree = treeFactorGraph(graph; root = :x1)

        @test_throws ErrorException graphFigure(graph; canvas = (padding = -1,))
        @test_throws ErrorException graphFigure(graph; canvas = (zoom = 0.0,))
        @test_throws ErrorException graphFigure(graph; style = (edgeStrokeWidth = 0,))
        @test_throws ErrorException graphFigure(graph; style = (edgeFill = "none",))
        @test_throws ErrorException graphFigure(graph; layout = (unaryFactorOffset = 0,))
        @test_throws ErrorException graphFigure(graph; layout = (rowSpacing = 0,))
        @test_throws ErrorException graphFigure(graph; layout = (columnSpacing = 0,))
        @test_throws ErrorException graphFigure(graph; node = (variableRadius = 0,))
        @test_throws ErrorException graphFigure(graph; label = (fontSize = 0,))
        @test_throws ErrorException graphFigure(graph; label = (tooltipDetail = :basic,))
        @test_throws ErrorException graphFigure(graph; label = (tooltipDetail = :data,))
        @test_throws ErrorException graphFigure(graph; canvas = (canvasPadding = 16,))
        @test_throws ErrorException graphFigure(graph; layout = (orientation = :diagonal,))
        @test_throws ErrorException graphFigure(graph; label = (labelPlacement = :inside,))
        @test_throws MethodError graphFigure(graph; highlightVariables = [:x1])
        @test_throws MethodError graphFigure(graph; highlightFactors = ["link_1_2"])
        @test_throws MethodError graphFigure(graph; highlightEdges = [(:x1, "link_1_2")])
        @test_throws MethodError graphFigure(graph; edgeOpacity = 0.5)
        @test_throws MethodError graphFigure(graph; highlightStrokeWidth = 0)
        @test_throws MethodError graphFigure(graph; canvasPadding = -1)
        @test_throws MethodError graphFigure(graph; zoom = 0.0)
        @test_throws MethodError graphFigure(graph; unaryFactorOffset = 0)
        @test_throws MethodError graphFigure(graph; rowSpacing = 0)
        @test_throws MethodError graphFigure(graph; columnSpacing = 0)
        @test_throws MethodError graphFigure(graph; nodeSpacing = 80)
        @test_throws MethodError graphFigure(graph; unaryFactorGap = 90)
        @test_throws MethodError graphFigure(graph; priorColumnGap = 210)
        @test_throws MethodError graphFigure(graph; margin = 48)
        @test_throws MethodError graphFigure(graph; nodeGap = 10)
        @test_throws MethodError graphFigure(tree; treeLayout = false)
    end
end
