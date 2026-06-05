include("setup.jl")

@testset verbose = true "Discrete factor graph" begin
    @testset "Constructs graph from variables and factors" begin
        variables = [
            DiscreteVariable(
                :x1,
                2;
                label = "x1",
                states = [:off, :on],
                probability = [0.7, 0.3]
            ),
            DiscreteVariable(:x2, 3; label = "x2", states = ["low", "mid", "high"])
        ]
        factors = [
            DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1", initialize = true),
            DiscreteFactor(
                :x1,
                :x2,
                [
                    1.0 0.2 0.4
                    0.1 1.1 0.3
                ];
                label = "x1_x2"
            )
        ]

        graph = factorGraph(variables, factors)

        @test graph isa DiscreteFactorGraph
        @test variableIndex(graph, :x1) == 1
        @test variableIndex(graph, "x2") == 2
        @test variableDimension(graph, :x2) == 3
        @test graph.variables[1].states == [:off, :on]
        @test stateIndex(graph, :x1, :on) == 2
        @test stateIndex(graph, "x2", "mid") == 2
        @test stateValue(graph, :x2, 3) == "high"
        @test factorIndex(graph, "x1_x2") == 2
        @test length(graph.edges) == 3
        @test graph.variableEdges[1] == [1, 2]
        @test graph.variableEdges[2] == [3]
        @test graph.factorEdges[1] == [1]
        @test graph.factorEdges[2] == [2, 3]
        @test edgeIndex(graph; variable = :x1, factor = "x1_x2") == 2
        @test edgeIndices(graph; variable = :x1) == [1, 2]
        @test edgeIndices(graph; factor = "x1_x2") == [2, 3]
    end

    @testset "Default states use integer references" begin
        graph = DiscreteFactorGraph(
            [DiscreteVariable(:x, 3; label = "x")],
            [DiscreteFactor(:x, [0.3, 0.4, 0.3])]
        )
        variable = graph.variables[1]

        @test variable.states == [1, 2, 3]
        @test stateIndex(graph, :x, 2) == 2
        @test stateValue(graph, :x, 3) == 3
    end

    @testset "Default discrete variable labels preserve IDs" begin
        graph = DiscreteFactorGraph(
            [DiscreteVariable(:x_1, 2)],
            [DiscreteFactor(:x_1, [0.3, 0.7])]
        )

        @test graph.variables[1].label == "x_1"
        @test graph.factors[1].label == "f1"

        explicit = DiscreteVariable(:x_2, 2; label = "x_2")
        @test explicit.label == "x_2"
    end

    @testset "Builds graph incrementally" begin
        graph = DiscreteFactorGraph()

        x1 = addVariable!(graph, :x1, 2; label = "x1")
        x2 = addVariable!(graph, DiscreteVariable(:x2, 2; label = "x2"))
        prior = addFactor!(graph, :x1, [0.6, 0.4]; label = "prior_x1")
        link = addFactor!(
            graph,
            DiscreteFactor(:x1, :x2, [1.0 0.2; 0.2 1.0]; label = "x1_x2")
        )

        @test x1 === graph.variables[1]
        @test x2 === graph.variables[2]
        @test prior === graph.factors[1]
        @test link === graph.factors[2]
        @test length(graph.edges) == 3
        @test graph.topologyVersion == 4
    end

    @testset "Discrete variable nodes can be used to construct factors" begin
        x1 = DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on])
        x2 = DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])

        prior = DiscreteFactor(x1, [0.8, 0.2]; label = "prior_x1", initialize = true)
        link = DiscreteFactor(x1, x2, [1.0 0.2; 0.1 0.9]; label = "link")
        graph = DiscreteFactorGraph([x1, x2], [prior, link])

        @test prior.variables == [:x1]
        @test link.variables == [:x1, :x2]
        @test variableIndex(graph, :x1) == 1
        @test variableDimension(graph, :x2) == 2
        @test stateIndex(graph, :x1, :on) == 2
        @test stateValue(graph, :x2, 2) == :high
        @test factorAxis(graph; factor = "link", variable = :x2) == 2
        @test factorAxes(graph, "link") == [1, 2]
        @test edgeIndex(graph; variable = :x2, factor = "link") == 3
        @test edgeIndices(graph; variable = :x1) == [1, 2]

        sumProductInference = sumproduct(graph)
        minSumInference = minsum(graph)

        @test length(marginal(graph, sumProductInference, :x1)) == 2
        @test marginalProbability(graph, sumProductInference, :x2, :high) isa Float64
        @test estimate(graph, minSumInference, :x1) in x1.states
        @test_throws ErrorException DiscreteFactor(GaussianVariable(:g1, 1), [0.8, 0.2])
    end

    @testset "Updates table and initialize flag" begin
        graph = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1"),
                DiscreteVariable(:x2, 2; label = "x2")
            ],
            [
                DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1"),
                DiscreteFactor(:x1, :x2, [1.0 0.2; 0.2 1.0]; label = "x1_x2")
            ]
        )

        updated = updateFactor!(
            graph;
            factor = "prior_x1",
            table = [0.2, 0.8],
            initialize = true
        )

        @test updated === graph.factors[1]
        @test updated.table == [0.2, 0.8]
        @test updated.initialize

        link = updateFactor!(graph, "x1_x2"; table = [0.9 0.1; 0.3 0.7])
        @test link.table == [0.9 0.1; 0.3 0.7]
        @test !link.initialize
    end

    @testset "Rejects invalid models atomically" begin
        @test_throws ErrorException factorGraph(
            [DiscreteVariable(:x1, 2; label = "x1")],
            [DiscreteFactor(:x1, [0.4, 0.3, 0.3]; label = "bad_prior")]
        )

        @test_throws ErrorException factorGraph(
            [DiscreteVariable(:x1, 2; label = "x1")],
            [DiscreteFactor(:x1, [-0.1, 1.1]; label = "bad_prior")]
        )

        @test_throws ErrorException DiscreteVariable(:x1, 2; states = [:off, :off])
        @test_throws ErrorException DiscreteVariable(:x1, 1; states = [""])
        @test_throws ErrorException DiscreteVariable(:x1, 2; states = [:off])
        @test_throws ErrorException DiscreteVariable(:x1, 0)
        @test_throws ErrorException DiscreteVariable(:x1, 2; label = "")
        @test_throws ErrorException DiscreteVariable(:x1, 2; states = :off)
        @test_throws ErrorException DiscreteVariable(:x1, 2; states = Symbol[])
        @test_throws ErrorException DiscreteVariable(:x1, 2; states = [1.0, 2.0])
        @test_throws ErrorException DiscreteVariable(:x1, 2; probability = [0.5])
        @test_throws ErrorException DiscreteVariable(:x1, 2; probability = [-0.1, 1.1])
        @test_throws ErrorException DiscreteVariable(:x1, 2; probability = [0.0, 0.0])
        @test_throws ErrorException DiscreteFactor()
        @test_throws ErrorException DiscreteFactor(:x1, [])
        @test_throws ErrorException DiscreteFactor(:x1, [-0.1, 1.1])

        @test_throws ErrorException factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1"),
                DiscreteVariable(:x2, 2; label = "x2")
            ],
            [
                DiscreteFactor(:x1, [0.6, 0.4]; label = "init_1", initialize = true),
                DiscreteFactor(:x1, [0.5, 0.5]; label = "init_2", initialize = true)
            ]
        )

        @test_throws ErrorException factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x"),
                DiscreteVariable(:x2, 2; label = "x")
            ],
            [
                DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1"),
                DiscreteFactor(:x2, [0.5, 0.5]; label = "prior_x2")
            ]
        )
        @test_throws ErrorException factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1"),
                DiscreteVariable(:x2, 2; label = "x2")
            ],
            [
                DiscreteFactor(:x1, [0.6, 0.4]; label = "prior"),
                DiscreteFactor(:x2, [0.5, 0.5]; label = "prior")
            ]
        )
        @test_throws ErrorException factorGraph(
            [DiscreteVariable(:x1, 2; label = "x1")],
            [DiscreteFactor(:x1, :x1, [0.8 0.2; 0.2 0.8]; label = "duplicate")]
        )
        @test_throws ErrorException factorGraph(
            [DiscreteVariable(:x1, 2; label = "x1")],
            [DiscreteFactor(:x1, :missing, [0.8 0.2; 0.2 0.8]; label = "missing")]
        )

        graph = factorGraph(
            [DiscreteVariable(:x1, 2; label = "x1")],
            [DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1")]
        )
        original = graph.factors[1]

        @test_throws ErrorException updateFactor!(graph; factor = "prior_x1", table = [1.0])
        @test graph.factors[1] === original
        @test_throws ErrorException factorIndex(graph, "missing")
        @test_throws ErrorException edgeIndex(graph; variable = :x1, factor = "missing")
        @test_throws ErrorException edgeIndex(graph; variable = :missing, factor = "prior_x1")
        @test_throws ErrorException edgeIndices(graph)
        @test_throws ErrorException stateIndex(graph, :x1, :missing)
        @test_throws ErrorException stateValue(graph, :x1, 0)
    end

    @testset "Constructs tree view" begin
        variables = [
            DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
            DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high]),
            DiscreteVariable(:x3, 2; label = "x3")
        ]
        factors = [
            DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1"),
            DiscreteFactor(:x1, :x2, [1.0 0.2; 0.2 1.0]; label = "x1_x2"),
            DiscreteFactor(:x2, :x3, [0.8 0.1; 0.1 0.8]; label = "x2_x3")
        ]

        tree = factorGraph(variables, factors; root = :x1)

        @test tree isa TreeFactorGraph{DiscreteFactorGraph}
        @test tree.rootVariableIndex == variableIndex(tree, :x1)
        @test length(tree.forwardOrder) == length(tree.graph.edges)
        @test length(tree.backwardOrder) == length(tree.graph.edges)
        @test sort(tree.forwardOrder) == collect(eachindex(tree.graph.edges))
        @test reverse(tree.forwardOrder) == tree.backwardOrder
        @test variableDimension(tree, :x2) == 2
        @test stateIndex(tree, :x1, :on) == 2
        @test stateValue(tree, :x2, 1) == :low
        @test factorIndex(tree, "x2_x3") == 3
        @test edgeIndex(tree; variable = :x2, factor = "x1_x2") == 3

        sameTree = treeFactorGraph(tree.graph; root = "x3")
        @test sameTree.rootVariableIndex == variableIndex(tree.graph, "x3")
    end

    @testset "Tree view refreshes after graph-only leaf addition" begin
        tree = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1"),
                DiscreteVariable(:x2, 2; label = "x2")
            ],
            [
                DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1"),
                DiscreteFactor(:x1, :x2, [1.0 0.2; 0.2 1.0]; label = "x1_x2")
            ];
            root = :x1
        )
        oldEdgeCount = length(tree.graph.edges)

        addVariable!(tree, :x3, 2; label = "x3", states = [:a, :b])
        @test length(tree.forwardOrder) == oldEdgeCount

        added = addFactor!(tree, :x2, :x3, [0.9 0.1; 0.1 0.9]; label = "x2_x3")

        @test added === tree.graph.factors[end]
        @test length(tree.graph.edges) == oldEdgeCount + 2
        @test length(tree.forwardOrder) == length(tree.graph.edges)
        @test length(tree.backwardOrder) == length(tree.graph.edges)
        @test sort(tree.forwardOrder) == collect(eachindex(tree.graph.edges))
    end

    @testset "Rejects invalid discrete tree views" begin
        cyclic = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1"),
                DiscreteVariable(:x2, 2; label = "x2")
            ],
            [
                DiscreteFactor(:x1, :x2, [1.0 0.2; 0.2 1.0]; label = "x1_x2_a"),
                DiscreteFactor(:x1, :x2, [0.8 0.1; 0.1 0.8]; label = "x1_x2_b")
            ]
        )

        disconnected = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1"),
                DiscreteVariable(:x2, 2; label = "x2")
            ],
            [
                DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1"),
                DiscreteFactor(:x2, [0.5, 0.5]; label = "prior_x2")
            ]
        )

        @test_throws ErrorException treeFactorGraph(cyclic; root = :x1)
        @test_throws ErrorException treeFactorGraph(disconnected; root = :x1)
    end
end
