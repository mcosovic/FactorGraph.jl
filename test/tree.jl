include("setup.jl")

function assertTreeMarginalsMatchWLS(
    graph::GaussianFactorGraph,
    inference;
    meanAtol::Float64 = 1e-7,
    covarianceAtol::Float64 = 1e-6
)
    expected = solveWLS(graph)

    @test length(marginals(inference)) == length(graph.variables)

    for variableIndex in eachindex(graph.variables)
        variable = graph.variables[variableIndex]
        marginalData = marginal(graph, inference, variable.id)

        @test isapprox(
            marginalData.mean,
            expected.variableMean[variableIndex];
            atol = meanAtol
        )
        @test isapprox(
            marginalData.covariance,
            expected.variableCovariance[variableIndex];
            atol = covarianceAtol
        )
    end

    return nothing
end

function robotLocalizationTree()
    variables = [
        GaussianVariable(:x1, 1),
        GaussianVariable(:x2, 1),
        GaussianVariable(:x3, 1),
        GaussianVariable(:x4, 1)
    ]

    factors = [
        GaussianFactor(:x1, 0.0, 1.0, 0.01; label = "prior"),
        GaussianFactor(:x1, :x2, 1.05, [-1.0 1.0], 0.04; label = "odom12"),
        GaussianFactor(:x2, :x3, 0.95, [-1.0 1.0], 0.04; label = "odom23"),
        GaussianFactor(:x3, :x4, 1.10, [-1.0 1.0], 0.04; label = "odom34"),
        GaussianFactor(:x4, 3.05, 1.0, 0.02; label = "landmark")
    ]

    return factorGraph(variables, factors; root = :x1)
end

@testset verbose = true "Tree GaussianFactor graph" begin
    @testset "Construction" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)

        @test tree.graph === graph
        @test tree.rootVariableIndex == variableIndex(graph, :x1)
        @test tree.variableParentEdge[tree.rootVariableIndex] === nothing
        @test length(tree.forwardOrder) == length(graph.edges)
        @test length(tree.backwardOrder) == length(graph.edges)
        @test tree.forwardOrder == reverse(tree.backwardOrder)
        @test sort(tree.forwardOrder) == collect(eachindex(graph.edges))
        @test tree.forwardIndex == 1
        @test tree.backwardIndex == 1

        schedule = forwardBackwardSchedule(tree)
        @test schedule.forwardIndex == 1
        @test schedule.backwardIndex == 1
    end

    @testset "Root keyword constructs tree view" begin
        graph = gaussianTreeTestGraph()
        tree = factorGraph(graph.variables, graph.factors; root = :x1)

        @test tree isa TreeFactorGraph
        @test tree.rootVariableIndex == variableIndex(tree.graph, :x1)
        @test tree.forwardOrder == treeFactorGraph(tree.graph; root = :x1).forwardOrder

        inference = canonical(tree; mean = 0.0, covariance = 1e6)
        forwardBackward!(tree, inference)
        assertTreeMarginalsMatchWLS(tree.graph, inference)

        reference = solveWLS(tree)
        @test reference.variableMean == solveWLS(tree.graph).variableMean
        @test marginalMean(tree, inference, :x1) == marginalMean(tree.graph, inference, :x1)
        @test marginalCovariance(tree, inference, :x1) ==
              marginalCovariance(tree.graph, inference, :x1)
        @test maxMeanError(tree, inference, reference) ==
              maxMeanError(tree.graph, inference, reference)
        @test maxCovarianceError(tree, inference, reference) ==
              maxCovarianceError(tree.graph, inference, reference)
        @test variableIndex(tree, :x1) == variableIndex(tree.graph, :x1)
        @test variableDimension(tree, :x2) == variableDimension(tree.graph, :x2)
        @test factorIndex(tree, "prior_x1") == factorIndex(tree.graph, "prior_x1")
        @test edgeIndices(tree; variable = :x1) == edgeIndices(tree.graph; variable = :x1)
        @test coefficientBlock(tree; factor = "factor_x1_x2", variable = :x2) ==
              coefficientBlock(tree.graph; factor = "factor_x1_x2", variable = :x2)
    end

    @testset "Root selection" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = "x3")

        @test tree.rootVariableIndex == variableIndex(graph, "x3")
        @test tree.variableParentEdge[tree.rootVariableIndex] === nothing
    end

    @testset "Rejects non-tree graph" begin
        @test_throws ErrorException treeFactorGraph(gaussianTestGraph(); root = :x1)
    end

    @testset "Controlled forward/backward steps" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = canonical(
            graph;
            mean = 0.0,
            covariance = 1e6
        )
        schedule = forwardBackwardSchedule(tree)

        for edgeId in tree.forwardOrder
            @test forwardStep!(tree, inference, schedule) == edgeId
        end

        @test forwardStep!(tree, inference, schedule) === nothing

        for edgeId in tree.backwardOrder
            @test backwardStep!(tree, inference, schedule) == edgeId
        end

        @test backwardStep!(tree, inference, schedule) === nothing

        reset!(schedule)
        @test schedule.forwardIndex == 1
        @test schedule.backwardIndex == 1
    end

    @testset "Controlled moment forward/backward steps" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = moment(
            graph;
            mean = 0.0,
            covariance = 1e6
        )
        schedule = forwardBackwardSchedule(tree)

        for edgeId in tree.forwardOrder
            @test forwardStep!(tree, inference, schedule) == edgeId
        end

        @test forwardStep!(tree, inference, schedule) === nothing

        for edgeId in tree.backwardOrder
            @test backwardStep!(tree, inference, schedule) == edgeId
        end

        @test backwardStep!(tree, inference, schedule) === nothing
        marginals!(graph, inference)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Tree cursor forward/backward steps" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = canonical(graph; mean = 0.0, covariance = 1e6)

        for edgeId in tree.forwardOrder
            @test forwardStep!(tree, inference) == edgeId
        end

        @test forwardStep!(tree, inference) === nothing

        for edgeId in tree.backwardOrder
            @test backwardStep!(tree, inference) == edgeId
        end

        @test backwardStep!(tree, inference) === nothing

        reset!(tree)
        @test tree.forwardIndex == 1
        @test tree.backwardIndex == 1
    end

    @testset "Selected forward/backward steps" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        stepped = canonical(graph; mean = 0.0, covariance = 1e6)
        selected = canonical(graph; mean = 0.0, covariance = 1e6)

        firstForwardEdge = tree.forwardOrder[1]
        @test forwardStep!(tree, stepped) == firstForwardEdge
        @test forwardStep!(tree, selected, firstForwardEdge) == firstForwardEdge
        @test selected.variableToFactor[firstForwardEdge].information ==
              stepped.variableToFactor[firstForwardEdge].information
        @test selected.variableToFactor[firstForwardEdge].precision ==
              stepped.variableToFactor[firstForwardEdge].precision
        @test selected.factorToVariable[firstForwardEdge].information ==
              stepped.factorToVariable[firstForwardEdge].information
        @test selected.factorToVariable[firstForwardEdge].precision ==
              stepped.factorToVariable[firstForwardEdge].precision

        reset!(tree)

        firstBackwardEdge = tree.backwardOrder[1]
        edge = graph.edges[firstBackwardEdge]
        variable = graph.variables[edge.variableIndex].id
        factor = graph.factors[edge.factorIndex].label

        @test backwardStep!(tree, stepped) == firstBackwardEdge
        @test backwardStep!(tree, selected; variable = variable, factor = factor) ==
              firstBackwardEdge
        @test selected.variableToFactor[firstBackwardEdge].information ==
              stepped.variableToFactor[firstBackwardEdge].information
        @test selected.variableToFactor[firstBackwardEdge].precision ==
              stepped.variableToFactor[firstBackwardEdge].precision
        @test selected.factorToVariable[firstBackwardEdge].information ==
              stepped.factorToVariable[firstBackwardEdge].information
        @test selected.factorToVariable[firstBackwardEdge].precision ==
              stepped.factorToVariable[firstBackwardEdge].precision
    end

    @testset "Tree freeze helpers delegate to Gaussian inference state" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = canonical(tree; mean = 0.0, covariance = 1e6)

        edgeId = tree.forwardOrder[1]
        edge = graph.edges[edgeId]
        variable = graph.variables[edge.variableIndex].id
        factor = graph.factors[edge.factorIndex].label
        variableMessage = copy(inference.variableToFactor[edgeId].information)
        factorMessage = copy(inference.factorToVariable[edgeId].information)

        freezeEdge!(tree, inference; variable = variable, factor = factor)
        @test isFrozenEdge(tree, inference; variable = variable, factor = factor)
        @test forwardStep!(tree, inference, edgeId) == edgeId
        @test inference.variableToFactor[edgeId].information == variableMessage
        @test inference.factorToVariable[edgeId].information == factorMessage

        unfreezeEdge!(tree, inference; variable = variable, factor = factor)
        @test !isFrozenEdge(tree, inference; variable = variable, factor = factor)
    end

    @testset "Tree damping helpers delegate to Gaussian inference state" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = canonical(tree; mean = 0.0, covariance = 1e6)

        edgeId = tree.forwardOrder[1]
        edge = graph.edges[edgeId]
        variable = graph.variables[edge.variableIndex].id
        factor = graph.factors[edge.factorIndex].label

        @test !areDampedEdges(tree, inference; variable = variable, factor = factor)

        dampEdges!(
            tree,
            inference;
            variable = variable,
            factor = factor,
            prob = 1.0,
            alpha = 0.25
        )
        @test areDampedEdges(tree, inference; variable = variable, factor = factor)
        @test inference.edgeDampingProb[edgeId] == 1.0
        @test inference.edgeDampingAlpha[edgeId] == 0.25

        undampEdges!(tree, inference; variable = variable, factor = factor)
        @test !areDampedEdges(tree, inference; variable = variable, factor = factor)

        dampEdges!(tree, inference; prob = 1.0, alpha = 0.3)
        @test all(inference.dampedEdges)
        @test all(==(0.3), inference.edgeDampingAlpha)

        undampEdges!(tree, inference)
        @test !any(inference.dampedEdges)
    end

    @testset "Tree warm-start addition refreshes view" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = moment(tree; mean = 0.0, covariance = 1e6)

        forwardBackward!(tree, inference)
        oldEdgeCount = length(tree.graph.edges)

        addVariable!(tree, inference, :x4, 1; label = "x4")
        @test length(tree.graph.variables) == 4
        @test length(tree.forwardOrder) == oldEdgeCount

        addFactor!(tree, inference, :x3, :x4, 0.8, [-1.0 1.0], 0.3; label = "x3_x4")

        @test length(tree.graph.edges) == oldEdgeCount + 2
        @test length(tree.forwardOrder) == length(tree.graph.edges)
        @test length(tree.backwardOrder) == length(tree.graph.edges)
        @test sort(tree.forwardOrder) == collect(eachindex(tree.graph.edges))

        forwardBackward!(tree, inference)
        assertTreeMarginalsMatchWLS(
            tree.graph,
            inference;
            meanAtol = 1e-6,
            covarianceAtol = 2e-6
        )
    end

    @testset "Tree warm-start GaussianVariable uses stored inference defaults" begin
        tree = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2")
            ],
            [
                GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "prior_x1"),
                GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.2; label = "x1_x2")
            ];
            root = :x1
        )
        fullTree = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2"),
                GaussianVariable(:x3, 1; label = "x3")
            ],
            [
                GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "prior_x1"),
                GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.2; label = "x1_x2"),
                GaussianFactor(:x2, :x3, 0.8, [-1.0 1.0], 0.3; label = "x2_x3")
            ];
            root = :x1
        )
        inference = canonical(tree; mean = 4.0, covariance = 9.0)
        full = canonical(fullTree; mean = 4.0, covariance = 9.0)

        addVariable!(tree, inference, :x3, 1; label = "x3")
        addFactor!(tree, inference, :x2, :x3, 0.8, [-1.0 1.0], 0.3; label = "x2_x3")

        x3Edge = tree.graph.variableEdges[variableIndex(tree, :x3)][1]

        @test isapprox(inference.initial[end].information, [4.0 / 9.0])
        @test isapprox(inference.initial[end].precision, [1.0 / 9.0;;])
        @test isapprox(inference.variableToFactor[x3Edge].information, [4.0 / 9.0])
        @test isapprox(inference.variableToFactor[x3Edge].precision, [1.0 / 9.0;;])
        @test isempty(inference.nextVariableToFactor)
        @test length(tree.forwardOrder) == length(tree.graph.edges)

        forwardBackward!(tree, inference)
        forwardBackward!(fullTree, full)

        for GaussianVariable in (:x1, :x2, :x3)
            @test isapprox(
                marginalMean(tree, inference, GaussianVariable),
                marginalMean(fullTree, full, GaussianVariable);
                atol = 1e-7
            )
            @test isapprox(
                marginalCovariance(tree, inference, GaussianVariable),
                marginalCovariance(fullTree, full, GaussianVariable);
                atol = 1e-7
            )
        end
    end

    @testset "Dynamic tree full sweep matches WLS after adding leaf" begin
        tree = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2"),
                GaussianVariable(:x3, 1; label = "x3")
            ],
            [
                GaussianFactor(:x1, 0.0, 1.0, 0.1; label = "prior_x1"),
                GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.2; label = "x1_x2"),
                GaussianFactor(:x2, :x3, 0.8, [-1.0 1.0], 0.2; label = "x2_x3")
            ];
            root = :x1
        )
        inference = canonical(tree; mean = 0.0, covariance = 1e6)

        forwardBackward!(tree, inference)
        addVariable!(tree, inference, :x4, 1; label = "x4")
        addFactor!(tree, inference, :x3, :x4, 1.1, [-1.0 1.0], 0.3; label = "x3_x4")

        forwardBackward!(tree, inference)
        assertTreeMarginalsMatchWLS(
            tree.graph,
            inference;
            meanAtol = 2e-6,
            covarianceAtol = 2e-6
        )
    end

    @testset "Robot-style incremental terminal update matches full sweep at new state" begin
        incrementalTree = robotLocalizationTree()
        fullTree = robotLocalizationTree()
        incremental = canonical(incrementalTree; mean = 0.0, covariance = 1e6)
        full = canonical(fullTree; mean = 0.0, covariance = 1e6)

        forwardBackward!(incrementalTree, incremental)
        forwardBackward!(fullTree, full)

        addVariable!(incrementalTree, incremental, :x5, 1)
        addFactor!(
            incrementalTree,
            incremental,
            :x4,
            :x5,
            0.90,
            [-1.0 1.0],
            0.04;
            label = "odom45"
        )
        addVariable!(fullTree, full, :x5, 1)
        addFactor!(
            fullTree,
            full,
            :x4,
            :x5,
            0.90,
            [-1.0 1.0],
            0.04;
            label = "odom45"
        )

        forwardStep!(incrementalTree, incremental; variable = :x5, factor = "odom45")
        backward!(incrementalTree, incremental)
        marginals!(incrementalTree, incremental)

        forwardBackward!(fullTree, full)

        incrementalX5 = marginal(incrementalTree.graph, incremental, :x5)
        fullX5 = marginal(fullTree.graph, full, :x5)

        @test isapprox(incrementalX5.mean, fullX5.mean; atol = 1e-7)
        @test isapprox(incrementalX5.covariance, fullX5.covariance; atol = 1e-7)
    end

    @testset "Forward and backward helpers" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = canonical(
            graph;
            mean = 0.0,
            covariance = 1e6
        )
        schedule = forwardBackwardSchedule(tree)

        @test forward!(tree, inference, schedule) === schedule
        @test schedule.forwardIndex == length(tree.forwardOrder) + 1
        @test schedule.backwardIndex == 1

        @test backward!(tree, inference, schedule) === schedule
        @test schedule.forwardIndex == length(tree.forwardOrder) + 1
        @test schedule.backwardIndex == length(tree.backwardOrder) + 1
        marginals!(graph, inference)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Moment marginals match WLS" begin
        graph = gaussianTreeTestGraph()
        inference = moment(
            graph;
            mean = 0.0,
            covariance = 1e6
        )

        gbp!(graph, inference; iterations = 6, schedule = :flooding)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Moment forward-backward matches WLS" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = moment(
            graph;
            mean = 0.0,
            covariance = 1e6
        )

        forwardBackward!(tree, inference)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Moment forward-backward with nondefault root" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = "x3")
        inference = moment(
            graph;
            mean = 0.0,
            covariance = 1e6
        )

        forwardBackward!(tree, inference)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Canonical marginals match WLS" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(
            graph;
            mean = 0.0,
            covariance = 1e6
        )

        gbp!(graph, inference; iterations = 6, schedule = :flooding)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Canonical forward-backward matches WLS" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = canonical(
            graph;
            mean = 0.0,
            covariance = 1e6
        )

        forwardBackward!(tree, inference)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Canonical forward-backward with nondefault root" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = "x3")
        inference = canonical(
            graph;
            mean = 0.0,
            covariance = 1e6
        )

        forwardBackward!(tree, inference)
        assertTreeMarginalsMatchWLS(graph, inference)
    end

    @testset "Gaussian tree updateFactor! delegates without changing topology" begin
        graph = gaussianTreeTestGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = canonical(tree; mean = 0.0, covariance = 1e6)
        oldForwardOrder = copy(tree.forwardOrder)
        oldBackwardOrder = copy(tree.backwardOrder)

        updateFactor!(tree, inference; factor = "prior_x1", mean = 0.75, covariance = 0.05)
        forwardBackward!(tree, inference)

        @test tree.forwardOrder == oldForwardOrder
        @test tree.backwardOrder == oldBackwardOrder
        assertTreeMarginalsMatchWLS(tree.graph, inference; meanAtol = 1e-6)
    end
end
