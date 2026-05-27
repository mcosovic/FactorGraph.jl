include("setup.jl")

if !isdefined(@__MODULE__, :messageSnapshots)
    function messageSnapshots(messages, edgeIds)
        return [
            (
                mean = copy(messages[edgeId].mean),
                covariance = copy(messages[edgeId].covariance)
            )
            for edgeId in edgeIds
        ]
    end
end

if !isdefined(@__MODULE__, :assertMessagesUnchanged)
    function assertMessagesUnchanged(messages, edgeIds, snapshots)
        for (snapshotIndex, edgeId) in pairs(collect(edgeIds))
            @test messages[edgeId].mean == snapshots[snapshotIndex].mean
            @test messages[edgeId].covariance == snapshots[snapshotIndex].covariance
        end

        return nothing
    end
end

@testset verbose = true "GaussianFactor additions" begin
    @testset "Gaussian inference common supertype" begin
        @test GaussianMomentInference <: GaussianInference
        @test GaussianCanonicalInference <: GaussianInference
        @test GaussianMinSumInference <: GaussianInference
        @test GaussianMomentInference <: GaussianSumProductInference
        @test GaussianCanonicalInference <: GaussianSumProductInference
        @test GaussianMomentInference <: AbstractSumProductInference
        @test GaussianCanonicalInference <: AbstractSumProductInference
        @test GaussianMinSumInference <: AbstractMinSumInference
    end

    @testset "Builds graph incrementally" begin
        graph = GaussianFactorGraph()

        x1 = addVariable!(graph, :x1, 1; label = "x1")
        x2 = addVariable!(graph, GaussianVariable(:x2, 1; label = "x2"))
        prior = addFactor!(graph, :x1, 1.0, 1.0, 0.1; label = "prior_x1")
        link = addFactor!(graph, :x1, :x2, -1.0, [1.0 -1.0], 0.2; label = "link")

        @test x1 === graph.variables[1]
        @test x2 === graph.variables[2]
        @test prior === graph.factors[1]
        @test link === graph.factors[2]
        @test length(graph.edges) == 3
        @test variableIndex(graph, :x1) == 1
        @test variableIndex(graph, "x2") == 2
        @test factorIndex(graph, "link") == 2

        inference = canonical(graph)
        gbp!(graph, inference; iterations = 20)
        @test isfinite(maxMeanError(graph, inference, solveWLS(graph)))
    end

    @testset "Gaussian addVariable convenience accepts components" begin
        graph = GaussianFactorGraph()
        variable = addVariable!(graph, :v, 2; label = "v", components = [:re, :im])

        @test variable.components == [:re, :im]
        @test componentIndex(graph, :v, :re) == 1
        @test componentIndex(graph, :v, :im) == 2
        @test componentValue(graph, :v, 2) == :im
    end

    @testset "Rejects duplicate incremental variables" begin
        graph = GaussianFactorGraph()
        addVariable!(graph, :x1, 1; label = "x1")

        @test_throws ErrorException addVariable!(graph, :x1, 1; label = "x1_duplicate")
        @test_throws ErrorException addVariable!(graph, :x2, 1; label = "x1")
        @test length(graph.variables) == 1
        @test isempty(graph.edges)
    end

    @testset "Adds variables after factors during graph-only construction" begin
        graph = GaussianFactorGraph()
        addVariable!(graph, :x1, 1; label = "x1")
        addFactor!(graph, :x1, 1.0, 1.0, 0.1; label = "prior_x1")

        x2 = addVariable!(graph, :x2, 1; label = "x2")
        link = addFactor!(graph, :x1, :x2, 0.5, [-1.0 1.0], 0.2; label = "link")

        @test x2 === graph.variables[2]
        @test link === graph.factors[2]
        @test length(graph.variables) == 2
        @test length(graph.factors) == 2
        @test length(graph.edges) == 3
        @test graph.variableEdges[variableIndex(graph, :x2)] == [3]

        inference = canonical(graph)
        gbp!(graph, inference; iterations = 20)
        @test isfinite(maxMeanError(graph, inference, solveWLS(graph)))
    end

    @testset "Adds labeled GaussianFactor and updates adjacency" begin
        graph = gaussianTreeTestGraph()

        oldFactorCount = length(graph.factors)
        oldEdgeCount = length(graph.edges)
        oldX3Edges = copy(graph.variableEdges[variableIndex(graph, :x3)])

        added = addFactor!(
            graph,
            GaussianFactor(:x3, 0.9, 1.0, 0.3; label = "prior_x3")
        )

        @test added === graph.factors[end]
        @test added.id == oldFactorCount + 1
        @test added.label == "prior_x3"
        @test added.mean == [0.9]
        @test length(graph.factors) == oldFactorCount + 1
        @test length(graph.edges) == oldEdgeCount + 1
        @test graph.factorEdges[end] == [oldEdgeCount + 1]
        @test graph.edges[end].factorIndex == oldFactorCount + 1
        @test graph.edges[end].variableIndex == variableIndex(graph, :x3)
        @test graph.variableEdges[variableIndex(graph, :x3)] ==
            [oldX3Edges; oldEdgeCount + 1]

        inference = canonical(graph; mean = 0.0, covariance = 1e6)
        gbp!(graph, inference; iterations = 20)
        @test length(marginals(inference)) == length(graph.variables)
    end

    @testset "Adds GaussianFactor from argument list with default label" begin
        graph = gaussianTreeTestGraph()
        nextFactorId = length(graph.factors) + 1

        added = addFactor!(graph, :x1, 0.1, 1.0, 0.4)

        @test added.id == nextFactorId
        @test added.label == "f_$nextFactorId"
        @test factorIndex(graph, added.label) == nextFactorId
    end

    @testset "Adds GaussianFactor and extends moment inference as warm start" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph; mean = 0.0, covariance = 1e6)
        gbp!(graph, inference; iterations = 5, schedule = :flooding)

        oldEdgeCount = length(graph.edges)
        oldFactorCount = length(graph.factors)
        oldVariableToFactor = messageSnapshots(inference.variableToFactor, 1:oldEdgeCount)
        oldFactorToVariable = messageSnapshots(inference.factorToVariable, 1:oldEdgeCount)

        added = addFactor!(
            graph,
            inference,
            :x3,
            0.9,
            1.0,
            0.3;
            label = "prior_x3"
        )

        @test added.id == oldFactorCount + 1
        @test length(graph.edges) == oldEdgeCount + 1
        @test length(inference.variableToFactor) == length(graph.edges)
        @test length(inference.factorToVariable) == length(graph.edges)
        @test length(inference.nextVariableToFactor) == oldEdgeCount
        @test length(inference.nextFactorToVariable) == oldEdgeCount
        @test length(inference.frozenEdges) == length(graph.edges)
        @test length(inference.dampedEdges) == length(graph.edges)
        @test length(inference.frozenFactors) == length(graph.factors)

        assertMessagesUnchanged(inference.variableToFactor, 1:oldEdgeCount, oldVariableToFactor)
        assertMessagesUnchanged(inference.factorToVariable, 1:oldEdgeCount, oldFactorToVariable)

        newEdge = graph.edges[end]
        initial = inference.initial[newEdge.variableIndex]
        @test inference.variableToFactor[end].mean == initial.mean
        @test inference.variableToFactor[end].covariance == initial.covariance
        @test inference.factorToVariable[end].mean == initial.mean
        @test inference.factorToVariable[end].covariance == initial.covariance

        gbp!(graph, inference; iterations = 20)
        @test length(marginals(inference)) == length(graph.variables)
        @test isfinite(maxMeanError(graph, inference, solveWLS(graph)))
    end

    @testset "Adds GaussianFactor and extends canonical inference as warm start" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(graph; mean = 0.0, covariance = 1e6)
        gbp!(graph, inference; iterations = 5)

        oldEdgeCount = length(graph.edges)
        oldFactorCount = length(graph.factors)
        oldInformation = [copy(message.information) for message in inference.variableToFactor]
        oldPrecision = [copy(message.precision) for message in inference.factorToVariable]

        added = addFactor!(
            graph,
            inference,
            GaussianFactor(:x3, 0.9, 1.0, 0.3; label = "prior_x3")
        )

        @test added.id == oldFactorCount + 1
        @test length(graph.edges) == oldEdgeCount + 1
        @test length(inference.variableToFactor) == length(graph.edges)
        @test length(inference.factorToVariable) == length(graph.edges)
        @test isempty(inference.nextVariableToFactor)
        @test isempty(inference.nextFactorToVariable)
        @test length(inference.frozenEdges) == length(graph.edges)
        @test length(inference.dampedEdges) == length(graph.edges)
        @test length(inference.frozenFactors) == length(graph.factors)

        for edgeId in 1:oldEdgeCount
            @test inference.variableToFactor[edgeId].information == oldInformation[edgeId]
            @test inference.factorToVariable[edgeId].precision == oldPrecision[edgeId]
        end

        newEdge = graph.edges[end]
        initial = inference.initial[newEdge.variableIndex]
        @test inference.variableToFactor[end].information == initial.information
        @test inference.variableToFactor[end].precision == initial.precision
        @test inference.factorToVariable[end].information == initial.information
        @test inference.factorToVariable[end].precision == initial.precision

        gbp!(graph, inference; iterations = 20)
        @test length(marginals(inference)) == length(graph.variables)
        @test isfinite(maxMeanError(graph, inference, solveWLS(graph)))
    end

    @testset "Adds GaussianVariable and GaussianFactor extending moment inference as warm start" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph; mean = 0.0, covariance = 1e6)
        gbp!(graph, inference; iterations = 5, schedule = :flooding)

        oldVariableCount = length(graph.variables)
        oldEdgeCount = length(graph.edges)
        oldVariableToFactor = messageSnapshots(inference.variableToFactor, 1:oldEdgeCount)

        addedVariable = addVariable!(graph, inference, :x4, 1; label = "x4")

        @test addedVariable === graph.variables[end]
        @test length(graph.variables) == oldVariableCount + 1
        @test variableIndex(graph, :x4) == oldVariableCount + 1
        @test graph.variableEdges[end] == Int[]
        @test length(inference.initial) == length(graph.variables)
        @test length(inference.marginal) == length(graph.variables)
        @test length(inference.frozenVariables) == length(graph.variables)
        assertMessagesUnchanged(inference.variableToFactor, 1:oldEdgeCount, oldVariableToFactor)

        addFactor!(graph, inference, :x3, :x4, 0.8, [-1.0 1.0], 0.3; label = "x3_x4")

        gbp!(graph, inference; iterations = 20)
        @test length(marginals(inference)) == length(graph.variables)
        @test isfinite(maxMeanError(graph, inference, solveWLS(graph)))
    end

    @testset "Added moment GaussianVariable uses stored inference defaults" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph; mean = 2.5, covariance = 7.0)

        addVariable!(graph, inference, :x4, 1; label = "x4", components = [:angle])

        @test graph.variables[end].components == [:angle]
        @test inference.initial[end].mean == [2.5]
        @test inference.initial[end].covariance == [7.0;;]
        @test inference.marginal[end].mean == [2.5]
        @test inference.marginal[end].covariance == [7.0;;]
    end

    @testset "Added moment GaussianVariable and edges use vector inference defaults" begin
        defaultMean = [1.0, -2.0]
        defaultCovariance = [3.0 0.2; 0.2 4.0]
        graph = factorGraph(
            [GaussianVariable(:x1, 2; label = "x1")],
            [GaussianFactor(:x1, [0.0, 0.0], Matrix{Float64}(I, 2, 2), [0.1, 0.1]; label = "prior_x1")]
        )
        inference = moment(
            graph;
            mean = defaultMean,
            covariance = defaultCovariance
        )

        addVariable!(graph, inference, :x2, 2; label = "x2")

        @test inference.initial[end].mean == defaultMean
        @test inference.initial[end].covariance == defaultCovariance
        @test inference.marginal[end].mean == defaultMean
        @test inference.marginal[end].covariance == defaultCovariance

        addFactor!(
            graph,
            inference,
            :x1,
            :x2,
            [0.3, -0.1],
            [
                1.0  0.0  -1.0   0.0
                0.0  1.0   0.0  -1.0
            ],
            [0.2, 0.3];
            label = "x1_x2"
        )

        newVariableEdge = graph.variableEdges[variableIndex(graph, :x2)][1]

        @test inference.variableToFactor[newVariableEdge].mean == defaultMean
        @test inference.variableToFactor[newVariableEdge].covariance == defaultCovariance
        @test inference.factorToVariable[newVariableEdge].mean == defaultMean
        @test inference.factorToVariable[newVariableEdge].covariance == defaultCovariance
        @test isempty(inference.nextVariableToFactor)
        @test isempty(inference.nextFactorToVariable)
    end

    @testset "Added moment GaussianVariable explicit prior overrides inference defaults" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph; mean = 2.5, covariance = 7.0)

        addVariable!(
            graph,
            inference,
            GaussianVariable(:x4, 1; label = "x4", mean = [9.0], covariance = [0.5;;])
        )

        @test inference.initial[end].mean == [9.0]
        @test inference.initial[end].covariance == [0.5;;]
        @test inference.marginal[end].mean == [9.0]
        @test inference.marginal[end].covariance == [0.5;;]
    end

    @testset "Moment inference uses initializing unary GaussianFactor as initial belief" begin
        graph = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1", mean = 9.0, covariance = 9.0),
                GaussianVariable(:x2, 1; label = "x2")
            ],
            [
                GaussianFactor(:x1, 1.5, 1.0, 0.25; label = "init_x1", initialize = true),
                GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.5; label = "link")
            ]
        )
        inference = moment(graph; mean = 0.0, covariance = 1e6)
        x1Index = variableIndex(graph, :x1)

        @test inference.initial[x1Index].mean == [1.5]
        @test inference.initial[x1Index].covariance == [0.25;;]

        for edgeId in graph.variableEdges[x1Index]
            @test inference.variableToFactor[edgeId].mean == [1.5]
            @test inference.variableToFactor[edgeId].covariance == [0.25;;]
            @test inference.factorToVariable[edgeId].mean == [1.5]
            @test inference.factorToVariable[edgeId].covariance == [0.25;;]
        end
    end

    @testset "Moment warm-start initializing unary GaussianFactor resets incident messages" begin
        graph = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2")
            ],
            [
                GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.5; label = "link")
            ]
        )
        inference = moment(graph; mean = 0.0, covariance = 1e6)

        addFactor!(
            graph,
            inference,
            :x1,
            2.0,
            1.0,
            0.4;
            label = "init_x1",
            initialize = true
        )

        x1Index = variableIndex(graph, :x1)

        @test isapprox(inference.initial[x1Index].mean, [2.0])
        @test isapprox(inference.initial[x1Index].covariance, [0.4;;])

        for edgeId in graph.variableEdges[x1Index]
            @test isapprox(inference.variableToFactor[edgeId].mean, [2.0])
            @test isapprox(inference.variableToFactor[edgeId].covariance, [0.4;;])
            @test isapprox(inference.factorToVariable[edgeId].mean, [2.0])
            @test isapprox(inference.factorToVariable[edgeId].covariance, [0.4;;])
        end
    end

    @testset "Adds GaussianVariable and GaussianFactor extending canonical inference as warm start" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(graph; mean = 0.0, covariance = 1e6)
        gbp!(graph, inference; iterations = 5)

        oldVariableCount = length(graph.variables)
        oldEdgeCount = length(graph.edges)
        oldInformation = [copy(message.information) for message in inference.variableToFactor]

        addedVariable = addVariable!(graph, inference, GaussianVariable(:x4, 1; label = "x4"))

        @test addedVariable === graph.variables[end]
        @test length(graph.variables) == oldVariableCount + 1
        @test variableIndex(graph, "x4") == oldVariableCount + 1
        @test graph.variableEdges[end] == Int[]
        @test length(inference.initial) == length(graph.variables)
        @test length(inference.marginal) == length(graph.variables)
        @test length(inference.frozenVariables) == length(graph.variables)

        for edgeId in 1:oldEdgeCount
            @test inference.variableToFactor[edgeId].information == oldInformation[edgeId]
        end

        addFactor!(graph, inference, :x3, :x4, 0.8, [-1.0 1.0], 0.3; label = "x3_x4")

        gbp!(graph, inference; iterations = 20)
        @test length(marginals(inference)) == length(graph.variables)
        @test isfinite(maxMeanError(graph, inference, solveWLS(graph)))
    end

    @testset "Added canonical GaussianVariable uses stored inference defaults" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(graph; mean = 2.5, covariance = 7.0)

        addVariable!(graph, inference, :x4, 1; label = "x4", components = [:angle])

        @test graph.variables[end].components == [:angle]
        @test isapprox(inference.initial[end].information, [2.5 / 7.0])
        @test isapprox(inference.initial[end].precision, [1.0 / 7.0;;])
        @test isapprox(inference.marginal[end].mean, [2.5])
        @test isapprox(inference.marginal[end].covariance, [7.0;;])
    end

    @testset "Added canonical GaussianVariable and edges use vector inference defaults" begin
        defaultMean = [1.0, -2.0]
        defaultCovariance = [3.0 0.2; 0.2 4.0]
        expectedPrecision = inv(Symmetric(defaultCovariance))
        expectedInformation = expectedPrecision * defaultMean
        graph = factorGraph(
            [GaussianVariable(:x1, 2; label = "x1")],
            [GaussianFactor(:x1, [0.0, 0.0], Matrix{Float64}(I, 2, 2), [0.1, 0.1]; label = "prior_x1")]
        )
        inference = canonical(
            graph;
            mean = defaultMean,
            covariance = defaultCovariance
        )

        addVariable!(graph, inference, :x2, 2; label = "x2")

        @test isapprox(inference.initial[end].information, expectedInformation)
        @test isapprox(inference.initial[end].precision, expectedPrecision)
        @test isapprox(inference.marginal[end].mean, defaultMean)
        @test isapprox(inference.marginal[end].covariance, defaultCovariance)

        addFactor!(
            graph,
            inference,
            :x1,
            :x2,
            [0.3, -0.1],
            [
                1.0  0.0  -1.0   0.0
                0.0  1.0   0.0  -1.0
            ],
            [0.2, 0.3];
            label = "x1_x2"
        )

        newVariableEdge = graph.variableEdges[variableIndex(graph, :x2)][1]

        @test isapprox(inference.variableToFactor[newVariableEdge].information, expectedInformation)
        @test isapprox(inference.variableToFactor[newVariableEdge].precision, expectedPrecision)
        @test isapprox(inference.factorToVariable[newVariableEdge].information, expectedInformation)
        @test isapprox(inference.factorToVariable[newVariableEdge].precision, expectedPrecision)
        @test isempty(inference.nextVariableToFactor)
        @test isempty(inference.nextFactorToVariable)
    end

    @testset "Added canonical GaussianVariable explicit prior overrides inference defaults" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(graph; mean = 2.5, covariance = 7.0)

        addVariable!(
            graph,
            inference,
            GaussianVariable(:x4, 1; label = "x4", mean = [9.0], covariance = [0.5;;])
        )

        @test isapprox(inference.initial[end].information, [18.0])
        @test isapprox(inference.initial[end].precision, [2.0;;])
        @test isapprox(inference.marginal[end].mean, [9.0])
        @test isapprox(inference.marginal[end].covariance, [0.5;;])
    end

    @testset "Canonical inference uses initializing unary GaussianFactor as initial belief" begin
        graph = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1", mean = 9.0, covariance = 9.0),
                GaussianVariable(:x2, 1; label = "x2")
            ],
            [
                GaussianFactor(:x1, 1.5, 1.0, 0.25; label = "init_x1", initialize = true),
                GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.5; label = "link")
            ]
        )
        inference = canonical(graph; mean = 0.0, covariance = 1e6)
        x1Index = variableIndex(graph, :x1)
        expectedInformation = [1.5 / 0.25]
        expectedPrecision = [1.0 / 0.25;;]

        @test isapprox(inference.initial[x1Index].information, expectedInformation)
        @test isapprox(inference.initial[x1Index].precision, expectedPrecision)
        @test isapprox(inference.marginal[x1Index].mean, [1.5])
        @test isapprox(inference.marginal[x1Index].covariance, [0.25;;])

        for edgeId in graph.variableEdges[x1Index]
            @test isapprox(inference.variableToFactor[edgeId].information, expectedInformation)
            @test isapprox(inference.variableToFactor[edgeId].precision, expectedPrecision)
            @test isapprox(inference.factorToVariable[edgeId].information, expectedInformation)
            @test isapprox(inference.factorToVariable[edgeId].precision, expectedPrecision)
        end
    end

    @testset "Canonical warm-start initializing unary GaussianFactor resets incident messages" begin
        graph = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2")
            ],
            [
                GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 0.5; label = "link")
            ]
        )
        inference = canonical(graph; mean = 0.0, covariance = 1e6)

        addFactor!(
            graph,
            inference,
            :x1,
            2.0,
            1.0,
            0.4;
            label = "init_x1",
            initialize = true
        )

        x1Index = variableIndex(graph, :x1)
        expectedInformation = [2.0 / 0.4]
        expectedPrecision = [1.0 / 0.4;;]

        @test isapprox(inference.initial[x1Index].information, expectedInformation)
        @test isapprox(inference.initial[x1Index].precision, expectedPrecision)
        @test isapprox(inference.marginal[x1Index].mean, [2.0])
        @test isapprox(inference.marginal[x1Index].covariance, [0.4;;])
        @test isempty(inference.nextVariableToFactor)
        @test isempty(inference.nextFactorToVariable)

        for edgeId in graph.variableEdges[x1Index]
            @test isapprox(inference.variableToFactor[edgeId].information, expectedInformation)
            @test isapprox(inference.variableToFactor[edgeId].precision, expectedPrecision)
            @test isapprox(inference.factorToVariable[edgeId].information, expectedInformation)
            @test isapprox(inference.factorToVariable[edgeId].precision, expectedPrecision)
        end
    end

    @testset "Rejects warm-start addition when inference no longer matches graph" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph)
        addFactor!(graph, :x3, 0.9, 1.0, 0.3; label = "prior_x3")

        oldFactorCount = length(graph.factors)
        oldEdgeCount = length(graph.edges)

        @test_throws ErrorException addFactor!(
            graph,
            inference,
            :x2,
            0.5,
            1.0,
            0.3;
            label = "prior_x2_extra"
        )
        @test length(graph.factors) == oldFactorCount
        @test length(graph.edges) == oldEdgeCount
    end

    @testset "Rejects stale inference after graph-only topology changes" begin
        graph = gaussianTreeTestGraph()
        momentInferenceState = moment(graph)
        canonicalInferenceState = canonical(graph)

        addVariable!(graph, :x4, 1; label = "x4")
        addFactor!(graph, :x3, :x4, 0.8, [-1.0 1.0], 0.3; label = "x3_x4")

        @test_throws ErrorException gbp!(graph, momentInferenceState; iterations = 1)
        @test_throws ErrorException marginals!(graph, canonicalInferenceState)
        @test_throws ErrorException marginalMean(graph, canonicalInferenceState, :x1)
        @test_throws ErrorException freezeFactor!(graph, momentInferenceState, "prior_x1")
        @test_throws ErrorException dampEdges!(graph, canonicalInferenceState; variable = :x1)
    end

    @testset "Rejects invalid additions atomically" begin
        graph = gaussianTreeTestGraph()
        oldFactorCount = length(graph.factors)
        oldEdgeCount = length(graph.edges)
        oldVariableEdges = deepcopy(graph.variableEdges)
        oldFactorEdges = deepcopy(graph.factorEdges)

        @test_throws ErrorException addFactor!(
            graph,
            GaussianFactor(:x1, 0.2, 1.0, 0.4; label = "prior_x1")
        )
        @test length(graph.factors) == oldFactorCount
        @test length(graph.edges) == oldEdgeCount
        @test graph.variableEdges == oldVariableEdges
        @test graph.factorEdges == oldFactorEdges

        @test_throws ErrorException addFactor!(
            graph,
            GaussianFactor(:unknown, 0.2, 1.0, 0.4; label = "unknown_prior")
        )
        @test length(graph.factors) == oldFactorCount
        @test length(graph.edges) == oldEdgeCount
        @test graph.variableEdges == oldVariableEdges
        @test graph.factorEdges == oldFactorEdges

        @test_throws ErrorException addFactor!(
            graph,
            GaussianFactor(:x1, [0.2, 0.3], 1.0, [0.4, 0.4]; label = "bad_dimensions")
        )
        @test length(graph.factors) == oldFactorCount
        @test length(graph.edges) == oldEdgeCount
        @test graph.variableEdges == oldVariableEdges
        @test graph.factorEdges == oldFactorEdges
    end
end

@testset verbose = true "Coefficient block lookups" begin
    graph = gaussianTreeTestGraph()

    @test coefficientBlock(graph; factor = "factor_x1_x2", variable = :x2) == [
        0.5  -0.3
        1.2   0.8
    ]
    @test coefficientBlock(graph; factor = "factor_x1_x2", variable = :x1) ==
        reshape([1.0, -0.4], 2, 1)
    @test coefficientBlocks(graph, "factor_x1_x2") == [
        reshape([1.0, -0.4], 2, 1),
        [
            0.5  -0.3
            1.2   0.8
        ]
    ]
    @test coefficientBlocks(graph, 2) == coefficientBlocks(
        graph,
        graph.factors[2]
    )
end

@testset verbose = true "GaussianFactor updates" begin
    @testset "Updates mean and covariance by label" begin
        graph = gaussianTreeTestGraph()
        factorIdx = factorIndex(graph, "prior_x1")
        oldFactor = graph.factors[factorIdx]
        oldWLS = solveWLS(graph)

        updated = updateFactor!(
            graph,
            "prior_x1";
            mean = [0.75],
            covariance = [0.05]
        )
        newWLS = solveWLS(graph)

        @test updated === graph.factors[factorIdx]
        @test updated.id == oldFactor.id
        @test updated.variables == oldFactor.variables
        @test updated.coefficient == oldFactor.coefficient
        @test updated.label == oldFactor.label
        @test updated.mean == [0.75]
        @test updated.covariance == reshape([0.05], 1, 1)
        @test newWLS.variableMean[1] != oldWLS.variableMean[1]
    end

    @testset "Updates mean by index" begin
        graph = gaussianTreeTestGraph()

        updated = updateFactor!(
            graph,
            1;
            mean = 0.5
        )

        @test updated.mean == [0.5]
        @test updated.covariance == reshape([0.15], 1, 1)
    end

    @testset "Updates by GaussianFactor keyword" begin
        graph = gaussianTreeTestGraph()

        updated = updateFactor!(
            graph;
            factor = "prior_x1",
            mean = 0.6,
            covariance = 0.2
        )

        @test updated === graph.factors[factorIndex(graph, "prior_x1")]
        @test updated.mean == [0.6]
        @test updated.covariance == reshape([0.2], 1, 1)
    end

    @testset "Initializes moment messages during warm-start update" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph)

        updated = updateFactor!(
            graph,
            inference;
            factor = "prior_x1",
            mean = 0.6,
            covariance = 0.2,
            initialize = true
        )

        variableIdx = variableIndex(graph, :x1)
        edgeIds = graph.variableEdges[variableIdx]

        @test updated.initialize
        @test isapprox(inference.initial[variableIdx].mean, [0.6])
        @test isapprox(inference.initial[variableIdx].covariance, reshape([0.2], 1, 1))
        @test isapprox(inference.marginal[variableIdx].mean, [0.6])
        @test isapprox(inference.marginal[variableIdx].covariance, reshape([0.2], 1, 1))

        for edgeId in edgeIds
            @test isapprox(inference.variableToFactor[edgeId].mean, [0.6])
            @test isapprox(
                inference.variableToFactor[edgeId].covariance,
                reshape([0.2], 1, 1)
            )
            @test isapprox(inference.factorToVariable[edgeId].mean, [0.6])
            @test isapprox(
                inference.factorToVariable[edgeId].covariance,
                reshape([0.2], 1, 1)
            )
        end
    end

    @testset "Initializes canonical messages during warm-start update" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(graph)

        updated = updateFactor!(
            graph,
            inference;
            factor = "prior_x1",
            mean = 0.6,
            covariance = 0.2,
            initialize = true
        )

        variableIdx = variableIndex(graph, :x1)
        edgeIds = graph.variableEdges[variableIdx]

        @test updated.initialize
        @test isapprox(inference.initial[variableIdx].information, [3.0])
        @test isapprox(inference.initial[variableIdx].precision, reshape([5.0], 1, 1))
        @test isapprox(inference.marginal[variableIdx].information, [3.0])
        @test isapprox(inference.marginal[variableIdx].precision, reshape([5.0], 1, 1))

        for edgeId in edgeIds
            @test isapprox(inference.variableToFactor[edgeId].information, [3.0])
            @test isapprox(
                inference.variableToFactor[edgeId].precision,
                reshape([5.0], 1, 1)
            )
            @test isapprox(inference.factorToVariable[edgeId].information, [3.0])
            @test isapprox(
                inference.factorToVariable[edgeId].precision,
                reshape([5.0], 1, 1)
            )
        end
    end

    @testset "Updates coefficient with fixed dimensions" begin
        graph = gaussianTreeTestGraph()
        factorIdx = factorIndex(graph, "factor_x1_x2")
        oldFactor = graph.factors[factorIdx]
        oldWLS = solveWLS(graph)
        newCoefficient = [0.9 0.4 -0.1; -0.2 1.0 0.6]

        updated = updateFactor!(
            graph;
            factor = "factor_x1_x2",
            coefficient = newCoefficient
        )
        newWLS = solveWLS(graph)

        @test updated === graph.factors[factorIdx]
        @test updated.id == oldFactor.id
        @test updated.variables == oldFactor.variables
        @test updated.mean == oldFactor.mean
        @test updated.covariance == oldFactor.covariance
        @test updated.label == oldFactor.label
        @test updated.coefficient == newCoefficient
        @test newWLS.variableMean != oldWLS.variableMean
    end

    @testset "Rejects invalid updates atomically" begin
        graph = gaussianTreeTestGraph()
        factorIdx = factorIndex(graph, "prior_x1")
        original = graph.factors[factorIdx]

        @test_throws ErrorException updateFactor!(
            graph,
            "prior_x1";
            mean = [0.4, 0.5]
        )
        @test graph.factors[factorIdx] === original

        @test_throws ErrorException updateFactor!(
            graph,
            "prior_x1";
            covariance = [0.0]
        )
        @test graph.factors[factorIdx] === original

        linkFactorIdx = factorIndex(graph, "factor_x1_x2")
        originalLink = graph.factors[linkFactorIdx]

        @test_throws ErrorException updateFactor!(
            graph,
            "factor_x1_x2";
            coefficient = [1.0 0.0]
        )
        @test graph.factors[linkFactorIdx] === originalLink

        @test_throws ErrorException updateFactor!(
            graph,
            "prior_x1";
            mean = [0.4, 0.5],
            coefficient = [1.0 0.0; 0.0 1.0],
            covariance = [0.4, 0.4]
        )
        @test graph.factors[factorIdx] === original

        @test_throws ErrorException updateFactor!(
            graph;
            factor = "factor_x1_x2",
            initialize = true
        )
        @test graph.factors[linkFactorIdx] === originalLink
    end
end
