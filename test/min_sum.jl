include("setup.jl")

function assertMinSumMatchesWLS(graph::GaussianFactorGraph, inference::GaussianMinSumInference)
    expected = solveWLS(graph)

    for variableIndexValue in eachindex(graph.variables)
        variable = graph.variables[variableIndexValue]

        @test estimate(graph, inference, variable.id) ≈
              expected.variableMean[variableIndexValue] atol = 1e-5
    end

    return nothing
end

function minSumMessageSnapshot(messages, edgeIds)
    return [
        (
            J = copy(messages[edgeId].J),
            h = copy(messages[edgeId].h),
            c = messages[edgeId].c
        )
        for edgeId in edgeIds
    ]
end

function assertMinSumMessagesUnchanged(messages, edgeIds, snapshots)
    for (snapshotIndex, edgeId) in pairs(edgeIds)
        @test messages[edgeId].J == snapshots[snapshotIndex].J
        @test messages[edgeId].h == snapshots[snapshotIndex].h
        @test messages[edgeId].c == snapshots[snapshotIndex].c
    end

    return nothing
end

@testset "Gaussian min-sum inference" begin
    @testset "Initial beliefs warm-start quadratic messages" begin
        variables = [
            GaussianVariable(:x1, 1; label = "x1", mean = 9.0, covariance = 3.0),
            GaussianVariable(:x2, 1; label = "x2")
        ]
        factors = [
            GaussianFactor(:x1, 5.0, 2.0, 8.0; label = "init_x1", initialize = true),
            GaussianFactor(:x1, :x2, 0.0, [1.0 -1.0], 1.0; label = "link")
        ]
        graph = factorGraph(variables, factors)
        inference = minsum(graph; mean = 2.0, covariance = 4.0)

        x1Message = QuadraticMessage([0.5;;], [1.25], 0.0)
        x2Message = QuadraticMessage([0.25;;], [0.5], 0.0)

        for edgeId in graph.variableEdges[variableIndex(graph, :x1)]
            @test inference.variableToFactor[edgeId].J ≈ x1Message.J
            @test inference.variableToFactor[edgeId].h ≈ x1Message.h
            @test inference.factorToVariable[edgeId].J ≈ x1Message.J
            @test inference.factorToVariable[edgeId].h ≈ x1Message.h
        end

        for edgeId in graph.variableEdges[variableIndex(graph, :x2)]
            @test inference.variableToFactor[edgeId].J ≈ x2Message.J
            @test inference.variableToFactor[edgeId].h ≈ x2Message.h
            @test inference.factorToVariable[edgeId].J ≈ x2Message.J
            @test inference.factorToVariable[edgeId].h ≈ x2Message.h
        end

        @test estimate(graph, inference, :x1) == [2.5]
        @test estimate(graph, inference, :x2) == [2.0]
    end

    @testset "Added min-sum GaussianVariable accepts components" begin
        graph = GaussianFactorGraph()
        inference = minsum(graph; mean = [1.0, -1.0], covariance = [2.0, 3.0])

        addVariable!(graph, inference, :x4, 2; label = "x4", components = [:re, :im])

        @test graph.variables[end].components == [:re, :im]
        @test componentIndex(graph, :x4, :re) == 1
        @test componentValue(graph, :x4, 2) == :im
        @test inference.estimate[end] == [1.0, -1.0]
    end

    @testset "Sequential GBP matches weighted least squares on a tree" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        gbp!(graph, inference; iterations = 10)

        assertMinSumMatchesWLS(graph, inference)
    end

    @testset "Flooding GBP matches weighted least squares on a tree" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        gbp!(graph, inference; iterations = 10, schedule = :flooding)

        assertMinSumMatchesWLS(graph, inference)
    end

    @testset "Broadcast sequential GBP matches weighted least squares on a tree" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        gbp!(graph, inference; iterations = 10, broadcast = true)

        assertMinSumMatchesWLS(graph, inference)
    end

    @testset "Broadcast flooding GBP matches weighted least squares on a tree" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        gbp!(graph, inference; iterations = 10, schedule = :flooding, broadcast = true)

        assertMinSumMatchesWLS(graph, inference)
    end

    @testset "Low-level schedule objects match weighted least squares on a tree" begin
        for scheduleFactory in (sequentialSchedule, floodingSchedule)
            graph = gaussianTreeTestGraph()
            inference = minsum(graph)
            schedule = scheduleFactory(graph, inference)

            for _ in 1:10
                messages!(graph, inference, schedule)
                estimates!(graph, inference)
            end

            assertMinSumMatchesWLS(graph, inference)
        end
    end

    @testset "Damping matches weighted least squares on a tree" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        dampEdges!(graph, inference; variable = :x2, factor = "factor_x1_x2")
        @test areDampedEdges(graph, inference; variable = :x2, factor = "factor_x1_x2")

        gbp!(graph, inference; iterations = 120, damping = true, prob = 1.0, alpha = 0.25)
        assertMinSumMatchesWLS(graph, inference)

        undampEdges!(graph, inference; variable = :x2, factor = "factor_x1_x2")
        @test !areDampedEdges(graph, inference; variable = :x2, factor = "factor_x1_x2")
    end

    @testset "Freezing preserves selected outgoing messages" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        gbp!(graph, inference; iterations = 2)

        factorIdx = factorIndex(graph, "factor_x1_x2")
        edgeIds = graph.factorEdges[factorIdx]
        snapshots = minSumMessageSnapshot(inference.factorToVariable, edgeIds)

        freezeFactor!(graph, inference, "factor_x1_x2")
        @test isFrozenFactor(graph, inference, "factor_x1_x2")

        gbp!(graph, inference; iterations = 3)
        assertMinSumMessagesUnchanged(inference.factorToVariable, edgeIds, snapshots)

        unfreezeFactor!(graph, inference, "factor_x1_x2")
        @test !isFrozenFactor(graph, inference, "factor_x1_x2")
    end

    @testset "Residual schedule updates requested count" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)
        schedule = residualSchedule(graph, inference; updateCount = 2)

        messages!(graph, inference, schedule)
        estimates!(graph, inference)

        @test length(schedule.lastUpdated) == 2
        @test maximum(schedule.residuals) > 0.0
    end

    @testset "Residual GBP matches weighted least squares on a tree" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        gbp!(
            graph,
            inference;
            iterations = 20,
            schedule = :residual,
            updateFraction = 1.0,
            tolerance = 1e-10
        )

        assertMinSumMatchesWLS(graph, inference)
    end

    @testset "Tree forward-backward matches weighted least squares" begin
        tree = treeFactorGraph(gaussianTreeTestGraph(); root = :x1)
        inference = minsum(tree)

        forwardBackward!(tree, inference)

        assertMinSumMatchesWLS(tree.graph, inference)
    end

    @testset "Tree forward and backward steps support min-sum" begin
        tree = treeFactorGraph(gaussianTreeTestGraph(); root = :x1)
        inference = minsum(tree)
        schedule = forwardBackwardSchedule(tree)

        for edgeId in tree.forwardOrder
            @test forwardStep!(tree, inference, schedule) == edgeId
        end

        @test forwardStep!(tree, inference, schedule) === nothing

        for edgeId in tree.backwardOrder
            @test backwardStep!(tree, inference, schedule) == edgeId
        end

        @test backwardStep!(tree, inference, schedule) === nothing

        estimates!(tree, inference)
        assertMinSumMatchesWLS(tree.graph, inference)
    end

    @testset "Dynamic graph updates extend the inference state" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                variables = [
                    GaussianVariable(:x1, 1; label = "x1"),
                    GaussianVariable(:x2, 1; label = "x2")
                ]
                factors = [
                    GaussianFactor(:x1, 1.0, 1.0, 0.1; label = "prior_x1"),
                    GaussianFactor(:x1, :x2, 0.5, [1.0 -1.0], 0.2; label = "link")
                ]
                graph = factorGraph(variables, factors)
                inference = minsum(graph)

                addFactor!(graph, inference, :x2, 1.2, 1.0, 0.3; label = "prior_x2")
                updateFactor!(graph, inference; factor = "prior_x2", mean = 1.1)

                gbp!(
                    graph,
                    inference;
                    iterations = 10,
                    schedule = gbpSchedule(schedule.flooding),
                    broadcast = schedule.broadcast
                )

                assertMinSumMatchesWLS(graph, inference)
            end
        end
    end
end
