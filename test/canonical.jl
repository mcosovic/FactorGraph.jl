include("setup.jl")

function canonicalTestGraph()
    return gaussianTestGraph()
end

function canonicalInference(
    graph::GaussianFactorGraph;
    broadcast::Bool,
    flooding::Bool
)
    return canonical(
        graph;
        mean = 0.0,
        covariance = 1e6,
    )
end

function assertCanonicalMarginalsAccurate(
    graph::GaussianFactorGraph,
    inference::GaussianCanonicalInference;
    meanAtol::Float64 = 1e-5,
    covarianceAtol::Float64 = 1e-10
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

        @test marginalMean(graph, inference, variable.id) == marginalData.mean
        @test marginalCovariance(graph, inference, variable.id) == marginalData.covariance
        @test size(marginalData.covariance) == (variable.dimension, variable.dimension)
        @test isapprox(
            marginalData.covariance,
            marginalData.covariance';
            atol = covarianceAtol
        )
        @test isposdef(Symmetric(marginalData.covariance))
    end

    return nothing
end

function assertCanonicalMarginalsAccurate(;
    broadcast::Bool,
    flooding::Bool,
    iterations::Int = 60,
    meanAtol::Float64 = 1e-5,
    covarianceAtol::Float64 = 1e-10
)
    graph = canonicalTestGraph()
    inference = canonicalInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(
        graph,
        inference;
        iterations = iterations,
        schedule = gbpSchedule(flooding),
        broadcast = broadcast
    )
    assertCanonicalMarginalsAccurate(
        graph,
        inference;
        meanAtol = meanAtol,
        covarianceAtol = covarianceAtol
    )

    return nothing
end

function messageSnapshots(messages, edgeIds)
    return [
        (
            information = copy(messages[edgeId].information),
            precision = copy(messages[edgeId].precision)
        )
        for edgeId in edgeIds
    ]
end

function assertMessagesUnchanged(messages, edgeIds, snapshots)
    for (snapshotIndex, edgeId) in pairs(edgeIds)
        @test messages[edgeId].information == snapshots[snapshotIndex].information
        @test messages[edgeId].precision == snapshots[snapshotIndex].precision
    end

    return nothing
end

function assertFreezeFactor(;
    broadcast::Bool,
    flooding::Bool,
    factor = "factor_x1_x2"
)
    graph = canonicalTestGraph()
    inference = canonicalInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(graph, inference; iterations = 2, schedule = gbpSchedule(flooding), broadcast = broadcast)

    factorIdx = factorIndex(graph, factor)
    edgeIds = graph.factorEdges[factorIdx]
    snapshots = messageSnapshots(inference.factorToVariable, edgeIds)

    freezeFactor!(graph, inference, factor)
    @test isFrozenFactor(graph, inference, factor)

    gbp!(graph, inference; iterations = 3, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMessagesUnchanged(inference.factorToVariable, edgeIds, snapshots)

    unfreezeFactor!(graph, inference, factor)
    @test !isFrozenFactor(graph, inference, factor)

    gbp!(graph, inference; iterations = 60, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertCanonicalMarginalsAccurate(graph, inference)

    return nothing
end

function assertFreezeVariable(;
    broadcast::Bool,
    flooding::Bool,
    variable::VariableRef = :x2
)
    graph = canonicalTestGraph()
    inference = canonicalInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(graph, inference; iterations = 2, schedule = gbpSchedule(flooding), broadcast = broadcast)

    variableIdx = variableIndex(graph, variable)
    edgeIds = graph.variableEdges[variableIdx]
    snapshots = messageSnapshots(inference.variableToFactor, edgeIds)

    freezeVariable!(graph, inference, variable)
    @test isFrozenVariable(graph, inference, variable)

    gbp!(graph, inference; iterations = 3, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMessagesUnchanged(inference.variableToFactor, edgeIds, snapshots)

    unfreezeVariable!(graph, inference, variable)
    @test !isFrozenVariable(graph, inference, variable)

    gbp!(graph, inference; iterations = 60, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertCanonicalMarginalsAccurate(graph, inference)

    return nothing
end

function assertFreezeEdge(;
    broadcast::Bool,
    flooding::Bool,
    variable::VariableRef = :x2,
    factor = "factor_x1_x2"
)
    graph = canonicalTestGraph()
    inference = canonicalInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(graph, inference; iterations = 2, schedule = gbpSchedule(flooding), broadcast = broadcast)

    edgeId = edgeIndex(graph; variable = variable, factor = factor)
    variableToFactorSnapshot = messageSnapshots(inference.variableToFactor, [edgeId])
    factorToVariableSnapshot = messageSnapshots(inference.factorToVariable, [edgeId])

    freezeEdge!(graph, inference; variable = variable, factor = factor)
    @test isFrozenEdge(graph, inference; variable = variable, factor = factor)

    gbp!(graph, inference; iterations = 3, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMessagesUnchanged(
        inference.variableToFactor,
        [edgeId],
        variableToFactorSnapshot
    )
    assertMessagesUnchanged(
        inference.factorToVariable,
        [edgeId],
        factorToVariableSnapshot
    )

    unfreezeEdge!(graph, inference; variable = variable, factor = factor)
    @test !isFrozenEdge(graph, inference; variable = variable, factor = factor)

    gbp!(graph, inference; iterations = 60, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertCanonicalMarginalsAccurate(graph, inference)

    return nothing
end

function assertGlobalDamping(;
    broadcast::Bool,
    flooding::Bool
)
    graph = canonicalTestGraph()
    inference = canonicalInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(
        graph,
        inference;
        iterations = 120,
        schedule = gbpSchedule(flooding),
        broadcast = broadcast,
        damping = true,
        prob = 1.0,
        alpha = 0.35
    )
    assertCanonicalMarginalsAccurate(graph, inference)

    return nothing
end

function assertDampedEdge(;
    broadcast::Bool,
    flooding::Bool,
    variable::VariableRef = :x2,
    factor = "factor_x1_x2"
)
    graph = canonicalTestGraph()
    inference = canonicalInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    @test !isDampedEdge(
        graph,
        inference;
        variable = variable,
        factor = factor
    )

    dampEdges!(
        graph,
        inference;
        variable = variable,
        factor = factor,
        prob = 1.0,
        alpha = 0.35
    )
    @test isDampedEdge(
        graph,
        inference;
        variable = variable,
        factor = factor
    )

    undampEdges!(graph, inference)
    @test !any(inference.dampedEdges)

    dampEdges!(graph, inference; prob = 1.0, alpha = 0.25)
    @test all(inference.dampedEdges)
    @test all(==(0.25), inference.edgeDampingAlpha)

    gbp!(
        graph,
        inference;
        iterations = 120,
        schedule = gbpSchedule(flooding),
        broadcast = broadcast
    )
    assertCanonicalMarginalsAccurate(graph, inference)

    undampEdges!(
        graph,
        inference;
        variable = variable,
        factor = factor
    )
    @test !isDampedEdge(
        graph,
        inference;
        variable = variable,
        factor = factor
    )

    return nothing
end

@testset verbose = true "Canonical form" begin
    for schedule in SCHEDULE_CASES
        @testset "Schedule: $(schedule.name)" begin
            assertCanonicalMarginalsAccurate(
                broadcast = schedule.broadcast,
                flooding = schedule.flooding
            )
        end
    end

    @testset "Freeze: factor" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertFreezeFactor(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Freeze: variable" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertFreezeVariable(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Freeze: edge" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertFreezeEdge(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Damping: global" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertGlobalDamping(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Damping: edge" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertDampedEdge(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end
end
