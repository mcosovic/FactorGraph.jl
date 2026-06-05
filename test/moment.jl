include("setup.jl")

function momentTestGraph()
    return gaussianTestGraph()
end

function momentInference(
    graph::GaussianFactorGraph;
    broadcast::Bool,
    flooding::Bool
)
    return moment(
        graph;
        mean = 0.0,
        covariance = 1e6,
    )
end

function assertMomentMarginalsAccurate(
    graph::GaussianFactorGraph,
    inference::GaussianMomentInference;
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

function assertMomentMarginalsAccurate(;
    broadcast::Bool,
    flooding::Bool,
    iterations::Int = 60,
    meanAtol::Float64 = 1e-5,
    covarianceAtol::Float64 = 1e-10
)
    graph = momentTestGraph()
    inference = momentInference(
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
    assertMomentMarginalsAccurate(
        graph,
        inference;
        meanAtol = meanAtol,
        covarianceAtol = covarianceAtol
    )

    return nothing
end

function momentMessageSnapshots(messages, edgeIds)
    return [
        (
            mean = copy(messages[edgeId].mean),
            covariance = copy(messages[edgeId].covariance)
        )
        for edgeId in edgeIds
    ]
end

function assertMomentMessagesUnchanged(messages, edgeIds, snapshots)
    for (snapshotIndex, edgeId) in pairs(edgeIds)
        @test messages[edgeId].mean == snapshots[snapshotIndex].mean
        @test messages[edgeId].covariance == snapshots[snapshotIndex].covariance
    end

    return nothing
end

function assertMomentFreezeFactor(;
    broadcast::Bool,
    flooding::Bool,
    factor = "factor_x1_x2"
)
    graph = momentTestGraph()
    inference = momentInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(graph, inference; iterations = 2, schedule = gbpSchedule(flooding), broadcast = broadcast)

    factorIdx = factorIndex(graph, factor)
    edgeIds = graph.factorEdges[factorIdx]
    snapshots = momentMessageSnapshots(inference.factorToVariable, edgeIds)

    freezeFactor!(graph, inference, factor)
    @test isFrozenFactor(graph, inference, factor)

    gbp!(graph, inference; iterations = 3, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMomentMessagesUnchanged(inference.factorToVariable, edgeIds, snapshots)

    unfreezeFactor!(graph, inference, factor)
    @test !isFrozenFactor(graph, inference, factor)

    gbp!(graph, inference; iterations = 60, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMomentMarginalsAccurate(graph, inference)

    return nothing
end

function assertMomentFreezeVariable(;
    broadcast::Bool,
    flooding::Bool,
    variable::VariableRef = :x2
)
    graph = momentTestGraph()
    inference = momentInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(graph, inference; iterations = 2, schedule = gbpSchedule(flooding), broadcast = broadcast)

    variableIdx = variableIndex(graph, variable)
    edgeIds = graph.variableEdges[variableIdx]
    snapshots = momentMessageSnapshots(inference.variableToFactor, edgeIds)

    freezeVariable!(graph, inference, variable)
    @test isFrozenVariable(graph, inference, variable)

    gbp!(graph, inference; iterations = 3, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMomentMessagesUnchanged(inference.variableToFactor, edgeIds, snapshots)

    unfreezeVariable!(graph, inference, variable)
    @test !isFrozenVariable(graph, inference, variable)

    gbp!(graph, inference; iterations = 60, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMomentMarginalsAccurate(graph, inference)

    return nothing
end

function assertMomentFreezeEdge(;
    broadcast::Bool,
    flooding::Bool,
    variable::VariableRef = :x2,
    factor = "factor_x1_x2"
)
    graph = momentTestGraph()
    inference = momentInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    gbp!(graph, inference; iterations = 2, schedule = gbpSchedule(flooding), broadcast = broadcast)

    edgeId = edgeIndex(graph; variable = variable, factor = factor)
    variableToFactorSnapshot = momentMessageSnapshots(inference.variableToFactor, [edgeId])
    factorToVariableSnapshot = momentMessageSnapshots(inference.factorToVariable, [edgeId])

    freezeEdge!(graph, inference; variable = variable, factor = factor)
    @test isFrozenEdge(graph, inference; variable = variable, factor = factor)

    gbp!(graph, inference; iterations = 3, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMomentMessagesUnchanged(
        inference.variableToFactor,
        [edgeId],
        variableToFactorSnapshot
    )
    assertMomentMessagesUnchanged(
        inference.factorToVariable,
        [edgeId],
        factorToVariableSnapshot
    )

    unfreezeEdge!(graph, inference; variable = variable, factor = factor)
    @test !isFrozenEdge(graph, inference; variable = variable, factor = factor)

    gbp!(graph, inference; iterations = 60, schedule = gbpSchedule(flooding), broadcast = broadcast)
    assertMomentMarginalsAccurate(graph, inference)

    return nothing
end

function assertMomentGlobalDamping(;
    broadcast::Bool,
    flooding::Bool
)
    graph = momentTestGraph()
    inference = momentInference(
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
    assertMomentMarginalsAccurate(graph, inference)

    return nothing
end

function assertMomentDampedEdge(;
    broadcast::Bool,
    flooding::Bool,
    variable::VariableRef = :x2,
    factor = "factor_x1_x2"
)
    graph = momentTestGraph()
    inference = momentInference(
        graph;
        broadcast = broadcast,
        flooding = flooding
    )

    @test !areDampedEdges(
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
    @test areDampedEdges(
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
    assertMomentMarginalsAccurate(graph, inference)

    undampEdges!(
        graph,
        inference;
        variable = variable,
        factor = factor
    )
    @test !areDampedEdges(
        graph,
        inference;
        variable = variable,
        factor = factor
    )

    return nothing
end

@testset verbose = true "Moment form" begin
    for schedule in SCHEDULE_CASES
        @testset "Schedule: $(schedule.name)" begin
            assertMomentMarginalsAccurate(
                broadcast = schedule.broadcast,
                flooding = schedule.flooding
            )
        end
    end

    @testset "Freeze: factor" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertMomentFreezeFactor(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Freeze: variable" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertMomentFreezeVariable(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Freeze: edge" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertMomentFreezeEdge(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Damping: global" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertMomentGlobalDamping(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end

    @testset "Damping: edge" begin
        for schedule in SCHEDULE_CASES
            @testset "schedule=$(schedule.name)" begin
                assertMomentDampedEdge(
                    broadcast = schedule.broadcast,
                    flooding = schedule.flooding
                )
            end
        end
    end
end
