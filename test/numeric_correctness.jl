include("setup.jl")

using Random

function numericGaussianReferenceGraph()
    variables = [
        GaussianVariable(:x1, 1; label = "x1"),
        GaussianVariable(:x2, 1; label = "x2")
    ]
    factors = [
        GaussianFactor(:x1, 1.0, 1.0, 0.25; label = "prior_x1"),
        GaussianFactor(:x2, 3.0, 1.0, 0.25; label = "prior_x2"),
        GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.5; label = "difference")
    ]

    return factorGraph(variables, factors)
end

function numericGaussianVectorGraph()
    variables = [
        GaussianVariable(:x1, 2; label = "x1"),
        GaussianVariable(:x2, 2; label = "x2")
    ]
    factors = [
        GaussianFactor(
            :x1,
            [1.0, -1.0],
            Matrix{Float64}(I, 2, 2),
            [0.5, 0.25];
            label = "prior_x1"
        ),
        GaussianFactor(
            :x2,
            [2.0, 0.5],
            Matrix{Float64}(I, 2, 2),
            [0.4, 0.3];
            label = "prior_x2"
        ),
        GaussianFactor(
            :x1, :x2,
            [1.2, -0.3],
            [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0],
            [0.2, 0.35];
            label = "difference"
        )
    ]

    return factorGraph(variables, factors)
end

function numericGaussianVectorChainGraph()
    variables = [
        GaussianVariable(:x1, 2; label = "x1", components = [:p, :v]),
        GaussianVariable(:x2, 2; label = "x2", components = [:p, :v]),
        GaussianVariable(:x3, 2; label = "x3", components = [:p, :v])
    ]
    factors = [
        GaussianFactor(
            :x1,
            [0.5, -0.2],
            Matrix{Float64}(I, 2, 2),
            [0.3 0.08; 0.08 0.4];
            label = "prior_x1"
        ),
        GaussianFactor(
            :x3,
            [2.2, 0.6],
            Matrix{Float64}(I, 2, 2),
            [0.45 -0.06; -0.06 0.35];
            label = "prior_x3"
        ),
        GaussianFactor(
            :x1, :x2,
            [0.9, 0.1],
            [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0],
            [0.25 0.04; 0.04 0.3];
            label = "x1_x2"
        ),
        GaussianFactor(
            :x2, :x3,
            [1.1, 0.2],
            [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0],
            [0.2 -0.03; -0.03 0.28];
            label = "x2_x3"
        )
    ]

    return factorGraph(variables, factors)
end

function numericGaussianLoopGraph()
    variables = [
        GaussianVariable(:x1, 1; label = "x1"),
        GaussianVariable(:x2, 1; label = "x2"),
        GaussianVariable(:x3, 1; label = "x3")
    ]
    factors = [
        GaussianFactor(:x1, 0.0, 1.0, 0.4; label = "prior_x1"),
        GaussianFactor(:x2, 1.0, 1.0, 0.5; label = "prior_x2"),
        GaussianFactor(:x3, 2.0, 1.0, 0.6; label = "prior_x3"),
        GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.25; label = "x1_x2"),
        GaussianFactor(:x2, :x3, 1.0, [-1.0 1.0], 0.3; label = "x2_x3"),
        GaussianFactor(:x1, :x3, 2.1, [-1.0 1.0], 0.35; label = "x1_x3")
    ]

    return factorGraph(variables, factors)
end

function numericDiscreteReferenceGraph()
    variables = [
        DiscreteVariable(:x1, 3; label = "x1", states = [:a, :b, :c]),
        DiscreteVariable(:x2, 2; label = "x2", states = [:off, :on]),
        DiscreteVariable(:x3, 2; label = "x3", states = [:low, :high])
    ]
    factors = [
        DiscreteFactor(:x1, [0.2, 0.5, 0.3]; label = "prior_x1"),
        DiscreteFactor(:x1, :x2, [0.9 0.1; 0.4 0.8; 0.2 0.7]; label = "x1_x2"),
        DiscreteFactor(:x2, :x3, [0.3 0.7; 0.6 0.4]; label = "x2_x3")
    ]

    return factorGraph(variables, factors)
end

function numericDiscreteLoopGraph()
    variables = [
        DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
        DiscreteVariable(:x2, 2; label = "x2", states = [:off, :on]),
        DiscreteVariable(:x3, 2; label = "x3", states = [:off, :on])
    ]
    factors = [
        DiscreteFactor(:x1, [0.55, 0.45]; label = "prior_x1"),
        DiscreteFactor(:x2, [0.35, 0.65]; label = "prior_x2"),
        DiscreteFactor(:x3, [0.6, 0.4]; label = "prior_x3"),
        DiscreteFactor(:x1, :x2, [1.2 0.8; 0.7 1.4]; label = "x1_x2"),
        DiscreteFactor(:x2, :x3, [1.1 0.9; 0.6 1.3]; label = "x2_x3"),
        DiscreteFactor(:x1, :x3, [1.3 0.75; 0.85 1.15]; label = "x1_x3")
    ]

    return factorGraph(variables, factors)
end

function numericDiscreteTernaryFactorGraph()
    variables = [
        DiscreteVariable(:x1, 2; label = "x1", states = [:a, :b]),
        DiscreteVariable(:x2, 3; label = "x2", states = [:low, :mid, :high]),
        DiscreteVariable(:x3, 2; label = "x3", states = [:off, :on])
    ]
    table = Array{Float64}(undef, 2, 3, 2)
    table[:, :, 1] .= [0.9 0.2 0.5; 0.3 0.8 0.4]
    table[:, :, 2] .= [0.4 0.7 0.1; 0.6 0.2 0.9]
    factors = [
        DiscreteFactor(:x1, [0.45, 0.55]; label = "prior_x1"),
        DiscreteFactor(:x2, [0.2, 0.5, 0.3]; label = "prior_x2"),
        DiscreteFactor(:x3, [0.65, 0.35]; label = "prior_x3"),
        DiscreteFactor(:x1, :x2, :x3, table; label = "ternary")
    ]

    return factorGraph(variables, factors)
end

function normalEquationReference(graph::GaussianFactorGraph)
    starts = Int[]
    stops = Int[]
    start = 1

    for variable in graph.variables
        stop = start + variable.dimension - 1

        push!(starts, start)
        push!(stops, stop)

        start = stop + 1
    end

    dimension = isempty(stops) ? 0 : last(stops)
    precision = zeros(dimension, dimension)
    information = zeros(dimension)

    for factorData in graph.factors
        H = zeros(length(factorData.mean), dimension)
        localStart = 1

        for variableRef in factorData.variables
            variableIdx = variableIndex(graph, variableRef)
            variableDimension = graph.variables[variableIdx].dimension
            localStop = localStart + variableDimension - 1

            H[:, starts[variableIdx]:stops[variableIdx]] .=
                factorData.coefficient[:, localStart:localStop]

            localStart = localStop + 1
        end

        factorPrecision = inv(Symmetric(factorData.covariance))

        precision .+= H' * factorPrecision * H
        information .+= H' * factorPrecision * factorData.mean
    end

    covariance = inv(Symmetric(precision))
    mean = covariance * information
    variableMean = [mean[starts[index]:stops[index]] for index in eachindex(graph.variables)]
    variableCovariance = [
        covariance[starts[index]:stops[index], starts[index]:stops[index]]
        for index in eachindex(graph.variables)
    ]

    return (; mean, covariance, variableMean, variableCovariance)
end

function exactDiscreteMarginalsAndMAP(graph::DiscreteFactorGraph)
    marginals = [zeros(variable.cardinality) for variable in graph.variables]
    bestWeight = -Inf
    bestAssignment = nothing
    total = 0.0
    ranges = [1:variable.cardinality for variable in graph.variables]

    for assignment in Iterators.product(ranges...)
        weight = 1.0

        for factorData in graph.factors
            indices = map(factorData.variables) do variableRef
                variableIdx = variableIndex(graph, variableRef)

                return assignment[variableIdx]
            end

            weight *= factorData.table[indices...]
        end

        total += weight

        for variableIdx in eachindex(marginals)
            marginals[variableIdx][assignment[variableIdx]] += weight
        end

        if weight > bestWeight
            bestWeight = weight
            bestAssignment = assignment
        end
    end

    for marginal in marginals
        marginal ./= total
    end

    mapEstimate = [
        graph.variables[index].states[bestAssignment[index]]
        for index in eachindex(graph.variables)
    ]

    return marginals, mapEstimate
end

function discreteAssignmentWeight(graph::DiscreteFactorGraph, states)
    assignment = [
        stateIndex(graph, graph.variables[index].id, states[index])
        for index in eachindex(graph.variables)
    ]
    weight = 1.0

    for factorData in graph.factors
        indices = map(factorData.variables) do variableRef
            variableIdx = variableIndex(graph, variableRef)

            return assignment[variableIdx]
        end

        weight *= factorData.table[indices...]
    end

    return weight
end

function bestDiscreteWeight(graph::DiscreteFactorGraph)
    bestWeight = -Inf
    ranges = [1:variable.cardinality for variable in graph.variables]

    for assignment in Iterators.product(ranges...)
        weight = 1.0

        for factorData in graph.factors
            indices = map(factorData.variables) do variableRef
                variableIdx = variableIndex(graph, variableRef)

                return assignment[variableIdx]
            end

            weight *= factorData.table[indices...]
        end

        bestWeight = max(bestWeight, weight)
    end

    return bestWeight
end

function assertGaussianMarginalsMatch(graph, inference, reference; atol = 1e-8)
    for variableIdx in eachindex(graph.variables)
        variable = graph.variables[variableIdx]

        @test marginalMean(graph, inference, variable.id) ≈
              reference.variableMean[variableIdx] atol = atol
        @test marginalCovariance(graph, inference, variable.id) ≈
              reference.variableCovariance[variableIdx] atol = atol
    end
end

function refreshGaussianSolution!(graph, inference)
    if inference isa GaussianMinSumInference
        estimates!(graph, inference)
    else
        marginals!(graph, inference)
    end

    return inference
end

function assertGaussianMeanSolution(graph, inference, reference; atol = 1e-8)
    for variableIdx in eachindex(graph.variables)
        variable = graph.variables[variableIdx]

        if inference isa GaussianMinSumInference
            @test estimate(graph, inference, variable.id) ≈
                  reference.variableMean[variableIdx] atol = atol
        else
            @test marginalMean(graph, inference, variable.id) ≈
                  reference.variableMean[variableIdx] atol = atol
        end
    end
end

function assertDiscreteMarginalsMatch(graph, inference, expectedMarginals; atol = 1e-8)
    for variableIdx in eachindex(graph.variables)
        variable = graph.variables[variableIdx]

        @test marginal(graph, inference, variable.id) ≈ expectedMarginals[variableIdx] atol = atol
    end

    return nothing
end

function assertDiscreteMAPMatch(graph, inference, expectedMAP)
    for variableIdx in eachindex(graph.variables)
        variable = graph.variables[variableIdx]

        @test estimate(graph, inference, variable.id) == expectedMAP[variableIdx]
    end

    return nothing
end

function storedMessageSnapshot(message)
    if message isa AbstractVector
        return (values = copy(message),)
    elseif hasproperty(message, :J)
        return (J = copy(message.J), h = copy(message.h), c = message.c)
    elseif hasproperty(message, :mean)
        return (
            mean = copy(message.mean),
            covariance = copy(message.covariance),
            precision = copy(message.precision),
            information = copy(message.information)
        )
    else
        return (
            information = copy(message.information),
            precision = copy(message.precision)
        )
    end
end

function storedMessageSnapshots(messages, edgeIds)
    return [storedMessageSnapshot(messages[edgeId]) for edgeId in edgeIds]
end

function assertStoredMessagesEqual(messages, edgeIds, snapshots)
    for (snapshotIdx, edgeId) in pairs(edgeIds)
        snapshot = snapshots[snapshotIdx]
        message = messages[edgeId]

        if hasproperty(snapshot, :values)
            @test message ≈ snapshot.values
        elseif hasproperty(snapshot, :J)
            @test message.J ≈ snapshot.J
            @test message.h ≈ snapshot.h
            @test message.c ≈ snapshot.c
        elseif hasproperty(snapshot, :mean)
            @test message.mean ≈ snapshot.mean
            @test message.covariance ≈ snapshot.covariance
            @test message.precision ≈ snapshot.precision
            @test message.information ≈ snapshot.information
        else
            @test message.information ≈ snapshot.information
            @test message.precision ≈ snapshot.precision
        end
    end

    return nothing
end

function assertAllFactorMessagesEqual(left, right)
    edgeIds = collect(eachindex(left.factorToVariable))
    snapshots = storedMessageSnapshots(right.factorToVariable, edgeIds)

    return assertStoredMessagesEqual(left.factorToVariable, edgeIds, snapshots)
end

function exerciseFreezeSchedule!(graph, inference, scheduleKind; damping::Bool = false)
    if scheduleKind === nothing
        return messages!(
            graph,
            inference;
            damping = damping,
            prob = 1.0,
            alpha = 0.25,
            rng = MersenneTwister(41)
        )
    elseif scheduleKind == :flooding
        return messages!(
            graph,
            inference;
            schedule = :flooding,
            damping = damping,
            prob = 1.0,
            alpha = 0.25,
            rng = MersenneTwister(41)
        )
    else
        return messages!(
            graph,
            inference;
            schedule = :residual,
            updateFraction = 1.0,
            damping = damping,
            prob = 1.0,
            alpha = 0.25,
            rng = MersenneTwister(41)
        )
    end
end

function runGaussianSchedule!(graph, inference, scheduleKind)
    if scheduleKind === nothing
        gbp!(graph, inference; iterations = 220, tolerance = 1e-12)
    elseif scheduleKind == :flooding
        gbp!(graph, inference; schedule = :flooding, iterations = 220, tolerance = 1e-12)
    else
        gbp!(
            graph,
            inference;
            schedule = :residual,
            updateFraction = 1.0,
            iterations = 220,
            tolerance = 1e-12
        )
    end

    return refreshGaussianSolution!(graph, inference)
end

function quietDisplayCall(f)
    displayBackend = TextDisplay(devnull)
    pushdisplay(displayBackend)

    try
        return redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                f()
            end
        end
    finally
        popdisplay(displayBackend)
    end
end

@testset verbose = true "Numerical correctness" begin
    @testset "Gaussian inference matches analytic posterior" begin
        graph = numericGaussianReferenceGraph()
        expectedMean = [[1.25], [2.75]]
        expectedFullMean = [1.25, 2.75]
        expectedFullCovariance = [0.1875 0.0625; 0.0625 0.1875]
        expectedCovariance = [reshape([0.1875], 1, 1), reshape([0.1875], 1, 1)]
        reference = solveWLS(graph)

        @test reference.mean ≈ expectedFullMean
        @test reference.covariance ≈ expectedFullCovariance
        @test reference.variableMean ≈ expectedMean
        @test reference.variableCovariance ≈ expectedCovariance

        for constructor in (moment, canonical)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            gbp!(graph, inference; iterations = 10, tolerance = 0.0)

            for variableIdx in eachindex(graph.variables)
                variable = graph.variables[variableIdx]

                @test marginalMean(graph, inference, variable.id) ≈ expectedMean[variableIdx]
                @test marginalCovariance(graph, inference, variable.id) ≈
                      expectedCovariance[variableIdx]
            end
        end

        minSumInference = minsum(graph)

        gbp!(graph, minSumInference; iterations = 10, tolerance = 0.0)

        @test estimate(graph, minSumInference, :x1) ≈ expectedMean[1]
        @test estimate(graph, minSumInference, :x2) ≈ expectedMean[2]
    end

    @testset "Vector Gaussian inference matches normal equations" begin
        graph = numericGaussianVectorGraph()
        reference = normalEquationReference(graph)
        wls = solveWLS(graph)

        @test wls.mean ≈ reference.mean
        @test wls.covariance ≈ reference.covariance
        @test wls.variableMean ≈ reference.variableMean
        @test wls.variableCovariance ≈ reference.variableCovariance

        for constructor in (moment, canonical)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            gbp!(graph, inference; iterations = 20, tolerance = 0.0)

            assertGaussianMarginalsMatch(graph, inference, reference)
        end

        minSumInference = minsum(graph)

        gbp!(graph, minSumInference; iterations = 20, tolerance = 0.0)

        for variableIdx in eachindex(graph.variables)
            variable = graph.variables[variableIdx]

            @test estimate(graph, minSumInference, variable.id) ≈
                  reference.variableMean[variableIdx]
        end
    end

    @testset "Correlated vector Gaussian chain matches WLS" begin
        graph = numericGaussianVectorChainGraph()
        reference = normalEquationReference(graph)
        wls = solveWLS(graph)

        @test wls.mean ≈ reference.mean
        @test wls.covariance ≈ reference.covariance
        @test wls.variableMean ≈ reference.variableMean
        @test wls.variableCovariance ≈ reference.variableCovariance

        for constructor in (moment, canonical)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            gbp!(graph, inference; iterations = 30, tolerance = 0.0)

            assertGaussianMarginalsMatch(graph, inference, reference)
        end

        minSumInference = minsum(graph)

        gbp!(graph, minSumInference; iterations = 30, tolerance = 0.0)
        assertGaussianMeanSolution(graph, minSumInference, reference)
    end

    @testset "Gaussian graph updates match fresh WLS solutions" begin
        for constructor in (moment, canonical, minsum)
            graph = numericGaussianReferenceGraph()
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            addVariable!(graph, inference, :x3, 1; label = "x3")
            addFactor!(graph, inference, :x2, :x3, 1.1, [-1.0 1.0], 0.4; label = "x2_x3")
            addFactor!(graph, inference, :x3, 3.4, 1.0, 0.6; label = "prior_x3")
            updateFactor!(graph, inference; factor = "prior_x2", mean = 2.6)

            reference = normalEquationReference(graph)

            gbp!(graph, inference; iterations = 200, tolerance = 0.0)

            assertGaussianMeanSolution(graph, inference, reference; atol = 1e-6)

            if inference isa GaussianSumProductInference
                assertGaussianMarginalsMatch(graph, inference, reference; atol = 1e-6)
            end
        end
    end

    @testset "Gaussian warm-start updates match fresh inference numerics" begin
        for constructor in (moment, canonical, minsum)
            warmGraph = numericGaussianVectorChainGraph()
            freshGraph = numericGaussianVectorChainGraph()
            warmInference = constructor(warmGraph; mean = 0.0, covariance = 100.0)

            gbp!(warmGraph, warmInference; iterations = 30, tolerance = 0.0)
            updateFactor!(
                warmGraph,
                warmInference,
                "x2_x3";
                mean = [0.8, 0.35],
                covariance = [0.22 -0.05; -0.05 0.31]
            )
            updateFactor!(
                freshGraph,
                "x2_x3";
                mean = [0.8, 0.35],
                covariance = [0.22 -0.05; -0.05 0.31]
            )

            freshInference = constructor(freshGraph; mean = 0.0, covariance = 100.0)
            reference = normalEquationReference(freshGraph)

            gbp!(warmGraph, warmInference; iterations = 40, tolerance = 0.0)
            gbp!(freshGraph, freshInference; iterations = 40, tolerance = 0.0)

            assertGaussianMeanSolution(warmGraph, warmInference, reference)
            assertGaussianMeanSolution(freshGraph, freshInference, reference)

            if warmInference isa GaussianSumProductInference
                assertGaussianMarginalsMatch(warmGraph, warmInference, reference)
                assertGaussianMarginalsMatch(freshGraph, freshInference, reference)
            end
        end
    end

    @testset "Discrete graph updates match brute-force references" begin
        for constructor in (sumproduct,)
            graph = numericDiscreteReferenceGraph()
            inference = constructor(graph)

            addVariable!(graph, inference, :x4, 2; label = "x4", states = [:cold, :warm])
            addFactor!(graph, inference, :x3, :x4, [0.8 0.2; 0.3 0.7]; label = "x3_x4")
            updateFactor!(
                graph,
                inference;
                factor = "prior_x1",
                table = [0.25, 0.35, 0.4],
                initialize = true
            )

            expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(graph)

            gbp!(graph, inference; iterations = 10, tolerance = 0.0)

            if inference isa DiscreteSumProductInference
                assertDiscreteMarginalsMatch(graph, inference, expectedMarginals)
            else
                assertDiscreteMAPMatch(graph, inference, expectedMAP)
            end
        end
    end

    @testset "Initializing unary updates reset to numeric references" begin
        for constructor in (moment, canonical, minsum)
            variables = [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2")
            ]
            factors = [
                GaussianFactor(:x1, 1.0, 1.0, 0.2; label = "prior_x1", initialize = true),
                GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.4; label = "x1_x2"),
                GaussianFactor(:x2, 2.0, 1.0, 0.3; label = "prior_x2")
            ]
            graph = factorGraph(variables, factors)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            updateFactor!(
                graph,
                inference,
                "prior_x1";
                mean = 1.4,
                covariance = 0.15,
                initialize = true
            )

            reference = normalEquationReference(graph)

            gbp!(graph, inference; iterations = 30, tolerance = 1e-12)

            assertGaussianMeanSolution(graph, inference, reference)
        end

        for constructor in (sumproduct, minsum)
            graph = numericDiscreteReferenceGraph()
            inference = constructor(graph)

            updateFactor!(
                graph,
                inference;
                factor = "prior_x1",
                table = [0.1, 0.2, 0.7],
                initialize = true
            )

            expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(graph)

            gbp!(graph, inference; iterations = 8, tolerance = 0.0)

            if inference isa DiscreteSumProductInference
                assertDiscreteMarginalsMatch(graph, inference, expectedMarginals)
            else
                assertDiscreteMAPMatch(graph, inference, expectedMAP)
            end
        end
    end

    @testset "Invalid updates leave numerical state unchanged" begin
        for constructor in (moment, canonical, minsum)
            graph = numericGaussianReferenceGraph()
            inference = constructor(graph; mean = 0.0, covariance = 100.0)
            reference = normalEquationReference(graph)

            @test_throws ErrorException addFactor!(
                graph,
                inference,
                :x2,
                :missing,
                0.5,
                [-1.0 1.0],
                0.4;
                label = "invalid"
            )
        @test_throws ErrorException updateFactor!(
            graph,
            inference,
            "difference";
            coefficient = [1.0 2.0 3.0]
        )
        @test_throws ErrorException updateFactor!(
            graph,
            inference,
            "prior_x1";
            covariance = 0.0
        )

        gbp!(graph, inference; iterations = 20, tolerance = 0.0)
        assertGaussianMeanSolution(graph, inference, reference)
    end

        graph = numericDiscreteReferenceGraph()
        inference = sumproduct(graph)
        expectedMarginals, _ = exactDiscreteMarginalsAndMAP(graph)

        @test_throws ErrorException addFactor!(
            graph,
            inference,
            :x3,
            :missing,
            [0.9 0.1; 0.2 0.8];
            label = "invalid"
        )
    @test_throws ErrorException updateFactor!(
        graph,
        inference;
        factor = "x1_x2",
        table = [0.2, 0.8]
    )
    @test_throws ErrorException updateFactor!(
        graph,
        inference;
        factor = "prior_x1",
        table = [0.0, 0.0, 0.0],
        initialize = true
    )
    @test_throws ErrorException updateFactor!(
        graph,
        inference;
        factor = "prior_x1",
        table = [0.2, -0.1, 0.9],
        initialize = true
    )

    gbp!(graph, inference; iterations = 8, tolerance = 0.0)
    assertDiscreteMarginalsMatch(graph, inference, expectedMarginals)
end

    @testset "Stale inference rejects numerical operations" begin
        gaussianGraph = numericGaussianReferenceGraph()
        gaussianInference = moment(gaussianGraph; mean = 0.0, covariance = 100.0)

        addVariable!(gaussianGraph, :x3, 1; label = "x3")

        @test_throws ErrorException messages!(gaussianGraph, gaussianInference)
        @test_throws ErrorException gbp!(gaussianGraph, gaussianInference; iterations = 1)
        @test_throws ErrorException marginals!(gaussianGraph, gaussianInference)
        @test_throws ErrorException marginalMean(gaussianGraph, gaussianInference, :x1)

        discreteGraph = numericDiscreteReferenceGraph()
        discreteInference = sumproduct(discreteGraph)

        addVariable!(discreteGraph, :x4, 2; label = "x4")

        @test_throws ErrorException messages!(discreteGraph, discreteInference)
        @test_throws ErrorException gbp!(discreteGraph, discreteInference; iterations = 1)
        @test_throws ErrorException marginals!(discreteGraph, discreteInference)
        @test_throws ErrorException marginal(discreteGraph, discreteInference, :x1)
    end

    @testset "Damping alpha extremes match direct message semantics" begin
        for constructor in (moment, canonical, minsum)
            graph = numericGaussianLoopGraph()
            undamped = constructor(graph; mean = 0.0, covariance = 100.0)
            alphaZero = constructor(graph; mean = 0.0, covariance = 100.0)
            alphaOne = constructor(graph; mean = 0.0, covariance = 100.0)

            messages!(graph, undamped)
            messages!(
                graph,
                alphaZero;
                damping = true,
                prob = 1.0,
                alpha = 0.0,
                rng = MersenneTwister(5)
            )

            assertAllFactorMessagesEqual(alphaZero, undamped)

            previous = storedMessageSnapshots(
                alphaOne.factorToVariable,
                collect(eachindex(alphaOne.factorToVariable))
            )

            messages!(
                graph,
                alphaOne;
                damping = true,
                prob = 1.0,
                alpha = 1.0,
                rng = MersenneTwister(5)
            )

            assertStoredMessagesEqual(
                alphaOne.factorToVariable,
                collect(eachindex(alphaOne.factorToVariable)),
                previous
            )
        end

        for constructor in (sumproduct, minsum)
            graph = numericDiscreteLoopGraph()
            undamped = constructor(graph)
            alphaZero = constructor(graph)
            alphaOne = constructor(graph)

            messages!(graph, undamped)
            messages!(
                graph,
                alphaZero;
                damping = true,
                prob = 1.0,
                alpha = 0.0,
                rng = MersenneTwister(5)
            )

            assertAllFactorMessagesEqual(alphaZero, undamped)

            previous = storedMessageSnapshots(
                alphaOne.factorToVariable,
                collect(eachindex(alphaOne.factorToVariable))
            )

            messages!(
                graph,
                alphaOne;
                damping = true,
                prob = 1.0,
                alpha = 1.0,
                rng = MersenneTwister(5)
            )

            assertStoredMessagesEqual(
                alphaOne.factorToVariable,
                collect(eachindex(alphaOne.factorToVariable)),
                previous
            )
        end
    end

    @testset "Manual one-step messages match closed-form references" begin
        gaussianGraph = factorGraph(
            [GaussianVariable(:x, 1; label = "x")],
            [GaussianFactor(:x, 4.0, 2.0, 0.5; label = "measurement")]
        )

        momentInference = moment(gaussianGraph; mean = 0.0, covariance = 100.0)
        canonicalInference = canonical(gaussianGraph; mean = 0.0, covariance = 100.0)
        minSumInference = minsum(gaussianGraph; mean = 0.0, covariance = 100.0)

        factorToVariableMessages!(gaussianGraph, momentInference)
        factorToVariableMessages!(gaussianGraph, canonicalInference)
        factorToVariableMessages!(gaussianGraph, minSumInference)

        @test momentInference.factorToVariable[1].mean ≈ [2.0]
        @test momentInference.factorToVariable[1].covariance ≈ reshape([0.125], 1, 1)
        @test canonicalInference.factorToVariable[1].precision ≈ reshape([8.0], 1, 1)
        @test canonicalInference.factorToVariable[1].information ≈ [16.0]
        @test minSumInference.factorToVariable[1].J ≈ reshape([8.0], 1, 1)
        @test minSumInference.factorToVariable[1].h ≈ [16.0]

        discreteGraph = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
                DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])
            ],
            [
                DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1", initialize = true),
                DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "x1_x2")
            ]
        )
        sumProductInference = sumproduct(discreteGraph)
        discreteMinSumInference = minsum(discreteGraph)
        linkX2Edge = edgeIndex(discreteGraph; variable = :x2, factor = "x1_x2")

        factorToVariableMessages!(discreteGraph, sumProductInference)
        factorToVariableMessages!(discreteGraph, discreteMinSumInference)

        @test sumProductInference.factorToVariable[linkX2Edge] ≈ [0.82, 0.34] ./ 1.16
        @test discreteMinSumInference.factorToVariable[linkX2Edge] ≈
              [0.0, -log(0.2 * 0.9) + log(0.8 * 1.0)]
    end

    @testset "Residual updateCount selects exactly one candidate" begin
        for constructor in (moment, canonical, minsum)
            graph = numericGaussianLoopGraph()
            inference = constructor(graph; mean = 0.0, covariance = 100.0)
            schedule = residualSchedule(graph, inference; updateCount = 1)

            messages!(graph, inference, schedule)

            @test length(schedule.lastUpdated) == 1
            @test only(schedule.lastUpdated) in eachindex(schedule.edgeIds)
        end

        for constructor in (sumproduct, minsum)
            graph = numericDiscreteLoopGraph()
            inference = constructor(graph)
            schedule = residualSchedule(graph, inference; updateCount = 1)

            messages!(graph, inference, schedule)

            @test length(schedule.lastUpdated) == 1
            @test only(schedule.lastUpdated) in eachindex(schedule.edgeIds)
        end
    end

    @testset "Invalid Gaussian models fail instead of producing silent numerics" begin
        underconstrained = factorGraph(
            [GaussianVariable(:x, 1; label = "x")],
            GaussianFactor[]
        )

        @test_throws ErrorException solveWLS(underconstrained)

        improper = factorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2")
            ],
            [GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.5; label = "difference")]
        )
        @test_throws ErrorException solveWLS(improper)
    end

    @testset "Gaussian results are invariant to model ordering" begin
        original = numericGaussianReferenceGraph()
        permuted = factorGraph(
            [
                GaussianVariable(:x2, 1; label = "x2"),
                GaussianVariable(:x1, 1; label = "x1")
            ],
            [
                GaussianFactor(:x1, :x2, 1.0, [-1.0 1.0], 0.5; label = "difference"),
                GaussianFactor(:x2, 3.0, 1.0, 0.25; label = "prior_x2"),
                GaussianFactor(:x1, 1.0, 1.0, 0.25; label = "prior_x1")
            ]
        )
        originalReference = normalEquationReference(original)
        permutedReference = normalEquationReference(permuted)

        for constructor in (moment, canonical, minsum)
            originalInference = constructor(original; mean = 0.0, covariance = 100.0)
            permutedInference = constructor(permuted; mean = 0.0, covariance = 100.0)

            gbp!(original, originalInference; iterations = 10, tolerance = 0.0)
            gbp!(permuted, permutedInference; iterations = 10, tolerance = 0.0)

            for variableId in (:x1, :x2)
                originalIdx = variableIndex(original, variableId)
                permutedIdx = variableIndex(permuted, variableId)

                @test originalReference.variableMean[originalIdx] ≈
                      permutedReference.variableMean[permutedIdx]

                if originalInference isa GaussianMinSumInference
                    @test estimate(original, originalInference, variableId) ≈
                          estimate(permuted, permutedInference, variableId)
                else
                    @test marginalMean(original, originalInference, variableId) ≈
                          marginalMean(permuted, permutedInference, variableId)
                    @test marginalCovariance(original, originalInference, variableId) ≈
                          marginalCovariance(permuted, permutedInference, variableId)
                end
            end
        end
    end

    @testset "Discrete zero probabilities match brute-force references" begin
        graph = factorGraph(
            [
                DiscreteVariable(:x1, 3; label = "x1", states = [:a, :b, :c]),
                DiscreteVariable(:x2, 2; label = "x2", states = [:off, :on]),
                DiscreteVariable(:x3, 2; label = "x3", states = [:low, :high])
            ],
            [
                DiscreteFactor(:x1, [0.0, 0.6, 0.4]; label = "prior_x1"),
                DiscreteFactor(:x1, :x2, [0.0 0.0; 0.7 0.3; 0.2 0.8]; label = "x1_x2"),
                DiscreteFactor(:x2, :x3, [0.9 0.1; 0.0 1.0]; label = "x2_x3")
            ]
        )
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(graph)
        sumProductInference = sumproduct(graph)
        minSumInference = minsum(graph)

        gbp!(graph, sumProductInference; iterations = 8, tolerance = 0.0)
        gbp!(graph, minSumInference; iterations = 8, tolerance = 0.0)

        assertDiscreteMarginalsMatch(graph, sumProductInference, expectedMarginals)
        assertDiscreteMAPMatch(graph, minSumInference, expectedMAP)
    end

    @testset "Discrete ternary factor matches brute-force references" begin
        graph = numericDiscreteTernaryFactorGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(graph)
        sumProductInference = sumproduct(graph)
        minSumInference = minsum(graph)

        gbp!(graph, sumProductInference; iterations = 8, tolerance = 0.0)
        gbp!(graph, minSumInference; iterations = 8, tolerance = 0.0)

        assertDiscreteMarginalsMatch(graph, sumProductInference, expectedMarginals)
        assertDiscreteMAPMatch(graph, minSumInference, expectedMAP)
    end

    @testset "Loopy Gaussian GBP converges to WLS mean" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            gbp!(graph, inference; iterations = 80, tolerance = 1e-12)

            for variableIdx in eachindex(graph.variables)
                variable = graph.variables[variableIdx]

                @test marginalMean(graph, inference, variable.id) ≈
                      reference.variableMean[variableIdx] atol = 1e-8
            end
        end

        minSumInference = minsum(graph)

        gbp!(graph, minSumInference; iterations = 80, tolerance = 1e-12)

        for variableIdx in eachindex(graph.variables)
            variable = graph.variables[variableIdx]

            @test estimate(graph, minSumInference, variable.id) ≈
                  reference.variableMean[variableIdx] atol = 1e-8
        end
    end

    @testset "Gaussian schedules recover WLS on a loopy graph" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical, minsum)
            floodingInference = constructor(graph; mean = 0.0, covariance = 100.0)

            gbp!(
                graph,
                floodingInference;
                iterations = 80,
                tolerance = 1e-12,
                schedule = :flooding
            )

            assertGaussianMeanSolution(graph, floodingInference, reference)

            for scheduleFactory in (sequentialSchedule, floodingSchedule)
                inference = constructor(graph; mean = 0.0, covariance = 100.0)
                schedule = scheduleFactory(graph, inference)

                for _ in 1:80
                    messages!(graph, inference, schedule)
                end

                refreshGaussianSolution!(graph, inference)
                assertGaussianMeanSolution(graph, inference, reference)
            end
        end
    end

    @testset "Tree dynamic updates match centralized references" begin
        for constructor in (moment, canonical, minsum)
            tree = treeFactorGraph(numericGaussianReferenceGraph(); root = :x1)
            inference = constructor(tree; mean = 0.0, covariance = 100.0)

            addVariable!(tree, inference, :x3, 1; label = "x3")
            addFactor!(tree, inference, :x2, :x3, 0.7, [-1.0 1.0], 0.35; label = "x2_x3")
            addFactor!(tree, inference, :x3, 3.2, 1.0, 0.5; label = "prior_x3")

            reference = normalEquationReference(tree.graph)

            gbp!(tree.graph, inference; iterations = 200, tolerance = 0.0)

            assertGaussianMeanSolution(tree.graph, inference, reference; atol = 1e-6)
        end

        for constructor in (sumproduct, minsum)
            tree = treeFactorGraph(numericDiscreteReferenceGraph(); root = :x1)
            inference = constructor(tree)

            addVariable!(tree, inference, :x4, 2; label = "x4", states = [:cold, :warm])
            addFactor!(tree, inference, :x3, :x4, [0.8 0.2; 0.3 0.7]; label = "x3_x4")

            expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(tree.graph)

            forwardBackward!(tree, inference)

            if inference isa DiscreteSumProductInference
                assertDiscreteMarginalsMatch(tree.graph, inference, expectedMarginals)
            else
                assertDiscreteMAPMatch(tree.graph, inference, expectedMAP)
            end
        end
    end

    @testset "Component lookup selects numerical marginal entries" begin
        graph = factorGraph(
            [
                GaussianVariable(:pose, 2; label = "pose", components = [:position, :velocity]),
                GaussianVariable(:bias, 1; label = "bias", components = [:offset])
            ],
            [
                GaussianFactor(:pose, [1.0, -0.5], Matrix{Float64}(I, 2, 2), [0.4, 0.3]; label = "prior_pose"),
                GaussianFactor(:bias, 0.2, 1.0, 0.5; label = "prior_bias"),
                GaussianFactor(
                    :pose, :bias,
                    [0.7, -0.2],
                    [1.0 0.5 -1.0; 0.2 1.0 0.4],
                    [0.25, 0.3];
                    label = "pose_bias"
                )
            ]
        )
        reference = normalEquationReference(graph)
        inference = moment(graph; mean = 0.0, covariance = 100.0)

        gbp!(graph, inference; iterations = 20, tolerance = 1e-12)

        positionIndex = componentIndex(graph, :pose, :position)
        velocityIndex = componentIndex(graph, :pose, :velocity)
        offsetIndex = componentIndex(graph, :bias, :offset)

        @test componentValue(graph, :pose, positionIndex) == :position
        @test componentValue(graph, :pose, velocityIndex) == :velocity
        @test componentValue(graph, :bias, offsetIndex) == :offset
        @test marginalMean(graph, inference, :pose)[positionIndex] ≈ reference.variableMean[1][1]
        @test marginalMean(graph, inference, :pose)[velocityIndex] ≈ reference.variableMean[1][2]
        @test marginalMean(graph, inference, :bias)[offsetIndex] ≈ reference.variableMean[2][1]
    end

    @testset "WLS comparison helpers report numerical agreement" begin
        graph = numericGaussianReferenceGraph()
        tree = treeFactorGraph(numericGaussianReferenceGraph(); root = :x1)

        for constructor in (moment, canonical)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)
            reference = solveWLS(graph)

            gbp!(graph, inference; iterations = 10, tolerance = 1e-12)

            @test maxMeanError(graph, inference, reference) <= 1e-8
            @test maxCovarianceError(graph, inference, reference) <= 1e-8
            @test quietDisplayCall() do
                compareMeanWithWLS(graph, inference, reference)
            end === nothing
            @test quietDisplayCall() do
                compareCovarianceWithWLS(graph, inference, reference)
            end === nothing

            treeInference = constructor(tree; mean = 0.0, covariance = 100.0)
            treeReference = solveWLS(tree)

            forwardBackward!(tree, treeInference)

            @test maxMeanError(tree, treeInference, treeReference) <= 1e-8
            @test maxCovarianceError(tree, treeInference, treeReference) <= 1e-8
            @test quietDisplayCall() do
                compareMeanWithWLS(tree, treeInference, treeReference)
            end === nothing
            @test quietDisplayCall() do
                compareCovarianceWithWLS(tree, treeInference, treeReference)
            end === nothing
        end
    end

    @testset "Partial residual schedules converge to numeric references" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical, minsum)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            FactorGraph.runResidualGBP!(
                graph,
                inference;
                iterations = 260,
                tolerance = 1e-12,
                updateCount = 1
            )

            assertGaussianMeanSolution(graph, inference, reference; atol = 1e-7)
        end

        discreteGraph = numericDiscreteReferenceGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(discreteGraph)

        for constructor in (sumproduct, minsum)
            inference = constructor(discreteGraph)

            FactorGraph.runResidualGBP!(
                discreteGraph,
                inference;
                iterations = 30,
                tolerance = 0.0,
                updateCount = 1
            )

            if inference isa DiscreteSumProductInference
                assertDiscreteMarginalsMatch(discreteGraph, inference, expectedMarginals)
            else
                assertDiscreteMAPMatch(discreteGraph, inference, expectedMAP)
            end
        end
    end

    @testset "Residual updateFraction schedules converge to numeric references" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical, minsum)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)
            schedule = residualSchedule(graph, inference; updateFraction = 0.25)

            for _ in 1:360
                residualStep!(graph, inference, schedule)
            end

            refreshGaussianSolution!(graph, inference)
            assertGaussianMeanSolution(graph, inference, reference; atol = 1e-7)
        end

        discreteGraph = numericDiscreteReferenceGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(discreteGraph)

        for constructor in (sumproduct, minsum)
            inference = constructor(discreteGraph)
            schedule = residualSchedule(discreteGraph, inference; updateFraction = 0.35)

            for _ in 1:35
                residualStep!(discreteGraph, inference, schedule)
            end

            if inference isa DiscreteSumProductInference
                marginals!(discreteGraph, inference)
                assertDiscreteMarginalsMatch(discreteGraph, inference, expectedMarginals)
            else
                estimates!(discreteGraph, inference)
                assertDiscreteMAPMatch(discreteGraph, inference, expectedMAP)
            end
        end
    end

    @testset "Residual Gaussian GBP reaches the same numerical solution" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical, minsum)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            FactorGraph.runResidualGBP!(
                graph,
                inference;
                iterations = 120,
                tolerance = 1e-12,
                updateFraction = 1.0
            )

            for variableIdx in eachindex(graph.variables)
                variable = graph.variables[variableIdx]

                if inference isa GaussianMinSumInference
                    @test estimate(graph, inference, variable.id) ≈
                          reference.variableMean[variableIdx] atol = 1e-8
                else
                    @test marginalMean(graph, inference, variable.id) ≈
                          reference.variableMean[variableIdx] atol = 1e-8
                end
            end
        end
    end

    @testset "Gaussian damping preserves the final numerical solution" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical, minsum)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            gbp!(
                graph,
                inference;
                iterations = 160,
                tolerance = 1e-12,
                damping = true,
                prob = 1.0,
                alpha = 0.35,
                rng = MersenneTwister(17)
            )

            for variableIdx in eachindex(graph.variables)
                variable = graph.variables[variableIdx]

                if inference isa GaussianMinSumInference
                    @test estimate(graph, inference, variable.id) ≈
                          reference.variableMean[variableIdx] atol = 1e-8
                else
                    @test marginalMean(graph, inference, variable.id) ≈
                          reference.variableMean[variableIdx] atol = 1e-8
                end
            end
        end
    end

    @testset "Gaussian edge damping preserves the final numerical solution" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical, minsum)
            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            dampEdges!(
                graph,
                inference;
                variable = :x1,
                factor = "x1_x3",
                prob = 1.0,
                alpha = 0.35
            )
            gbp!(
                graph,
                inference;
                iterations = 160,
                tolerance = 1e-12,
                damping = false,
                prob = 0.0,
                rng = MersenneTwister(29)
            )

            assertGaussianMeanSolution(graph, inference, reference)
        end
    end

    @testset "Gaussian freeze controls recover after unfreeze" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        freezeCases = (
            (
                freeze = inference -> freezeFactor!(graph, inference, "x1_x3"),
                unfreeze = inference -> unfreezeFactor!(graph, inference, "x1_x3")
            ),
            (
                freeze = inference -> freezeVariable!(graph, inference, :x2),
                unfreeze = inference -> unfreezeVariable!(graph, inference, :x2)
            ),
            (
                freeze = inference -> freezeEdge!(
                    graph,
                    inference;
                    variable = :x1,
                    factor = "x1_x3"
                ),
                unfreeze = inference -> unfreezeEdge!(
                    graph,
                    inference;
                    variable = :x1,
                    factor = "x1_x3"
                )
            )
        )

        for constructor in (moment, canonical, minsum)
            for freezeCase in freezeCases
                inference = constructor(graph; mean = 0.0, covariance = 100.0)

                freezeCase.freeze(inference)
                freezeCase.unfreeze(inference)
                gbp!(graph, inference; iterations = 120, tolerance = 1e-12)

                assertGaussianMeanSolution(graph, inference, reference)
            end

            inference = constructor(graph; mean = 0.0, covariance = 100.0)

            freezeEdge!(graph, inference; variable = :x1, factor = "x1_x3")
            gbp!(graph, inference; iterations = 5, tolerance = 0.0)
            unfreezeEdge!(graph, inference; variable = :x1, factor = "x1_x3")
            gbp!(graph, inference; iterations = 120, tolerance = 1e-12)

            assertGaussianMeanSolution(graph, inference, reference)
        end
    end

    @testset "Gaussian freeze preserves stored messages until unfreeze" begin
        graph = numericGaussianLoopGraph()
        scheduleCases = (
            (schedule = nothing, damping = false),
            (schedule = :flooding, damping = false),
            (schedule = :residual, damping = false),
            (schedule = nothing, damping = true),
            (schedule = :flooding, damping = true),
            (schedule = :residual, damping = true)
        )

        for constructor in (moment, canonical, minsum)
            for scheduleCase in scheduleCases
                inference = constructor(graph; mean = 0.0, covariance = 100.0)
                factorEdgeIds = edgeIndices(graph; factor = "x1_x3")

                gbp!(graph, inference; iterations = 2, tolerance = 0.0)

                factorSnapshots =
                    storedMessageSnapshots(inference.factorToVariable, factorEdgeIds)

                freezeFactor!(graph, inference, "x1_x3")
                exerciseFreezeSchedule!(
                    graph,
                    inference,
                    scheduleCase.schedule;
                    damping = scheduleCase.damping
                )
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    factorEdgeIds,
                    factorSnapshots
                )

                unfreezeFactor!(graph, inference, "x1_x3")
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    factorEdgeIds,
                    factorSnapshots
                )

                variableEdgeIds = edgeIndices(graph; variable = :x2)
                variableSnapshots =
                    storedMessageSnapshots(inference.variableToFactor, variableEdgeIds)

                freezeVariable!(graph, inference, :x2)
                exerciseFreezeSchedule!(
                    graph,
                    inference,
                    scheduleCase.schedule;
                    damping = scheduleCase.damping
                )
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    variableEdgeIds,
                    variableSnapshots
                )

                unfreezeVariable!(graph, inference, :x2)
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    variableEdgeIds,
                    variableSnapshots
                )

                edgeId = edgeIndex(graph; variable = :x1, factor = "x1_x3")
                variableToFactorSnapshot =
                    storedMessageSnapshots(inference.variableToFactor, [edgeId])
                factorToVariableSnapshot =
                    storedMessageSnapshots(inference.factorToVariable, [edgeId])

                freezeEdge!(graph, inference; variable = :x1, factor = "x1_x3")
                exerciseFreezeSchedule!(
                    graph,
                    inference,
                    scheduleCase.schedule;
                    damping = scheduleCase.damping
                )
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    [edgeId],
                    variableToFactorSnapshot
                )
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    [edgeId],
                    factorToVariableSnapshot
                )

                unfreezeEdge!(graph, inference; variable = :x1, factor = "x1_x3")
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    [edgeId],
                    variableToFactorSnapshot
                )
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    [edgeId],
                    factorToVariableSnapshot
                )
            end
        end
    end

    @testset "Discrete freeze preserves stored messages until unfreeze" begin
        graph = numericDiscreteLoopGraph()
        scheduleCases = (
            (schedule = nothing, damping = false),
            (schedule = :flooding, damping = false),
            (schedule = :residual, damping = false),
            (schedule = nothing, damping = true),
            (schedule = :flooding, damping = true),
            (schedule = :residual, damping = true)
        )

        for constructor in (sumproduct, minsum)
            for scheduleCase in scheduleCases
                inference = constructor(graph)
                factorEdgeIds = edgeIndices(graph; factor = "x1_x3")

                gbp!(graph, inference; iterations = 2, tolerance = 0.0)

                factorSnapshots =
                    storedMessageSnapshots(inference.factorToVariable, factorEdgeIds)

                freezeFactor!(graph, inference, "x1_x3")
                exerciseFreezeSchedule!(
                    graph,
                    inference,
                    scheduleCase.schedule;
                    damping = scheduleCase.damping
                )
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    factorEdgeIds,
                    factorSnapshots
                )

                unfreezeFactor!(graph, inference, "x1_x3")
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    factorEdgeIds,
                    factorSnapshots
                )

                variableEdgeIds = edgeIndices(graph; variable = :x2)
                variableSnapshots =
                    storedMessageSnapshots(inference.variableToFactor, variableEdgeIds)

                freezeVariable!(graph, inference, :x2)
                exerciseFreezeSchedule!(
                    graph,
                    inference,
                    scheduleCase.schedule;
                    damping = scheduleCase.damping
                )
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    variableEdgeIds,
                    variableSnapshots
                )

                unfreezeVariable!(graph, inference, :x2)
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    variableEdgeIds,
                    variableSnapshots
                )

                edgeId = edgeIndex(graph; variable = :x1, factor = "x1_x3")
                variableToFactorSnapshot =
                    storedMessageSnapshots(inference.variableToFactor, [edgeId])
                factorToVariableSnapshot =
                    storedMessageSnapshots(inference.factorToVariable, [edgeId])

                freezeEdge!(graph, inference; variable = :x1, factor = "x1_x3")
                exerciseFreezeSchedule!(
                    graph,
                    inference,
                    scheduleCase.schedule;
                    damping = scheduleCase.damping
                )
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    [edgeId],
                    variableToFactorSnapshot
                )
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    [edgeId],
                    factorToVariableSnapshot
                )

                unfreezeEdge!(graph, inference; variable = :x1, factor = "x1_x3")
                assertStoredMessagesEqual(
                    inference.variableToFactor,
                    [edgeId],
                    variableToFactorSnapshot
                )
                assertStoredMessagesEqual(
                    inference.factorToVariable,
                    [edgeId],
                    factorToVariableSnapshot
                )
            end
        end
    end

    @testset "Gaussian forward-backward matches analytic posterior" begin
        tree = treeFactorGraph(numericGaussianReferenceGraph(); root = :x1)
        expectedMean = [[1.25], [2.75]]
        expectedCovariance = [reshape([0.1875], 1, 1), reshape([0.1875], 1, 1)]

        for constructor in (moment, canonical)
            inference = constructor(tree; mean = 0.0, covariance = 100.0)

            forwardBackward!(tree, inference)

            for variableIdx in eachindex(tree.graph.variables)
                variable = tree.graph.variables[variableIdx]

                @test marginalMean(tree, inference, variable.id) ≈ expectedMean[variableIdx]
                @test marginalCovariance(tree, inference, variable.id) ≈
                      expectedCovariance[variableIdx]
            end
        end
    end

    @testset "Discrete sum-product matches brute-force marginals" begin
        graph = numericDiscreteReferenceGraph()
        tree = treeFactorGraph(graph; root = :x1)
        expectedMarginals, _ = exactDiscreteMarginalsAndMAP(graph)
        graphInference = sumproduct(graph)
        treeInference = sumproduct(tree)

        gbp!(graph, graphInference; iterations = 6, tolerance = 0.0)
        forwardBackward!(tree, treeInference)

        for variableIdx in eachindex(graph.variables)
            variable = graph.variables[variableIdx]

            @test marginal(graph, graphInference, variable.id) ≈ expectedMarginals[variableIdx]
            @test marginal(tree, treeInference, variable.id) ≈ expectedMarginals[variableIdx]
        end

        @test marginalProbability(tree, treeInference, :x1, :b) ≈ expectedMarginals[1][2]
    end

    @testset "Discrete min-sum matches brute-force MAP" begin
        graph = numericDiscreteReferenceGraph()
        tree = treeFactorGraph(graph; root = :x1)
        _, expectedMAP = exactDiscreteMarginalsAndMAP(graph)
        graphInference = minsum(graph)
        treeInference = minsum(tree)

        gbp!(graph, graphInference; iterations = 6, tolerance = 0.0)
        forwardBackward!(tree, treeInference)

        for variableIdx in eachindex(graph.variables)
            variable = graph.variables[variableIdx]

            @test estimate(graph, graphInference, variable.id) == expectedMAP[variableIdx]
            @test estimate(tree, treeInference, variable.id) == expectedMAP[variableIdx]
        end
    end

    @testset "Gaussian schedule matrix converges to the same WLS solution" begin
        graph = numericGaussianLoopGraph()
        reference = normalEquationReference(graph)

        for constructor in (moment, canonical, minsum)
            for scheduleKind in (nothing, :flooding, :residual)
                inference = constructor(graph; mean = 0.0, covariance = 100.0)

                runGaussianSchedule!(graph, inference, scheduleKind)
                assertGaussianMeanSolution(graph, inference, reference; atol = 1e-6)
            end
        end
    end

    @testset "Tree forward-backward matches centralized numeric references" begin
        for graph in (numericGaussianReferenceGraph(), numericGaussianVectorGraph())
            reference = normalEquationReference(graph)

            for root in (:x1, :x2)
                tree = treeFactorGraph(graph; root = root)

                for constructor in (moment, canonical)
                    inference = constructor(tree; mean = 0.0, covariance = 100.0)

                    forwardBackward!(tree, inference)
                    assertGaussianMarginalsMatch(tree.graph, inference, reference)
                end

                minSumInference = minsum(tree)

                forwardBackward!(tree, minSumInference)
                assertGaussianMeanSolution(tree.graph, minSumInference, reference)
            end
        end

        graph = numericDiscreteReferenceGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(graph)

        for root in (:x1, :x2, :x3)
            tree = treeFactorGraph(graph; root = root)
            sumProductInference = sumproduct(tree)
            minSumInference = minsum(tree)

            forwardBackward!(tree, sumProductInference)
            forwardBackward!(tree, minSumInference)

            assertDiscreteMarginalsMatch(tree.graph, sumProductInference, expectedMarginals)
            assertDiscreteMAPMatch(tree.graph, minSumInference, expectedMAP)
        end
    end

    @testset "Tree step-by-step sweeps match full numeric sweeps" begin
        gaussianGraph = numericGaussianVectorChainGraph()
        gaussianReference = normalEquationReference(gaussianGraph)

        for constructor in (moment, canonical, minsum)
            tree = treeFactorGraph(gaussianGraph; root = :x2)
            stepped = constructor(tree; mean = 0.0, covariance = 100.0)
            full = constructor(tree; mean = 0.0, covariance = 100.0)

            while forwardStep!(tree, stepped) !== nothing
            end

            while backwardStep!(tree, stepped) !== nothing
            end

            refreshGaussianSolution!(tree.graph, stepped)
            forwardBackward!(tree, full)

            assertGaussianMeanSolution(tree.graph, stepped, gaussianReference)
            assertGaussianMeanSolution(tree.graph, full, gaussianReference)

            if stepped isa GaussianSumProductInference
                assertGaussianMarginalsMatch(tree.graph, stepped, gaussianReference)
                assertGaussianMarginalsMatch(tree.graph, full, gaussianReference)
            end
        end

        discreteGraph = numericDiscreteReferenceGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(discreteGraph)

        for constructor in (sumproduct, minsum)
            tree = treeFactorGraph(discreteGraph; root = :x2)
            stepped = constructor(tree)
            full = constructor(tree)

            while forwardStep!(tree, stepped) !== nothing
            end

            while backwardStep!(tree, stepped) !== nothing
            end

            if stepped isa DiscreteSumProductInference
                marginals!(tree, stepped)
            else
                estimates!(tree, stepped)
            end

            forwardBackward!(tree, full)

            if stepped isa DiscreteSumProductInference
                assertDiscreteMarginalsMatch(tree.graph, stepped, expectedMarginals)
                assertDiscreteMarginalsMatch(tree.graph, full, expectedMarginals)
            else
                assertDiscreteMAPMatch(tree.graph, stepped, expectedMAP)
                assertDiscreteMAPMatch(tree.graph, full, expectedMAP)
            end
        end
    end

    @testset "Selected tree edge sweeps match full numeric sweeps" begin
        gaussianGraph = numericGaussianVectorChainGraph()
        gaussianReference = normalEquationReference(gaussianGraph)

        for constructor in (moment, canonical, minsum)
            tree = treeFactorGraph(gaussianGraph; root = :x2)
            selected = constructor(tree; mean = 0.0, covariance = 100.0)
            full = constructor(tree; mean = 0.0, covariance = 100.0)

            for edgeId in tree.forwardOrder
                @test forwardStep!(tree, selected, edgeId) == edgeId
            end

            for edgeId in tree.backwardOrder
                @test backwardStep!(tree, selected, edgeId) == edgeId
            end

            refreshGaussianSolution!(tree.graph, selected)
            forwardBackward!(tree, full)

            assertGaussianMeanSolution(tree.graph, selected, gaussianReference)
            assertGaussianMeanSolution(tree.graph, full, gaussianReference)

            if selected isa GaussianSumProductInference
                assertGaussianMarginalsMatch(tree.graph, selected, gaussianReference)
                assertGaussianMarginalsMatch(tree.graph, full, gaussianReference)
            end
        end

        discreteGraph = numericDiscreteReferenceGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(discreteGraph)

        for constructor in (sumproduct, minsum)
            tree = treeFactorGraph(discreteGraph; root = :x2)
            selected = constructor(tree)
            full = constructor(tree)

            for edgeId in tree.forwardOrder
                @test forwardStep!(tree, selected, edgeId) == edgeId
            end

            for edgeId in tree.backwardOrder
                @test backwardStep!(tree, selected, edgeId) == edgeId
            end

            if selected isa DiscreteSumProductInference
                marginals!(tree, selected)
            else
                estimates!(tree, selected)
            end

            forwardBackward!(tree, full)

            if selected isa DiscreteSumProductInference
                assertDiscreteMarginalsMatch(tree.graph, selected, expectedMarginals)
                assertDiscreteMarginalsMatch(tree.graph, full, expectedMarginals)
            else
                assertDiscreteMAPMatch(tree.graph, selected, expectedMAP)
                assertDiscreteMAPMatch(tree.graph, full, expectedMAP)
            end
        end
    end

    @testset "Tree root refresh preserves numeric results" begin
        gaussianGraph = numericGaussianVectorChainGraph()
        gaussianReference = normalEquationReference(gaussianGraph)
        gaussianTree = treeFactorGraph(gaussianGraph; root = :x1)

        for root in (:x1, :x2, :x3)
            refresh!(gaussianTree; root = root)

            for constructor in (moment, canonical, minsum)
                inference = constructor(gaussianTree; mean = 0.0, covariance = 100.0)

                forwardBackward!(gaussianTree, inference)
                assertGaussianMeanSolution(gaussianTree.graph, inference, gaussianReference)

                if inference isa GaussianSumProductInference
                    assertGaussianMarginalsMatch(gaussianTree.graph, inference, gaussianReference)
                end
            end
        end

        discreteGraph = numericDiscreteReferenceGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(discreteGraph)
        discreteTree = treeFactorGraph(discreteGraph; root = :x1)

        for root in (:x1, :x2, :x3)
            refresh!(discreteTree; root = root)

            for constructor in (sumproduct, minsum)
                inference = constructor(discreteTree)

                forwardBackward!(discreteTree, inference)

                if inference isa DiscreteSumProductInference
                    assertDiscreteMarginalsMatch(discreteTree.graph, inference, expectedMarginals)
                else
                    assertDiscreteMAPMatch(discreteTree.graph, inference, expectedMAP)
                end
            end
        end
    end

    @testset "Tree and graph-only topology changes reject stale inference" begin
        gaussianTree = treeFactorGraph(numericGaussianReferenceGraph(); root = :x1)
        gaussianInference = canonical(gaussianTree; mean = 0.0, covariance = 100.0)
        gaussianSchedule = sequentialSchedule(gaussianTree, gaussianInference)

        addVariable!(gaussianTree.graph, :x3, 1; label = "x3")
        addFactor!(gaussianTree.graph, :x2, :x3, 0.7, [-1.0 1.0], 0.3; label = "x2_x3")
        refresh!(gaussianTree)

        @test_throws ErrorException messages!(gaussianTree, gaussianInference, gaussianSchedule)
        @test_throws ErrorException forwardBackward!(gaussianTree, gaussianInference)
        @test_throws ErrorException marginalMean(gaussianTree, gaussianInference, :x1)

        discreteTree = treeFactorGraph(numericDiscreteReferenceGraph(); root = :x1)
        discreteInference = sumproduct(discreteTree)
        discreteSchedule = floodingSchedule(discreteTree, discreteInference)

        addVariable!(discreteTree.graph, :x4, 2; label = "x4", states = [:cold, :warm])
        addFactor!(
            discreteTree.graph,
            :x3,
            :x4,
            [0.8 0.2; 0.3 0.7];
            label = "x3_x4"
        )
        refresh!(discreteTree)

        @test_throws ErrorException messages!(discreteTree, discreteInference, discreteSchedule)
        @test_throws ErrorException forwardBackward!(discreteTree, discreteInference)
        @test_throws ErrorException marginal(discreteTree, discreteInference, :x1)
    end

    @testset "Edge damping probability extremes match direct message semantics" begin
        for constructor in (moment, canonical, minsum)
            graph = numericGaussianLoopGraph()
            undamped = constructor(graph; mean = 0.0, covariance = 100.0)
            neverDamped = constructor(graph; mean = 0.0, covariance = 100.0)
            alwaysDamped = constructor(graph; mean = 0.0, covariance = 100.0)
            edgeId = edgeIndex(graph; variable = :x1, factor = "x1_x3")

            messages!(graph, undamped)
            dampEdges!(
                graph,
                neverDamped;
                variable = :x1,
                factor = "x1_x3",
                prob = 0.0,
                alpha = 1.0
            )
            messages!(graph, neverDamped; rng = MersenneTwister(9))
            assertAllFactorMessagesEqual(neverDamped, undamped)

            previous = storedMessageSnapshots(alwaysDamped.factorToVariable, [edgeId])

            dampEdges!(
                graph,
                alwaysDamped;
                variable = :x1,
                factor = "x1_x3",
                prob = 1.0,
                alpha = 1.0
            )
            messages!(graph, alwaysDamped; rng = MersenneTwister(9))
            assertStoredMessagesEqual(alwaysDamped.factorToVariable, [edgeId], previous)
        end
    end

    @testset "Discrete min-sum tie estimates have optimal brute-force weight" begin
        variables = [
            DiscreteVariable(:x1, 2; label = "x1", states = [:a, :b]),
            DiscreteVariable(:x2, 2; label = "x2", states = [:a, :b])
        ]
        factors = [
            DiscreteFactor(:x1, [0.5, 0.5]; label = "prior_x1"),
            DiscreteFactor(:x2, [0.5, 0.5]; label = "prior_x2"),
            DiscreteFactor(:x1, :x2, [2.0 0.5; 0.5 2.0]; label = "tie")
        ]
        graph = factorGraph(variables, factors)
        tree = treeFactorGraph(graph; root = :x1)
        bestWeight = bestDiscreteWeight(graph)

        for graphLike in (graph, tree)
            inference = minsum(graphLike)

            forwardBackward!(tree, inference)

            states = [estimate(graphLike, inference, variable.id) for variable in graph.variables]

            @test states in ([:a, :a], [:b, :b])
            @test discreteAssignmentWeight(graph, states) ≈ bestWeight
        end
    end

    @testset "Lookup helpers remain aligned after inference-aware updates" begin
        graph = numericGaussianReferenceGraph()
        inference = moment(graph; mean = 0.0, covariance = 100.0)

        addVariable!(
            graph,
            inference,
            GaussianVariable(:x3, 2; label = "pose", components = [:p, :v])
        )
        addFactor!(
            graph,
            inference,
            :x3,
            [2.0, -0.5],
            Matrix{Float64}(I, 2, 2),
            [0.2, 0.3];
            label = "prior_pose"
        )

        reference = normalEquationReference(graph)

        gbp!(graph, inference; iterations = 30, tolerance = 0.0)

        poseIndex = variableIndex(graph, :x3)
        positionIndex = componentIndex(graph, :x3, :p)
        velocityIndex = componentIndex(graph, "pose", :v)

        @test marginalMean(graph, inference, "pose") ≈ reference.variableMean[poseIndex]
        @test marginalMean(graph, inference, :x3)[positionIndex] ≈
              reference.variableMean[poseIndex][positionIndex]
        @test marginalMean(graph, inference, "pose")[velocityIndex] ≈
              reference.variableMean[poseIndex][velocityIndex]
    end

    @testset "Loopy discrete inference remains close to brute-force reference" begin
        graph = numericDiscreteLoopGraph()
        expectedMarginals, expectedMAP = exactDiscreteMarginalsAndMAP(graph)
        sumProductInference = sumproduct(graph)
        minSumInference = minsum(graph)

        gbp!(graph, sumProductInference; iterations = 60, tolerance = 1e-12)
        gbp!(graph, minSumInference; iterations = 60, tolerance = 1e-12)

        for variableIdx in eachindex(graph.variables)
            variable = graph.variables[variableIdx]

            @test marginal(graph, sumProductInference, variable.id) ≈
                  expectedMarginals[variableIdx] atol = 0.04
            @test estimate(graph, minSumInference, variable.id) == expectedMAP[variableIdx]
        end
    end
end
