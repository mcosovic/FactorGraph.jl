include("setup.jl")

function discreteInferenceGraph()
    variables = [
        DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
        DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])
    ]
    factors = [
        DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1", initialize = true),
        DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "x1_x2")
    ]

    return factorGraph(variables, factors)
end

function ternaryDiscreteInferenceGraph()
    table = cat(
        [0.9 0.2; 0.3 0.4],
        [0.1 0.8; 0.7 0.6];
        dims = 3
    )

    variables = [
        DiscreteVariable(:x1, 2; label = "x1", states = [:a, :b]),
        DiscreteVariable(:x2, 2; label = "x2", states = [:c, :d]),
        DiscreteVariable(:x3, 2; label = "x3", states = [:e, :f])
    ]
    factors = [
        DiscreteFactor(:x1, [0.6, 0.4]; label = "prior_x1", initialize = true),
        DiscreteFactor(:x1, :x2, :x3, table; label = "ternary")
    ]

    return factorGraph(variables, factors)
end

function bruteForceDiscreteMarginals(graph::DiscreteFactorGraph)
    marginals = [zeros(variable.cardinality) for variable in graph.variables]
    total = 0.0
    ranges = [1:variable.cardinality for variable in graph.variables]

    for assignment in Iterators.product(ranges...)
        weight = 1.0

        for factorData in graph.factors
            factorIndexTuple = map(factorData.variables) do variableRef
                variableIdx = variableIndex(graph, variableRef)

                return assignment[variableIdx]
            end

            weight *= factorData.table[factorIndexTuple...]
        end

        total += weight

        for variableIdx in eachindex(marginals)
            marginals[variableIdx][assignment[variableIdx]] += weight
        end
    end

    for marginal in marginals
        marginal ./= total
    end

    return marginals
end

function bruteForceDiscreteMAP(graph::DiscreteFactorGraph)
    bestWeight = -Inf
    bestAssignment = nothing
    ranges = [1:variable.cardinality for variable in graph.variables]

    for assignment in Iterators.product(ranges...)
        weight = 1.0

        for factorData in graph.factors
            factorIndexTuple = map(factorData.variables) do variableRef
                variableIdx = variableIndex(graph, variableRef)

                return assignment[variableIdx]
            end

            weight *= factorData.table[factorIndexTuple...]
        end

        if weight > bestWeight
            bestWeight = weight
            bestAssignment = assignment
        end
    end

    return StateRef[
        graph.variables[index].states[bestAssignment[index]]
        for index in eachindex(graph.variables)
    ]
end

@testset verbose = true "Discrete sum-product inference" begin
    @testset "Constructs inference state" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)

        @test inference isa DiscreteSumProductInference
        @test inference isa DiscreteInference
        @test inference isa AbstractSumProductInference
        @test inference.initial[1] == [0.8, 0.2]
        @test inference.initial[2] == [0.5, 0.5]
        @test length(inference.variableToFactor) == length(graph.edges)
        @test length(inference.factorToVariable) == length(graph.edges)
        @test marginal(graph, inference, :x1) == [0.8, 0.2]
        @test marginalProbability(graph, inference, :x1, :on) == 0.2
    end

    @testset "Updates messages manually" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)

        factorToVariableMessages!(graph, inference)
        variableToFactorMessages!(graph, inference)
        marginals!(graph, inference)

        priorEdge = edgeIndex(graph; variable = :x1, factor = "prior_x1")
        linkX2Edge = edgeIndex(graph; variable = :x2, factor = "x1_x2")

        @test inference.factorToVariable[priorEdge] == [0.8, 0.2]
        @test inference.factorToVariable[linkX2Edge] ≈ [0.82, 0.34] ./ 1.16
        @test inference.variableToFactor[linkX2Edge] == [0.5, 0.5]
    end

    @testset "Runs gbp! to brute-force marginals on a tree" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)
        expected = bruteForceDiscreteMarginals(graph)

        gbp!(graph, inference; iterations = 5, tolerance = 0.0)

        @test marginal(graph, inference, :x1) ≈ expected[1]
        @test marginal(graph, inference, :x2) ≈ expected[2]
        @test marginals(inference) === inference.marginal
    end

    @testset "Handles ternary factor tables" begin
        graph = ternaryDiscreteInferenceGraph()
        inference = sumproduct(graph)
        expected = bruteForceDiscreteMarginals(graph)

        gbp!(graph, inference; iterations = 5, tolerance = 0.0)

        @test marginal(graph, inference, :x1) ≈ expected[1]
        @test marginal(graph, inference, :x2) ≈ expected[2]
        @test marginal(graph, inference, :x3) ≈ expected[3]
    end

    @testset "Applies global and selected damping" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)
        linkX2Edge = edgeIndex(graph; variable = :x2, factor = "x1_x2")

        factorToVariableMessages!(
            graph,
            inference;
            damping = true,
            prob = 1.0,
            alpha = 0.5
        )

        @test inference.factorToVariable[linkX2Edge] ≈
              0.5 .* [0.5, 0.5] .+ 0.5 .* ([0.82, 0.34] ./ 1.16)

        dampEdges!(graph, inference; variable = :x2, factor = "x1_x2", prob = 1.0, alpha = 0.5)
        @test areDampedEdges(graph, inference; variable = :x2, factor = "x1_x2")

        undampEdges!(graph, inference; variable = :x2, factor = "x1_x2")
        @test !areDampedEdges(graph, inference; variable = :x2, factor = "x1_x2")

        dampEdges!(graph, inference; prob = 1.0, alpha = 0.25)
        @test all(inference.dampedEdges)
        @test all(==(0.25), inference.edgeDampingAlpha)

        undampEdges!(graph, inference)
        @test !any(inference.dampedEdges)
    end

    @testset "Updates an initializing unary factor with inference" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)

        updateFactor!(
            graph,
            inference;
            factor = "prior_x1",
            table = [0.2, 0.8],
            initialize = true
        )

        @test inference.initial[1] == [0.2, 0.8]

        gbp!(graph, inference; iterations = 5)
        expected = bruteForceDiscreteMarginals(graph)

        @test marginal(graph, inference, :x1) ≈ expected[1]
        @test marginal(graph, inference, :x2) ≈ expected[2]
    end

    @testset "Freezes selected updates" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)
        priorEdge = edgeIndex(graph; variable = :x1, factor = "prior_x1")
        linkX1Edge = edgeIndex(graph; variable = :x1, factor = "x1_x2")

        freezeFactor!(graph, inference, "prior_x1")
        @test isFrozenFactor(graph, inference, "prior_x1")

        updateFactor!(graph; factor = "prior_x1", table = [0.2, 0.8], initialize = true)
        factorToVariableMessages!(graph, inference)

        @test inference.factorToVariable[priorEdge] == [0.8, 0.2]

        unfreezeFactor!(graph, inference, "prior_x1")
        @test !isFrozenFactor(graph, inference, "prior_x1")

        freezeVariable!(graph, inference, :x1)
        @test isFrozenVariable(graph, inference, :x1)

        variableMessage = copy(inference.variableToFactor[linkX1Edge])
        variableToFactorMessages!(graph, inference)

        @test inference.variableToFactor[linkX1Edge] == variableMessage

        unfreezeVariable!(graph, inference, :x1)
        @test !isFrozenVariable(graph, inference, :x1)

        freezeEdge!(graph, inference; variable = :x1, factor = "x1_x2")
        @test isFrozenEdge(graph, inference; variable = :x1, factor = "x1_x2")

        unfreezeEdge!(graph, inference; variable = :x1, factor = "x1_x2")
        @test !isFrozenEdge(graph, inference; variable = :x1, factor = "x1_x2")
    end

    @testset "Extends inference after graph additions" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)
        oldVariableToFactor = deepcopy(inference.variableToFactor)
        oldFactorToVariable = deepcopy(inference.factorToVariable)

        addVariable!(graph, inference, :x3, 2; label = "x3", states = [:a, :b])
        addFactor!(graph, inference, :x2, :x3, [0.9 0.1; 0.1 0.9]; label = "x2_x3")

        @test length(inference.initial) == length(graph.variables)
        @test length(inference.variableToFactor) == length(graph.edges)
        @test length(inference.factorToVariable) == length(graph.edges)
        @test isempty(inference.nextVariableToFactor)
        @test isempty(inference.nextFactorToVariable)
        @test inference.variableToFactor[1:length(oldVariableToFactor)] == oldVariableToFactor
        @test inference.factorToVariable[1:length(oldFactorToVariable)] == oldFactorToVariable

        gbp!(graph, inference; iterations = 6)
        expected = bruteForceDiscreteMarginals(graph)

        @test marginal(graph, inference, :x1) ≈ expected[1]
        @test marginal(graph, inference, :x2) ≈ expected[2]
        @test marginal(graph, inference, :x3) ≈ expected[3]
    end

    @testset "Tree forward-backward matches brute-force marginals" begin
        tree = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
                DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high]),
                DiscreteVariable(:x3, 2; label = "x3")
            ],
            [
                DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1"),
                DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "x1_x2"),
                DiscreteFactor(:x2, :x3, [0.9 0.3; 0.2 0.8]; label = "x2_x3")
            ];
            root = :x1
        )
        inference = sumproduct(tree)
        expected = bruteForceDiscreteMarginals(tree.graph)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)
        previousMarginals = deepcopy(inference.marginal)

        schedule = forwardBackward!(tree, inference)

        @test schedule.forwardIndex == length(tree.forwardOrder) + 1
        @test schedule.backwardIndex == length(tree.backwardOrder) + 1
        @test marginal(tree, inference, :x1) ≈ expected[1]
        @test marginal(tree, inference, :x2) ≈ expected[2]
        @test marginal(tree, inference, :x3) ≈ expected[3]
        @test marginalProbability(tree, inference, :x1, :on) ≈ expected[1][2]
        @test maxMessageChange(
            tree,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
        @test maxMarginalChange(tree, inference, previousMarginals) > 0.0
    end

    @testset "Tree forward and backward steps" begin
        graph = discreteInferenceGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = sumproduct(graph)
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

    @testset "Tree cursor and selected steps" begin
        graph = discreteInferenceGraph()
        tree = treeFactorGraph(graph; root = :x1)
        stepped = sumproduct(graph)
        selected = sumproduct(graph)

        firstForwardEdge = tree.forwardOrder[1]
        @test forwardStep!(tree, stepped) == firstForwardEdge
        @test forwardStep!(tree, selected, firstForwardEdge) == firstForwardEdge
        @test selected.variableToFactor[firstForwardEdge] ==
              stepped.variableToFactor[firstForwardEdge]
        @test selected.factorToVariable[firstForwardEdge] ==
              stepped.factorToVariable[firstForwardEdge]

        reset!(tree)

        firstBackwardEdge = tree.backwardOrder[1]
        edge = graph.edges[firstBackwardEdge]
        variable = graph.variables[edge.variableIndex].id
        factor = graph.factors[edge.factorIndex].label

        @test backwardStep!(tree, stepped) == firstBackwardEdge
        @test backwardStep!(tree, selected; variable = variable, factor = factor) ==
              firstBackwardEdge
        @test selected.variableToFactor[firstBackwardEdge] ==
              stepped.variableToFactor[firstBackwardEdge]
        @test selected.factorToVariable[firstBackwardEdge] ==
              stepped.factorToVariable[firstBackwardEdge]
    end

    @testset "Tree freeze helpers delegate to discrete inference state" begin
        graph = discreteInferenceGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = sumproduct(tree)

        freezeFactor!(tree, inference, "x1_x2")
        @test isFrozenFactor(tree, inference, "x1_x2")
        unfreezeFactor!(tree, inference, "x1_x2")
        @test !isFrozenFactor(tree, inference, "x1_x2")

        freezeVariable!(tree, inference, :x1)
        @test isFrozenVariable(tree, inference, :x1)
        unfreezeVariable!(tree, inference, :x1)
        @test !isFrozenVariable(tree, inference, :x1)

        freezeEdge!(tree, inference; variable = :x1, factor = "x1_x2")
        @test isFrozenEdge(tree, inference; variable = :x1, factor = "x1_x2")
        unfreezeEdge!(tree, inference; variable = :x1, factor = "x1_x2")
        @test !isFrozenEdge(tree, inference; variable = :x1, factor = "x1_x2")

        mapInference = minsum(tree)
        freezeEdge!(tree, mapInference; variable = :x1, factor = "x1_x2")
        @test isFrozenEdge(tree, mapInference; variable = :x1, factor = "x1_x2")
    end

    @testset "Tree damping helpers delegate to discrete inference state" begin
        graph = discreteInferenceGraph()
        tree = treeFactorGraph(graph; root = :x1)
        inference = sumproduct(tree)

        edgeId = edgeIndex(graph; variable = :x1, factor = "x1_x2")

        @test !areDampedEdges(tree, inference; variable = :x1, factor = "x1_x2")

        dampEdges!(
            tree,
            inference;
            variable = :x1,
            factor = "x1_x2",
            prob = 1.0,
            alpha = 0.25
        )
        @test areDampedEdges(tree, inference; variable = :x1, factor = "x1_x2")
        @test inference.edgeDampingProb[edgeId] == 1.0
        @test inference.edgeDampingAlpha[edgeId] == 0.25

        undampEdges!(tree, inference; variable = :x1, factor = "x1_x2")
        @test !areDampedEdges(tree, inference; variable = :x1, factor = "x1_x2")

        dampEdges!(tree, inference; prob = 1.0, alpha = 0.3)
        @test all(inference.dampedEdges)
        @test all(==(0.3), inference.edgeDampingAlpha)

        undampEdges!(tree, inference)
        @test !any(inference.dampedEdges)

        mapInference = minsum(tree)
        dampEdges!(
            tree,
            mapInference;
            variable = :x1,
            factor = "x1_x2",
            prob = 1.0
        )
        @test areDampedEdges(tree, mapInference; variable = :x1, factor = "x1_x2")
    end

    @testset "Tree warm-start addition refreshes view" begin
        tree = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1"),
                DiscreteVariable(:x2, 2; label = "x2")
            ],
            [
                DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1"),
                DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "x1_x2")
            ];
            root = :x1
        )
        inference = sumproduct(tree)

        forwardBackward!(tree, inference)
        oldEdgeCount = length(tree.graph.edges)

        addVariable!(tree, inference, :x3, 2; label = "x3", states = [:a, :b])
        @test length(tree.forwardOrder) == oldEdgeCount

        addFactor!(tree, inference, :x2, :x3, [0.9 0.1; 0.1 0.9]; label = "x2_x3")

        @test length(tree.graph.edges) == oldEdgeCount + 2
        @test length(tree.forwardOrder) == length(tree.graph.edges)
        @test length(tree.backwardOrder) == length(tree.graph.edges)

        forwardBackward!(tree, inference)
        expected = bruteForceDiscreteMarginals(tree.graph)

        @test marginal(tree, inference, :x1) ≈ expected[1]
        @test marginal(tree, inference, :x2) ≈ expected[2]
        @test marginal(tree, inference, :x3) ≈ expected[3]
    end

    @testset "Detects stale topology" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)

        addVariable!(graph, :x3, 2; label = "x3")

        @test_throws ErrorException gbp!(graph, inference; iterations = 1)
    end

    @testset "Reports message and marginal changes" begin
        graph = discreteInferenceGraph()
        inference = sumproduct(graph)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)
        previousMarginals = deepcopy(inference.marginal)

        gbp!(graph, inference; iterations = 1)

        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
        @test maxMarginalChange(graph, inference, previousMarginals) > 0.0
    end
end

@testset verbose = true "Discrete min-sum MAP inference" begin
    @testset "Constructs inference state" begin
        graph = discreteInferenceGraph()
        inference = minsum(graph)

        @test inference isa DiscreteMinSumInference
        @test inference isa DiscreteInference
        @test inference isa AbstractMinSumInference
        @test length(inference.variableToFactor) == length(graph.edges)
        @test length(inference.factorToVariable) == length(graph.edges)
        @test estimate(graph, inference, :x1) == :off
        @test estimates(graph, inference) === inference.estimate
    end

    @testset "Runs schedules to brute-force MAP" begin
        graph = discreteInferenceGraph()
        expected = bruteForceDiscreteMAP(graph)

        for schedule in (:sequential, :flooding, :residual)
            inference = minsum(graph)
            gbp!(
                graph,
                inference;
                iterations = 10,
                tolerance = 0.0,
                schedule = schedule,
                updateFraction = schedule == :residual ? 1.0 : nothing
            )

            @test estimates(graph, inference) == expected
        end
    end

    @testset "Handles ternary factor tables" begin
        graph = ternaryDiscreteInferenceGraph()
        inference = minsum(graph)
        expected = bruteForceDiscreteMAP(graph)

        gbp!(graph, inference; iterations = 5, tolerance = 0.0)

        @test estimates(graph, inference) == expected
    end

    @testset "Supports damping and freezing" begin
        graph = discreteInferenceGraph()
        inference = minsum(graph)
        priorEdge = edgeIndex(graph; variable = :x1, factor = "prior_x1")

        dampEdges!(graph, inference; variable = :x1, factor = "prior_x1", prob = 1.0)
        @test areDampedEdges(graph, inference; variable = :x1, factor = "prior_x1")

        freezeFactor!(graph, inference, "prior_x1")
        @test isFrozenFactor(graph, inference, "prior_x1")

        oldMessage = copy(inference.factorToVariable[priorEdge])
        updateFactor!(graph; factor = "prior_x1", table = [0.2, 0.8], initialize = true)
        factorToVariableMessages!(graph, inference)

        @test inference.factorToVariable[priorEdge] == oldMessage

        unfreezeFactor!(graph, inference, "prior_x1")
        @test !isFrozenFactor(graph, inference, "prior_x1")

        undampEdges!(graph, inference; variable = :x1, factor = "prior_x1")
        @test !areDampedEdges(graph, inference; variable = :x1, factor = "prior_x1")
    end

    @testset "Tree forward-backward matches brute-force MAP" begin
        tree = factorGraph(
            [
                DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
                DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high]),
                DiscreteVariable(:x3, 2; label = "x3")
            ],
            [
                DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1"),
                DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "x1_x2"),
                DiscreteFactor(:x2, :x3, [0.9 0.3; 0.2 0.8]; label = "x2_x3")
            ];
            root = :x1
        )
        inference = minsum(tree)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)
        previousEstimates = copy(inference.estimate)
        expected = bruteForceDiscreteMAP(tree.graph)

        forwardBackward!(tree, inference)

        @test estimates(tree, inference) == expected
        @test estimate(tree, inference, :x1) == expected[1]
        @test maxMessageChange(
            tree,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
        @test maxEstimateChange(tree, inference, previousEstimates) >= 0.0
    end

    @testset "Extends inference after graph additions" begin
        graph = discreteInferenceGraph()
        inference = minsum(graph)

        addVariable!(graph, inference, :x3, 2; label = "x3", states = [:a, :b])
        addFactor!(graph, inference, :x2, :x3, [0.9 0.1; 0.1 0.9]; label = "x2_x3")

        @test length(inference.initial) == length(graph.variables)
        @test length(inference.variableToFactor) == length(graph.edges)
        @test isempty(inference.nextVariableToFactor)

        gbp!(graph, inference; iterations = 10)

        @test estimates(graph, inference) == bruteForceDiscreteMAP(graph)
    end
end
