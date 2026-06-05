include("setup.jl")

function residualDiscreteGraph()
    return factorGraph(
        [
            DiscreteVariable(:x1, 2; label = "x1", states = [:off, :on]),
            DiscreteVariable(:x2, 2; label = "x2", states = [:low, :high])
        ],
        [
            DiscreteFactor(:x1, [0.8, 0.2]; label = "prior_x1", initialize = true),
            DiscreteFactor(:x1, :x2, [1.0 0.2; 0.1 0.9]; label = "x1_x2")
        ]
    )
end

function residualDiscreteMarginals(graph::DiscreteFactorGraph)
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

@testset verbose = true "Residual Gaussian scheduling" begin
    @testset "Stateless schedules validate matching inference" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph; mean = 0.0, covariance = 1e6)
        tree = treeFactorGraph(graph; root = :x1)

        @test sequentialSchedule(graph, inference) isa SequentialSchedule
        @test floodingSchedule(graph, inference) isa FloodingSchedule
        @test sequentialSchedule(tree, inference) isa SequentialSchedule
        @test floodingSchedule(tree, inference) isa FloodingSchedule

        addVariable!(graph, GaussianVariable(:x4, 1; label = "x4"))

        @test_throws ErrorException sequentialSchedule(graph, inference)
        @test_throws ErrorException floodingSchedule(graph, inference)
    end

    @testset "Builds and validates residual schedule" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph; mean = 0.0, covariance = 1e6)

        schedule = residualSchedule(graph, inference; updateFraction = 0.25)

        @test length(schedule.edgeIds) == 2 * length(graph.edges)
        @test length(schedule.directions) == length(schedule.edgeIds)
        @test length(schedule.residuals) == length(schedule.edgeIds)
        @test schedule.updateFraction == 0.25
        @test schedule.updateCount === nothing

        @test_throws ErrorException residualSchedule(graph, inference; updateFraction = 0.0)
        @test_throws ErrorException residualSchedule(graph, inference; updateFraction = 1.2)
        @test_throws ErrorException residualSchedule(graph, inference; updateCount = 0)
        @test_throws ErrorException residualSchedule(
            graph,
            inference;
            updateFraction = 0.2,
            updateCount = 2
        )
    end

    @testset "Updates top fraction of moment residuals" begin
        graph = gaussianTreeTestGraph()
        inference = moment(graph; mean = 0.0, covariance = 1e6)
        schedule = residualSchedule(graph, inference; updateFraction = 0.5)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)

        messages!(graph, inference, schedule)

        @test length(schedule.lastUpdated) ==
              ceil(Int, 0.5 * 2 * length(graph.edges))
        @test maximum(schedule.residuals) > 0.0
        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
    end

    @testset "Updates requested count of canonical residuals" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(graph; mean = 0.0, covariance = 1e6)
        schedule = residualSchedule(graph, inference; updateCount = 3)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)

        messages!(graph, inference, schedule)

        @test length(schedule.lastUpdated) == 3
        @test maximum(schedule.residuals) > 0.0
        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
    end

    @testset "Keyword residual messages match manual Gaussian schedules" begin
        graph = gaussianTreeTestGraph()

        for keywordInference in (
            moment(graph; mean = 0.0, covariance = 1e6),
            canonical(graph; mean = 0.0, covariance = 1e6),
            minsum(graph)
        )
            manualInference = deepcopy(keywordInference)
            schedule = residualSchedule(graph, manualInference; updateCount = 3)

            messages!(graph, keywordInference; schedule = :residual, updateCount = 3)
            messages!(graph, manualInference, schedule)

            @test maxMessageChange(graph, keywordInference, manualInference) == 0.0
        end
    end

    @testset "Manual Gaussian sequential schedule matches keyword" begin
        graph = gaussianTreeTestGraph()

        for keywordInference in (
            moment(graph; mean = 0.0, covariance = 1e6),
            canonical(graph; mean = 0.0, covariance = 1e6),
            minsum(graph)
        )
            manualInference = deepcopy(keywordInference)
            schedule = sequentialSchedule(graph, manualInference)

            messages!(graph, keywordInference; schedule = :sequential)
            messages!(graph, manualInference, schedule)

            @test maxMessageChange(graph, keywordInference, manualInference) == 0.0
        end
    end

    @testset "Manual Gaussian tree sequential schedule delegates to graph" begin
        tree = treeFactorGraph(gaussianTreeTestGraph(); root = :x1)
        automatic = moment(tree; mean = 0.0, covariance = 1e6)
        manual = moment(tree; mean = 0.0, covariance = 1e6)
        schedule = sequentialSchedule(tree, manual)

        messages!(tree.graph, automatic; schedule = :sequential)
        messages!(tree, manual, schedule)

        @test maxMessageChange(tree, automatic, manual) == 0.0
    end

    @testset "Residual Gaussian tree schedules and steps delegate to graph" begin
        tree = treeFactorGraph(gaussianTreeTestGraph(); root = :x1)

        for constructor in (
            graph -> moment(graph; mean = 0.0, covariance = 1e6),
            graph -> canonical(graph; mean = 0.0, covariance = 1e6),
            minsum
        )
            treeInference = constructor(tree)
            graphInference = constructor(tree.graph)
            treeSchedule = residualSchedule(tree, treeInference; updateCount = 3)
            graphSchedule = residualSchedule(tree.graph, graphInference; updateCount = 3)

            residualStep!(tree, treeInference, treeSchedule)
            residualStep!(tree.graph, graphInference, graphSchedule)

            @test maxMessageChange(tree, treeInference, graphInference) == 0.0
            @test treeSchedule.lastUpdated == graphSchedule.lastUpdated
        end
    end

    @testset "Residual Gaussian tree messages delegate to graph" begin
        tree = treeFactorGraph(gaussianTreeTestGraph(); root = :x1)
        treeInference = moment(tree; mean = 0.0, covariance = 1e6)
        graphInference = moment(tree.graph; mean = 0.0, covariance = 1e6)
        treeSchedule = residualSchedule(tree, treeInference; updateCount = 3)
        graphSchedule = residualSchedule(tree.graph, graphInference; updateCount = 3)

        messages!(tree, treeInference, treeSchedule)
        messages!(tree.graph, graphInference, graphSchedule)

        @test maxMessageChange(tree, treeInference, graphInference) == 0.0
        @test treeSchedule.lastUpdated == graphSchedule.lastUpdated
    end

    @testset "Residual Gaussian tree GBP delegates to graph" begin
        tree = treeFactorGraph(gaussianTreeTestGraph(); root = :x1)

        for constructor in (
            graph -> moment(graph; mean = 0.0, covariance = 1e6),
            graph -> canonical(graph; mean = 0.0, covariance = 1e6),
            minsum
        )
            treeInference = constructor(tree)
            graphInference = constructor(tree.graph)

            treeSchedule = FactorGraph.runResidualGBP!(
                tree,
                treeInference;
                iterations = 2,
                updateCount = 3
            )
            graphSchedule = FactorGraph.runResidualGBP!(
                tree.graph,
                graphInference;
                iterations = 2,
                updateCount = 3
            )

            @test maxMessageChange(tree, treeInference, graphInference) == 0.0
            @test treeSchedule.lastUpdated == graphSchedule.lastUpdated
        end
    end

    @testset "Residual Gaussian updates honor per-edge damping" begin
        graph = gaussianTreeTestGraph()

        for inference in (
            moment(graph; mean = 0.0, covariance = 1e6),
            canonical(graph; mean = 0.0, covariance = 1e6),
            minsum(graph)
        )
            dampEdges!(graph, inference; prob = 1.0, alpha = 0.25)
            schedule = residualSchedule(graph, inference; updateFraction = 1.0)

            residualStep!(graph, inference, schedule; prob = 0.0, alpha = 0.8)

            @test all(inference.dampedEdges)
            @test length(schedule.lastUpdated) == 2 * length(graph.edges)
        end
    end

    @testset "Residual GBP produces Gaussian tree MAP means" begin
        graph = gaussianTreeTestGraph()
        inference = canonical(graph; mean = 0.0, covariance = 1e6)
        reference = solveWLS(graph)

        gbp!(
            graph,
            inference;
            schedule = :residual,
            iterations = 30,
            updateFraction = 1.0,
            tolerance = 1e-8
        )

        for variableIndex in eachindex(graph.variables)
            variable = graph.variables[variableIndex]

            @test isapprox(
                marginalMean(graph, inference, variable.id),
                reference.variableMean[variableIndex];
                atol = 1e-5
            )
        end
    end

    @testset "gbp! selects Gaussian schedules" begin
        graph = gaussianTreeTestGraph()
        residualInference = canonical(graph; mean = 0.0, covariance = 1e6)
        floodingInference = moment(graph; mean = 0.0, covariance = 1e6)

        gbp!(
            graph,
            residualInference;
            iterations = 5,
            schedule = :residual,
            updateCount = 3
        )

        @test maxMessageChange(
            graph,
            residualInference,
            canonical(graph; mean = 0.0, covariance = 1e6)
        ) > 0.0

        gbp!(graph, floodingInference; iterations = 1, schedule = :flooding)

        @test length(floodingInference.nextVariableToFactor) == length(graph.edges)
        @test length(floodingInference.nextFactorToVariable) == length(graph.edges)
        @test_throws ErrorException gbp!(graph, floodingInference; schedule = :unknown)
    end

    @testset "Residual scheduling rejects broadcast mode" begin
        graph = gaussianTreeTestGraph()

        for inference in (
            moment(graph; mean = 0.0, covariance = 1e6),
            canonical(graph; mean = 0.0, covariance = 1e6),
            minsum(graph)
        )
            schedule = residualSchedule(graph, inference; updateCount = 2)

            @test_throws ErrorException messages!(
                graph,
                inference,
                schedule;
                broadcast = true
            )
            @test_throws ErrorException gbp!(
                graph,
                inference;
                schedule = :residual,
                broadcast = true
            )
        end
    end
end

@testset verbose = true "Residual discrete scheduling" begin
    @testset "Builds and validates residual schedule" begin
        graph = residualDiscreteGraph()
        inference = sumproduct(graph)

        schedule = residualSchedule(graph, inference; updateFraction = 0.25)

        @test length(schedule.edgeIds) == 2 * length(graph.edges)
        @test length(schedule.directions) == length(schedule.edgeIds)
        @test length(schedule.residuals) == length(schedule.edgeIds)
        @test schedule.updateFraction == 0.25
        @test schedule.updateCount === nothing

        @test_throws ErrorException residualSchedule(graph, inference; updateFraction = 0.0)
        @test_throws ErrorException residualSchedule(graph, inference; updateFraction = 1.2)
        @test_throws ErrorException residualSchedule(graph, inference; updateCount = 0)
        @test_throws ErrorException residualSchedule(
            graph,
            inference;
            updateFraction = 0.2,
            updateCount = 2
        )
    end

    @testset "Updates requested count of discrete residuals" begin
        graph = residualDiscreteGraph()
        inference = sumproduct(graph)
        schedule = residualSchedule(graph, inference; updateCount = 2)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)

        messages!(graph, inference, schedule)

        @test length(schedule.lastUpdated) == 2
        @test maximum(schedule.residuals) > 0.0
        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
    end

    @testset "Rejects stale discrete residual schedules" begin
        graph = residualDiscreteGraph()
        inference = sumproduct(graph)
        schedule = residualSchedule(graph, inference; updateCount = 2)

        addVariable!(graph, inference, :x3, 2)
        addFactor!(graph, inference, :x2, :x3, [0.9 0.1; 0.1 0.9])

        @test_throws ErrorException messages!(graph, inference, schedule)
    end

    @testset "Keyword residual messages match manual discrete schedules" begin
        graph = residualDiscreteGraph()

        for keywordInference in (sumproduct(graph), minsum(graph))
            manualInference = deepcopy(keywordInference)
            schedule = residualSchedule(graph, manualInference; updateCount = 2)

            messages!(graph, keywordInference; schedule = :residual, updateCount = 2)
            messages!(graph, manualInference, schedule)

            @test maxMessageChange(graph, keywordInference, manualInference) == 0.0
        end
    end

    @testset "Manual discrete sequential schedule matches keyword" begin
        graph = residualDiscreteGraph()

        for keywordInference in (sumproduct(graph), minsum(graph))
            manualInference = deepcopy(keywordInference)
            schedule = sequentialSchedule(graph, manualInference)

            messages!(graph, keywordInference; schedule = :sequential)
            messages!(graph, manualInference, schedule)

            @test maxMessageChange(graph, keywordInference, manualInference) == 0.0
        end
    end

    @testset "Manual discrete tree sequential schedule delegates to graph" begin
        tree = treeFactorGraph(residualDiscreteGraph(); root = :x1)

        for constructor in (sumproduct, minsum)
            automatic = constructor(tree)
            manual = constructor(tree)
            schedule = sequentialSchedule(tree, manual)

            messages!(tree.graph, automatic; schedule = :sequential)
            messages!(tree, manual, schedule)

            @test maxMessageChange(tree, automatic, manual) == 0.0
        end
    end

    @testset "Residual discrete tree steps delegate to graph" begin
        tree = treeFactorGraph(residualDiscreteGraph(); root = :x1)

        for constructor in (sumproduct, minsum)
            treeInference = constructor(tree)
            graphInference = constructor(tree.graph)
            treeSchedule = residualSchedule(tree, treeInference; updateCount = 2)
            graphSchedule = residualSchedule(tree.graph, graphInference; updateCount = 2)

            residualStep!(tree, treeInference, treeSchedule)
            residualStep!(tree.graph, graphInference, graphSchedule)

            @test maxMessageChange(tree, treeInference, graphInference) == 0.0
            @test treeSchedule.lastUpdated == graphSchedule.lastUpdated
        end
    end

    @testset "Residual discrete tree GBP delegates to graph" begin
        tree = treeFactorGraph(residualDiscreteGraph(); root = :x1)

        for constructor in (sumproduct, minsum)
            treeInference = constructor(tree)
            graphInference = constructor(tree.graph)

            treeSchedule = FactorGraph.runResidualGBP!(
                tree,
                treeInference;
                iterations = 2,
                updateCount = 2
            )
            graphSchedule = FactorGraph.runResidualGBP!(
                tree.graph,
                graphInference;
                iterations = 2,
                updateCount = 2
            )

            @test maxMessageChange(tree, treeInference, graphInference) == 0.0
            @test treeSchedule.lastUpdated == graphSchedule.lastUpdated
        end
    end

    @testset "Residual discrete updates honor per-edge damping" begin
        graph = residualDiscreteGraph()

        for inference in (sumproduct(graph), minsum(graph))
            dampEdges!(graph, inference; prob = 1.0, alpha = 0.25)
            schedule = residualSchedule(graph, inference; updateFraction = 1.0)

            residualStep!(graph, inference, schedule; prob = 0.0, alpha = 0.8)

            @test all(inference.dampedEdges)
            @test length(schedule.lastUpdated) == 2 * length(graph.edges)
        end
    end

    @testset "Residual GBP matches brute-force marginals on a tree" begin
        graph = residualDiscreteGraph()
        inference = sumproduct(graph)
        expected = residualDiscreteMarginals(graph)

        gbp!(
            graph,
            inference;
            schedule = :residual,
            iterations = 10,
            updateFraction = 1.0,
            tolerance = 1e-10
        )

        @test marginal(graph, inference, :x1) ≈ expected[1]
        @test marginal(graph, inference, :x2) ≈ expected[2]
    end

    @testset "Residual tree overload delegates to underlying discrete graph" begin
        graph = treeFactorGraph(residualDiscreteGraph(); root = :x1)
        inference = sumproduct(graph)
        schedule = residualSchedule(graph, inference; updateCount = 2)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)

        messages!(graph, inference, schedule)
        marginals!(graph, inference)

        @test length(schedule.lastUpdated) == 2
        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) ==
              maxMessageChange(
                  graph.graph,
                  inference,
                  previousVariableMessages,
                  previousFactorMessages
              )
    end

    @testset "Residual tree overload supports discrete min-sum" begin
        graph = treeFactorGraph(residualDiscreteGraph(); root = :x1)
        inference = minsum(graph)
        schedule = residualSchedule(graph, inference; updateCount = 2)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)

        messages!(graph, inference, schedule)
        estimates!(graph, inference)

        @test length(schedule.lastUpdated) == 2
        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) ==
              maxMessageChange(
                  graph.graph,
                  inference,
                  previousVariableMessages,
                  previousFactorMessages
              )
    end

    @testset "gbp! selects discrete schedules" begin
        graph = residualDiscreteGraph()
        residualInference = sumproduct(graph)
        floodingInference = sumproduct(graph)
        initialInference = sumproduct(graph)

        gbp!(
            graph,
            residualInference;
            iterations = 5,
            schedule = :residual,
            updateCount = 2
        )

        @test maxMessageChange(graph, residualInference, sumproduct(graph)) > 0.0

        gbp!(graph, floodingInference; iterations = 1, schedule = :flooding)

        @test length(floodingInference.nextVariableToFactor) == length(graph.edges)
        @test length(floodingInference.nextFactorToVariable) == length(graph.edges)
        @test_throws ErrorException gbp!(graph, floodingInference; schedule = :unknown)

        @test gbp!(
            graph,
            initialInference;
            iterations = 0,
            schedule = :residual,
            updateCount = 2
        ) === nothing

        @test maxMessageChange(graph, initialInference, sumproduct(graph)) == 0.0
        @test_throws ErrorException gbp!(
            graph,
            sumproduct(graph);
            iterations = -1,
            schedule = :residual
        )
    end

    @testset "Manual discrete flooding schedule matches gbp!" begin
        graph = residualDiscreteGraph()
        automatic = sumproduct(graph)
        manual = sumproduct(graph)
        schedule = floodingSchedule(graph, manual)

        gbp!(graph, automatic; iterations = 1, schedule = :flooding)
        messages!(graph, manual, schedule)
        marginals!(graph, manual)

        @test maxMessageChange(graph, manual, automatic) == 0.0
        @test maxMarginalChange(graph, manual, automatic) == 0.0
    end

    @testset "Manual discrete min-sum flooding schedule matches gbp!" begin
        graph = residualDiscreteGraph()
        automatic = minsum(graph)
        manual = minsum(graph)
        schedule = floodingSchedule(graph, manual)

        gbp!(graph, automatic; iterations = 1, schedule = :flooding)
        messages!(graph, manual, schedule)
        estimates!(graph, manual)

        @test maxMessageChange(graph, manual, automatic) == 0.0
        @test maxEstimateChange(graph, manual, automatic) == 0.0
    end
end
