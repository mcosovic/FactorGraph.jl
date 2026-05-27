include("setup.jl")

@testset verbose = true "Diagnostics" begin
    @testset "Moment message and marginal changes" begin
        graph = gaussianTestGraph()
        inference = moment(graph; mean = 0.0, covariance = 1e6)

        unchangedVariableMessages = deepcopy(inference.variableToFactor)
        unchangedFactorMessages = deepcopy(inference.factorToVariable)
        unchangedMarginals = deepcopy(inference.marginal)

        @test maxMessageChange(
            graph,
            inference,
            unchangedVariableMessages,
            unchangedFactorMessages
        ) == 0.0
        @test maxVariableMessageChange(graph, inference, unchangedVariableMessages) == 0.0
        @test maxFactorMessageChange(graph, inference, unchangedFactorMessages) == 0.0
        @test maxMarginalChange(graph, inference, unchangedMarginals) == 0.0

        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)
        previousMarginals = deepcopy(inference.marginal)
        messages!(graph, inference)

        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
        @test maxMarginalChange(graph, inference, previousMarginals) == 0.0

        previousInference = deepcopy(inference)
        previousMarginals = deepcopy(inference.marginal)
        marginals!(graph, inference)

        @test maxMarginalChange(graph, inference, previousMarginals) > 0.0
        @test maxMarginalChange(graph, inference, previousInference) ==
              maxMarginalChange(graph, inference, previousMarginals)
    end

    @testset "Canonical message and marginal changes" begin
        graph = gaussianTestGraph()
        inference = canonical(graph; mean = 0.0, covariance = 1e6)

        unchangedVariableMessages = deepcopy(inference.variableToFactor)
        unchangedFactorMessages = deepcopy(inference.factorToVariable)
        unchangedMarginals = deepcopy(inference.marginal)

        @test maxMessageChange(
            graph,
            inference,
            unchangedVariableMessages,
            unchangedFactorMessages
        ) == 0.0
        @test maxVariableMessageChange(graph, inference, unchangedVariableMessages) == 0.0
        @test maxFactorMessageChange(graph, inference, unchangedFactorMessages) == 0.0
        @test maxMarginalChange(graph, inference, unchangedMarginals) == 0.0

        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)
        previousMarginals = deepcopy(inference.marginal)
        messages!(graph, inference)

        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
        @test maxMarginalChange(graph, inference, previousMarginals) == 0.0

        previousInference = deepcopy(inference)
        previousMarginals = deepcopy(inference.marginal)
        marginals!(graph, inference)

        @test maxMarginalChange(graph, inference, previousMarginals) > 0.0
        @test maxMarginalChange(graph, inference, previousInference) ==
              maxMarginalChange(graph, inference, previousMarginals)
    end

    @testset "Gaussian min-sum message and estimate changes" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        unchangedVariableMessages = deepcopy(inference.variableToFactor)
        unchangedFactorMessages = deepcopy(inference.factorToVariable)
        unchangedEstimates = [copy(estimate) for estimate in inference.estimate]

        @test maxMessageChange(
            graph,
            inference,
            unchangedVariableMessages,
            unchangedFactorMessages
        ) == 0.0
        @test maxVariableMessageChange(graph, inference, unchangedVariableMessages) == 0.0
        @test maxFactorMessageChange(graph, inference, unchangedFactorMessages) == 0.0
        @test maxEstimateChange(graph, inference, unchangedEstimates) == 0.0

        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)
        previousEstimates = [copy(estimate) for estimate in inference.estimate]

        messages!(graph, inference)

        @test maxMessageChange(
            graph,
            inference,
            previousVariableMessages,
            previousFactorMessages
        ) > 0.0
        @test maxEstimateChange(graph, inference, previousEstimates) == 0.0

        previousEstimates = [copy(estimate) for estimate in inference.estimate]
        estimates!(graph, inference)

        @test maxEstimateChange(graph, inference, previousEstimates) > 0.0
    end

    @testset "GBP tolerance stops on marginal change for moment inference" begin
        graph = gaussianTreeTestGraph()
        automatic = moment(graph; mean = 0.0, covariance = 1e6)
        manual = moment(graph; mean = 0.0, covariance = 1e6)
        tolerance = 1e-4

        manualIterations = 0

        for _ in 1:100
            previousMarginals = deepcopy(manual.marginal)
            messages!(graph, manual)
            marginals!(graph, manual)
            manualIterations += 1

            if maxMarginalChange(graph, manual, previousMarginals) <= tolerance
                break
            end
        end

        gbp!(graph, automatic; iterations = 100, tolerance = tolerance)

        @test manualIterations < 100

        for GaussianVariable in (:x1, :x2, :x3)
            @test isapprox(
                marginalMean(graph, automatic, GaussianVariable),
                marginalMean(graph, manual, GaussianVariable)
            )
            @test isapprox(
                marginalCovariance(graph, automatic, GaussianVariable),
                marginalCovariance(graph, manual, GaussianVariable)
            )
        end
    end

    @testset "GBP tolerance stops on marginal change for canonical inference" begin
        graph = gaussianTreeTestGraph()
        automatic = canonical(graph; mean = 0.0, covariance = 1e6)
        manual = canonical(graph; mean = 0.0, covariance = 1e6)
        tolerance = 1e-4

        manualIterations = 0

        for _ in 1:100
            previousMarginals = deepcopy(manual.marginal)
            messages!(graph, manual)
            marginals!(graph, manual)
            manualIterations += 1

            if maxMarginalChange(graph, manual, previousMarginals) <= tolerance
                break
            end
        end

        gbp!(graph, automatic; iterations = 100, tolerance = tolerance)

        @test manualIterations < 100

        for GaussianVariable in (:x1, :x2, :x3)
            @test isapprox(
                marginalMean(graph, automatic, GaussianVariable),
                marginalMean(graph, manual, GaussianVariable)
            )
            @test isapprox(
                marginalCovariance(graph, automatic, GaussianVariable),
                marginalCovariance(graph, manual, GaussianVariable)
            )
        end
    end

    @testset "GBP tolerance validation" begin
        graph = gaussianTreeTestGraph()

        @test_throws ErrorException gbp!(graph, moment(graph); tolerance = -1e-6)
        @test_throws ErrorException gbp!(graph, canonical(graph); tolerance = "small")
    end

    @testset "Tree overloads delegate to the underlying graph" begin
        graph = factorGraph(
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
        inference = canonical(graph)
        previousVariableMessages = deepcopy(inference.variableToFactor)
        previousFactorMessages = deepcopy(inference.factorToVariable)
        previousMarginals = deepcopy(inference.marginal)

        forwardBackward!(graph, inference)

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
        @test maxMarginalChange(graph, inference, previousMarginals) ==
              maxMarginalChange(graph.graph, inference, previousMarginals)
    end

    @testset "Rejects mismatched inference forms" begin
        graph = gaussianTestGraph()
        momentInferenceState = moment(graph)
        canonicalInferenceState = canonical(graph)

        @test_throws MethodError maxMessageChange(
            graph,
            momentInferenceState,
            canonicalInferenceState
        )
        @test_throws MethodError maxMarginalChange(
            graph,
            momentInferenceState,
            canonicalInferenceState
        )
    end
end
