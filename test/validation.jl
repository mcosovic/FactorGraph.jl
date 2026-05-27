include("setup.jl")

@testset verbose = true "Validation helpers" begin
    @testset "Residuals use current marginal means" begin
        graph = GaussianFactorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1"),
                GaussianVariable(:x2, 1; label = "x2")
            ],
            [
                GaussianFactor(:x1, 1.0, 1.0, 0.25; label = "prior_x1"),
                GaussianFactor(:x2, 2.0, 1.0, 0.36; label = "prior_x2"),
                GaussianFactor(:x1, :x2, 0.8, [-1.0 1.0], 0.16; label = "difference")
            ]
        )
        inference = canonical(graph; mean = 0.0, covariance = 1e6)

        gbp!(graph, inference; iterations = 6)

        result = residuals(graph, inference)
        normalized = normalizedResiduals(graph, inference)
        x1 = marginalMean(graph, inference, :x1)
        x2 = marginalMean(graph, inference, :x2)

        @test length(result) == 3
        @test result[1].factor == "prior_x1"
        @test result[1].index == 1
        @test result[1].value == [1.0] - x1
        @test result[2].value == [2.0] - x2
        @test result[3].value == [0.8] - [-1.0 1.0] * vcat(x1, x2)

        @test normalized[1].value == result[1].value ./ 0.5
        @test normalized[2].value == result[2].value ./ 0.6
        @test normalized[3].value == result[3].value ./ 0.4
    end

    @testset "Residuals support tree graph views" begin
        graph = gaussianTreeTestGraph()
        tree = factorGraph(graph.variables, graph.factors; root = :x1)
        inference = canonical(tree; mean = 0.0, covariance = 1e6)

        forwardBackward!(tree, inference)

        @test residuals(tree, inference) == residuals(tree.graph, inference)
        @test normalizedResiduals(tree, inference) ==
              normalizedResiduals(tree.graph, inference)
    end
end
