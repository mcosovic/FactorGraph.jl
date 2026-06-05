include("setup.jl")

@testset verbose = true "Scalar inputs" begin
    @testset "One-dimensional graph accepts scalar GaussianFactor arguments" begin
        graph = GaussianFactorGraph(
            [
                GaussianVariable(:x1, 1; label = "x1", mean = 0.0, covariance = 1.0)
            ],
            [
                GaussianFactor(:x1, 0.25, 1.0, 0.15; label = "prior_x1")
            ]
        )

        factorData = graph.factors[1]

        @test graph.variables[1].mean == [0.0]
        @test graph.variables[1].covariance == reshape([1.0], 1, 1)
        @test graph.variables[1].components == [1]
        @test factorData.mean == [0.25]
        @test factorData.coefficient == reshape([1.0], 1, 1)
        @test factorData.covariance == reshape([0.15], 1, 1)
        @test solveWLS(graph).variableMean[1] == [0.25]

        vectorFactor = GaussianFactor(:x1, [0.25], [1.0], [0.15])
        @test vectorFactor.mean == [0.25]
        @test vectorFactor.coefficient == reshape([1.0], 1, 1)
        @test vectorFactor.covariance == reshape([0.15], 1, 1)

        isotropicFactor = GaussianFactor(:x1, [0.25, 0.5], Matrix{Float64}(I, 2, 2), 0.15)
        @test isotropicFactor.covariance == 0.15 * Matrix{Float64}(I, 2, 2)

        @test_throws ErrorException GaussianFactor(:x1, [0.25, 0.5], [1.0, 2.0], [0.15])
        @test_throws ErrorException GaussianFactor(:x1, 0.0, 1.0, "bad")
    end

    @testset "Gaussian variables accept named components" begin
        graph = GaussianFactorGraph(
            [
                GaussianVariable(
                    :x1,
                    2;
                    label = "x1",
                    components = [:position, :velocity]
                )
            ],
            [
                GaussianFactor(:x1, [0.25, 1.0], Matrix{Float64}(I, 2, 2), [0.15, 0.2])
            ]
        )

        variable = graph.variables[1]

        @test variable.components == [:position, :velocity]
        @test componentIndex(graph, :x1, :position) == 1
        @test componentIndex(graph, :x1, :velocity) == 2
        @test componentValue(graph, :x1, 2) == :velocity
        @test_throws ErrorException GaussianVariable(:bad, 2; components = [:x])
        @test_throws ErrorException GaussianVariable(:bad, 2; components = [:x, :x])
        @test_throws ErrorException GaussianVariable(:bad, 1; components = [""])
        @test_throws ErrorException GaussianVariable(:bad, 0)
        @test_throws ErrorException GaussianVariable(:bad, 1; label = "")
        @test_throws ErrorException GaussianVariable(:bad, 2; components = :x)
        @test_throws ErrorException GaussianVariable(:bad, 2; components = Symbol[])
        @test_throws ErrorException GaussianVariable(:bad, 2; components = [1.0, 2.0])
        @test_throws ErrorException GaussianVariable(:bad, 2; mean = [0.0])
        @test_throws ErrorException GaussianVariable(:bad, 2; covariance = [1.0 0.0; 0.0 -1.0])
    end

    @testset "Default Gaussian variable labels preserve IDs" begin
        graph = GaussianFactorGraph(
            [GaussianVariable(:x_1, 1)],
            [GaussianFactor(:x_1, 0.25, 1.0, 0.15)]
        )

        @test graph.variables[1].label == "x_1"
        @test graph.factors[1].label == "f1"

        explicit = GaussianVariable(:x_2, 1; label = "x_2")
        @test explicit.label == "x_2"
    end

    @testset "Gaussian variable nodes can be used to construct factors" begin
        x1 = GaussianVariable(:x1, 1; label = "x1")
        x2 = GaussianVariable(:x2, 1; label = "x2")

        prior = GaussianFactor(x1, 0.25, 1.0, 0.15; label = "prior_x1")
        link = GaussianFactor(x1, x2, 0.0, [1.0 -1.0], 0.5; label = "link")
        graph = GaussianFactorGraph([x1, x2], [prior, link])

        @test prior.variables == [:x1]
        @test link.variables == [:x1, :x2]
        @test variableIndex(graph, :x1) == 1
        @test variableDimension(graph, :x2) == 1
        @test coefficientBlock(graph; factor = "link", variable = :x2) == reshape([-1.0], 1, 1)
        @test edgeIndex(graph; variable = :x2, factor = "link") == 3
        @test edgeIndices(graph; variable = :x1) == [1, 2]

        momentInference = moment(graph)
        canonicalInference = canonical(graph)
        minSumInference = minsum(graph)

        @test length(marginal(graph, momentInference, :x1).mean) == 1
        @test length(marginal(graph, canonicalInference, :x2).mean) == 1
        @test length(estimate(graph, minSumInference, :x1)) == 1
        @test_throws ErrorException GaussianFactor(DiscreteVariable(:d1, 2), 0.0, 1.0, 0.1)
    end

    @testset "Inference defaults accept scalar mean and covariance" begin
        graph = gaussianTreeTestGraph()

        momentInference = moment(graph; mean = 0.2, covariance = 10.0)
        canonicalInference = canonical(graph; mean = 0.2, covariance = 10.0)

        for variableIndex in eachindex(graph.variables)
            dimension = graph.variables[variableIndex].dimension

            @test momentInference.initial[variableIndex].mean == fill(0.2, dimension)
            @test momentInference.initial[variableIndex].covariance ==
                10.0 * Matrix{Float64}(I, dimension, dimension)
            @test canonicalInference.initial[variableIndex].precision ==
                0.1 * Matrix{Float64}(I, dimension, dimension)
        end
    end
end
