include("setup.jl")

@testset verbose = true "Manual message updates" begin
    @testset "Moment messages then marginals matches gbp!" begin
        graph = gaussianTestGraph()
        automatic = moment(graph; mean = 0.0, covariance = 1e6)
        manual = moment(graph; mean = 0.0, covariance = 1e6)

        gbp!(graph, automatic; iterations = 20)

        for _ in 1:20
            messages!(graph, manual)
        end
        marginals!(graph, manual)

        for variableIndex in eachindex(graph.variables)
            @test manual.marginal[variableIndex].mean ==
                automatic.marginal[variableIndex].mean
            @test manual.marginal[variableIndex].covariance ==
                automatic.marginal[variableIndex].covariance
        end
    end

    @testset "Moment directional messages then marginals matches gbp!" begin
        graph = gaussianTestGraph()
        automatic = moment(graph; mean = 0.0, covariance = 1e6)
        manual = moment(graph; mean = 0.0, covariance = 1e6)

        gbp!(graph, automatic; iterations = 20)

        for _ in 1:20
            factorToVariableMessages!(graph, manual)
            variableToFactorMessages!(graph, manual)
        end
        marginals!(graph, manual)

        for variableIndex in eachindex(graph.variables)
            @test manual.marginal[variableIndex].mean ==
                automatic.marginal[variableIndex].mean
            @test manual.marginal[variableIndex].covariance ==
                automatic.marginal[variableIndex].covariance
        end
    end

    @testset "Canonical messages then marginals matches gbp!" begin
        graph = gaussianTestGraph()
        automatic = canonical(graph; mean = 0.0, covariance = 1e6)
        manual = canonical(graph; mean = 0.0, covariance = 1e6)

        gbp!(graph, automatic; iterations = 20)

        for _ in 1:20
            messages!(graph, manual)
        end
        marginals!(graph, manual)

        for variableIndex in eachindex(graph.variables)
            @test manual.marginal[variableIndex].mean ==
                automatic.marginal[variableIndex].mean
            @test manual.marginal[variableIndex].covariance ==
                automatic.marginal[variableIndex].covariance
        end
    end

    @testset "Canonical directional messages then marginals matches gbp!" begin
        graph = gaussianTestGraph()
        automatic = canonical(graph; mean = 0.0, covariance = 1e6)
        manual = canonical(graph; mean = 0.0, covariance = 1e6)

        gbp!(graph, automatic; iterations = 20)

        for _ in 1:20
            factorToVariableMessages!(graph, manual)
            variableToFactorMessages!(graph, manual)
        end
        marginals!(graph, manual)

        for variableIndex in eachindex(graph.variables)
            @test manual.marginal[variableIndex].mean ==
                automatic.marginal[variableIndex].mean
            @test manual.marginal[variableIndex].covariance ==
                automatic.marginal[variableIndex].covariance
        end
    end

    @testset "Canonical directional broadcast messages then marginals matches gbp!" begin
        graph = gaussianTestGraph()
        automatic = canonical(graph; mean = 0.0, covariance = 1e6)
        manual = canonical(graph; mean = 0.0, covariance = 1e6)

        gbp!(graph, automatic; iterations = 20, broadcast = true)

        for _ in 1:20
            factorToVariableMessages!(graph, manual; broadcast = true)
            variableToFactorMessages!(graph, manual; broadcast = true)
        end
        marginals!(graph, manual)

        for variableIndex in eachindex(graph.variables)
            @test manual.marginal[variableIndex].mean ==
                automatic.marginal[variableIndex].mean
            @test manual.marginal[variableIndex].covariance ==
                automatic.marginal[variableIndex].covariance
        end
    end
end
