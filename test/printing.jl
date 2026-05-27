include("setup.jl")

function capturedPrint(f)
    mktemp() do path, output
        redirect_stdout(output) do
            f()
        end

        flush(output)

        return read(path, String)
    end
end

function printingDiscreteGraph()
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

@testset verbose = true "Printing helpers" begin
    @testset "Prints empty graphs" begin
        gaussianOutput = capturedPrint() do
            printGraph(GaussianFactorGraph())
        end
        discreteOutput = capturedPrint() do
            printGraph(DiscreteFactorGraph())
        end

        @test occursin("Number of variable nodes: 0", gaussianOutput)
        @test occursin("Number of variable nodes: 0", discreteOutput)
    end

    @testset "Prints graph models and edges" begin
        gaussianGraph = gaussianTreeTestGraph()
        discreteGraph = printingDiscreteGraph()

        gaussianModelOutput = capturedPrint() do
            printModel(gaussianGraph)
        end
        discreteModelOutput = capturedPrint() do
            printModel(discreteGraph)
        end
        edgeOutput = capturedPrint() do
            printEdges(gaussianGraph)
        end

        @test occursin("GaussianVariable nodes", gaussianModelOutput)
        @test occursin("GaussianFactor nodes", gaussianModelOutput)
        @test occursin("DiscreteVariable nodes", discreteModelOutput)
        @test occursin("DiscreteFactor nodes", discreteModelOutput)
        @test occursin("x1 <-> prior_x1", edgeOutput)
    end

    @testset "Prints Gaussian marginals, estimates, and WLS results" begin
        graph = gaussianTreeTestGraph()
        momentInference = moment(graph)
        canonicalInference = canonical(graph)
        minSumInference = minsum(graph)

        gbp!(graph, momentInference; iterations = 3)
        gbp!(graph, canonicalInference; iterations = 3)
        gbp!(graph, minSumInference; iterations = 3)

        momentOutput = capturedPrint() do
            printMarginal(graph, momentInference; variable = :x1)
        end
        canonicalOutput = capturedPrint() do
            printMarginal(graph, canonicalInference; variable = :x1)
        end
        estimateOutput = capturedPrint() do
            printEstimate(graph, minSumInference; variable = :x1)
        end
        wlsOutput = capturedPrint() do
            printWLS(graph, solveWLS(graph))
        end

        @test occursin("moment form", momentOutput)
        @test occursin("mean =", momentOutput)
        @test occursin("canonical form", canonicalOutput)
        @test occursin("precision =", canonicalOutput)
        @test occursin("MAP estimate", estimateOutput)
        @test occursin("WLS result", wlsOutput)
        @test occursin("Estimate for GaussianVariable x1", wlsOutput)
    end

    @testset "Prints discrete sum-product messages and marginals" begin
        graph = printingDiscreteGraph()
        inference = sumproduct(graph)

        messages!(graph, inference)
        marginals!(graph, inference)

        messageOutput = capturedPrint() do
            printMessages(graph, inference; variable = :x1)
        end
        marginalOutput = capturedPrint() do
            printMarginal(graph, inference; variable = :x1)
        end

        @test occursin("sum-product", messageOutput)
        @test occursin("probability =", messageOutput)
        @test occursin("sum-product form", marginalOutput)
        @test occursin("probability =", marginalOutput)
    end

    @testset "Prints Gaussian min-sum messages" begin
        graph = gaussianTreeTestGraph()
        inference = minsum(graph)

        messages!(graph, inference)

        variableOutput = capturedPrint() do
            printMessages(graph, inference; variable = :x1)
        end
        factorOutput = capturedPrint() do
            printMessages(graph, inference; factor = "prior_x1")
        end

        @test occursin("min-sum", variableOutput)
        @test occursin("J =", variableOutput)
        @test occursin("h =", factorOutput)
    end

    @testset "Prints discrete min-sum messages and estimates" begin
        graph = printingDiscreteGraph()
        inference = minsum(graph)

        messages!(graph, inference)
        estimates!(graph, inference)

        variableOutput = capturedPrint() do
            printMessages(graph, inference; variable = :x1)
        end
        factorOutput = capturedPrint() do
            printMessages(graph, inference; factor = "prior_x1")
        end
        estimateOutput = capturedPrint() do
            printEstimate(graph, inference; variable = :x1)
        end

        @test occursin("cost =", variableOutput)
        @test occursin("cost =", factorOutput)
        @test occursin("estimate =", estimateOutput)
    end
end
