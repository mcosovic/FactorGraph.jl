using LinearAlgebra
using Test

if !isdefined(@__MODULE__, :FactorGraph) || !isdefined(FactorGraph, :treeFactorGraph)
    include(joinpath(@__DIR__, "..", "src", "FactorGraph.jl"))
end

using .FactorGraph

if !isdefined(@__MODULE__, :SCHEDULE_CASES)
    const SCHEDULE_CASES = [
        (name = "standard", broadcast = false, flooding = false),
        (name = "broadcast", broadcast = true, flooding = false),
        (name = "standard + flooding", broadcast = false, flooding = true),
        (name = "broadcast + flooding", broadcast = true, flooding = true)
    ]
end

if !isdefined(@__MODULE__, :gbpSchedule)
    gbpSchedule(flooding::Bool) = flooding ? :flooding : nothing
end

if !isdefined(@__MODULE__, :gaussianTestGraph)
    function gaussianTestGraph()
        I2 = Matrix{Float64}(I, 2, 2)

        variables = [
            GaussianVariable(:x1, 1; label = "x1"),
            GaussianVariable(:x2, 2; label = "x2"),
            GaussianVariable(:x3, 4; label = "x3"),
            GaussianVariable(:x4, 3; label = "x4")
        ]

        factors = [
            GaussianFactor(:x1, [0.25], 1.0, [0.15]; label = "prior_x1"),
            GaussianFactor(:x2, [1.0, -0.5], I2, [0.2, 0.2]; label = "prior_x2"),

            GaussianFactor(
                :x1, :x2,
                [0.4, -0.2],
                [
                    1.0   0.5  -0.3
                   -0.4   1.2   0.8
                ],
                [0.3, 0.25];
                label = "factor_x1_x2"
            ),

            GaussianFactor(
                :x3, :x4,
                [1.0, -0.5, 0.2, 0.4],
                [
                    1.0   0.2  -0.1   0.0   0.5  -0.3   0.1
                    0.0   1.1   0.4  -0.2  -0.6   0.7   0.3
                   -0.3   0.0   0.8   1.2   0.2   0.1  -0.5
                    0.2  -0.4   0.1   1.0   0.3  -0.2   0.6
                ],
                [0.4, 0.35, 0.45, 0.5];
                label = "factor_x3_x4"
            ),

            GaussianFactor(
                :x1, :x2, :x3,
                [0.1, -0.3, 0.7, -0.2],
                [
                    1.0   0.2  -0.1   0.5   0.0  -0.4   0.3
                   -0.2   1.0   0.6   0.0  -0.5   0.8  -0.1
                    0.4  -0.3   1.2  -0.6   0.7   0.2   0.5
                    0.3   0.4  -0.2   0.1  -0.3   0.6   1.1
                ],
                [0.5, 0.45, 0.55, 0.6];
                label = "factor_x1_x2_x3"
            ),

            GaussianFactor(
                :x1, :x2, :x3, :x4,
                [0.2, -0.1, 0.5, -0.4],
                [
                    1.0   0.2  -0.1   0.3   0.0  -0.2   0.4   0.5  -0.3   0.1
                   -0.3   1.1   0.5   0.0  -0.4   0.7  -0.2  -0.1   0.8   0.2
                    0.2  -0.6   1.0  -0.5   0.3   0.1   0.9   0.4   0.0  -0.7
                   -0.1   0.4   0.2   0.8  -0.6   0.5   0.0  -0.3   1.2   0.6
                ],
                [0.6, 0.5, 0.55, 0.65];
                label = "factor_x1_x2_x3_x4"
            )
        ]

        return GaussianFactorGraph(variables, factors)
    end
end

if !isdefined(@__MODULE__, :gaussianTreeTestGraph)
    function gaussianTreeTestGraph()
        variables = [
            GaussianVariable(:x1, 1; label = "x1"),
            GaussianVariable(:x2, 2; label = "x2"),
            GaussianVariable(:x3, 1; label = "x3")
        ]

        factors = [
            GaussianFactor(:x1, [0.25], 1.0, [0.15]; label = "prior_x1"),
            GaussianFactor(
                :x1, :x2,
                [0.4, -0.2],
                [
                    1.0   0.5  -0.3
                   -0.4   1.2   0.8
                ],
                [0.3, 0.25];
                label = "factor_x1_x2"
            ),
            GaussianFactor(
                :x2, :x3,
                [0.1, -0.3],
                [
                    0.8  -0.2   1.0
                    0.1   0.7  -0.4
                ],
                [0.5, 0.45];
                label = "factor_x2_x3"
            )
        ]

        return GaussianFactorGraph(variables, factors)
    end
end
