using GaussBP
using SparseArrays
using Test

@testset "GaussBP" begin
    results, system = gbp("data33_14.csv"; max = 1000, damp = 50, bump = 1000, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, algorithm = "gbp", out = ["no", "error", "wls", "display"])
    @test maximum(results.means ./ results.meansWLS) < 1.0000009

    results, system = gbp("data33_14.csv"; max = 1000, damp = 50, bump = 1000, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, algorithm = "kahan", out = ["no", "error", "wls", "display"])
    @test maximum(results.means ./ results.meansWLS) < 1.00000000000003

    results, system = gbp("data33_14.csv"; max = 1000, damp = 50, bump = 1000, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, algorithm = "efficient", out = ["no", "error", "wls", "display"])
    @test maximum(results.means ./ results.meansWLS) < 1.0000009
end
