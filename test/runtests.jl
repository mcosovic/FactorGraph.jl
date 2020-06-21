using GaussianBP
using SparseArrays
using Test

@testset "GaussianBP" begin
    Xbp, system = bp("data33_14.csv"; max = 1000, damp = 50, bump = 1000, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, method = "passing", algorithm = "sum", wls = "builtin")
    W = spdiagm(0 =>  @. 1.0 / sqrt(system.v))
    H = W * system.J
    Xwls = (H' * H) \ (H' * W * system.z)
    @test maximum(Xbp ./ Xwls) < 1.00000009

    Xbp, system = bp("data33_14.csv"; max = 1000, damp = 50, bump = 1000, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, method = "passing", algorithm = "kahan", wls = "builtin")
    @test maximum(Xbp ./ Xwls) < 1.000000000000003

    Xbp, system = bp("data33_14.csv"; max = 1000, damp = 50, bump = 100, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, method = "passing", algorithm = "kahan", wls = "builtin")
    @test maximum(Xbp ./ Xwls) < 1.000000000000003

    Xbp, system = bp("data33_14.csv"; max = 1000, damp = 1000, bump = 100, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, method = "passing", algorithm = "kahan", wls = "builtin")
    @test maximum(Xbp ./ Xwls) < 1.000000000000003

    Xbp, system = bp("data33_14.csv"; max = 1000, damp = 50, bump = 100, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, method = "recursion", algorithm = "sum", wls = "builtin")
    @test maximum(Xbp ./ Xwls) < 1.00000009
end
