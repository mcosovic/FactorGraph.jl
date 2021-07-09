using GaussBP
using SparseArrays
using Test

@testset "GaussBP" begin
    results, system = gbp("data33_14.csv"; max = 1000, damp = 50, bump = 1000, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, algorithm = "gbp", out = ["no", "error", "wls", "display"])
    @test maximum(results.means ./ results.meansWLS) < 1.0000009

    # Xbp, system = bp("data33_14.csv"; max = 1000, damp = 50, bump = 1000, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, method = "passing", algorithm = "kahan", wls = "builtin")
    # @test maximum(Xbp ./ Xwls) < 1.00000000000003

    # Xbp, system = bp("data33_14.csv"; max = 1000, damp = 50, bump = 100, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, method = "passing", algorithm = "kahan", wls = "builtin")
    # @test maximum(Xbp ./ Xwls) < 1.00000000000003

    # Xbp, system = bp("data33_14.csv"; max = 1000, damp = 1000, bump = 100, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, method = "passing", algorithm = "kahan", wls = "builtin")
    # @test maximum(Xbp ./ Xwls) < 1.00000000000003

    # Xbp, system = bp("data33_14.csv"; max = 1000, damp = 50, bump = 100, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, method = "recursion", algorithm = "sum", wls = "builtin")
    # @test maximum(Xbp ./ Xwls) < 1.0000009

    # H = [1.5 0.0 2.0;
    #      0.0 3.1 4.6;
    #      2.6 8.1 0.5]

    # b = [0.8; 4.1; 2.2];
    # v = [1.0; 1.0; 1.0]     


    # results, system = gbp(H, b, v; max = 300, out = ["error", "wls", "display"], algorithm = "kahan")

end
