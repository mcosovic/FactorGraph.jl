using GaussBP
using SparseArrays
using Test

@testset "StaticGBP" begin
    results, system = gbp("data33_14.xlsx"; max = 1000, damp = 50, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, algorithm = "vanilla", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.0000000009

    results, system = gbp("data33_14.h5"; max = 1000, damp = 50, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, algorithm = "kahan", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.00000000000003

    results, system = gbp("data33_14.h5"; max = 1000, damp = 50, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, algorithm = "efficient", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.0000000009

    results, system = gbp("data33_14.h5"; max = 1000, damp = 50, bump = 500, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e5, algorithm = "efficient", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.0000000009

    results, system = gbp("data33_14.h5"; max = 500, damp = 50, algorithm = "vanilla", out = ["iterate", "wls"])
    @test maximum(results.gbp.mean[:, 500] ./ results.wls.mean) < 1.0000000009

    results, system = gbp("data33_14.h5"; max = 500, damp = 50, algorithm = "vanilla", out = ["wls", "error"])
    @test abs(results.gbp.rmse[1] - results.wls.rmse[1]) < 1e-10
    @test abs(results.gbp.mae[1] - results.wls.mae[1]) < 1e-10
    @test abs(results.gbp.wrss[1] - results.wls.wrss[1]) < 1e-10
    @test results.wls.rmseGBPWLS[1] < 1e-10
    @test results.wls.maeGBPWLS[1] < 1e-10

    H = [3.0 2.0 -1.0; -2.0 2.0 1.0; 1.0 1.0 1.0]
    z = [6.0; 3.0; 4.0] 
    v = [1.0; 1.0; 1.0]    

    results, system = gbp(H, z, v; max = 1000, damp = 50, algorithm = "vanilla", out = "wls")
    @test results.gbp.mean ≈ [1.0; 2.0; 1.0]
end

@testset "DynamicGBP" begin
    H = [3.0 2.0 -1.0; -2.0 2.0 1.0; 1.0 1.0  1.0]
    z = [16.0; 13.0; 14.0] 
    v = [10.0; 20.0; 30.0]    
    d = [10 1 6.0 1.0; 200 3 4.0 1.0; 100 2 3.0 1.0]

    results, system = gbp(H, z, v, d; max = 1000, damp = 50, algorithm = "vanillaDynamic", out = ["wls", "display"])
    @test results.gbp.mean[:, end] ≈ [1.0; 2.0; 1.0]

    results, system = gbp(H, z, v, d; max = 1000, damp = 50, algorithm = "efficientDynamic")
    @test results.gbp.mean[:, end] ≈ [1.0; 2.0; 1.0]

    results, system = gbp(H, z, v, d; max = 1000, damp = 50, algorithm = "kahanDynamic")
    @test results.gbp.mean[:, end] ≈ [1.0; 2.0; 1.0]

    results, system = gbp("dataDynamic33_14.xlsx"; max = 1000, algorithm = "kahanDynamic", out = ["wls", "error"])
    @test maximum(abs.(results.gbp.rmse - results.wls.rmse)) < 1e-8
    @test maximum(abs.(results.gbp.mae - results.wls.mae)) < 1e-8
    @test maximum(abs.(results.gbp.wrss - results.wls.wrss)) < 1e-8
    @test maximum(results.wls.rmseGBPWLS) < 1e-8
    @test maximum(results.wls.maeGBPWLS) < 1e-8

    results, system = gbp("dataDynamic33_14.h5"; max = 1000, algorithm = "kahanDynamic", out = ["wls", "error"])
    @test maximum(abs.(results.gbp.rmse[end] - results.wls.rmse[end])) < 1e-8

    results, system = gbp("dataDynamic33_14.h5"; max = 1000, algorithm = "kahanDynamic", out = ["wls", "error", "iterate"])
    @test maximum(abs.(results.gbp.rmse[end] - results.wls.rmse[end])) < 1e-8
end