using GaussBP
using SparseArrays
using Test

@testset "StaticGBP" begin
    results, system = gbp("data33_14.xlsx"; max = 1000, damp = 50, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, algorithm = "vanilla", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.0000000009

    results, system = gbp("data33_14.h5"; max = 1000, damp = 50, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e60, algorithm = "kahan", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.0000000000003

    results, system = gbp("data33_14.h5"; max = 1000, damp = 50, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e10, algorithm = "efficient", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.00000000003

    results, system = gbp("data33_14.h5"; max = 1000, damp = 50, bump = 500, prob = 0.6, alpha = 0.4, mean = 0.0, variance = 1e10, algorithm = "efficient", out = ["error", "wls", "display"])
    @test maximum(results.gbp.mean ./ results.wls.mean) < 1.00000000003

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
    results, system = gbp(H, z, v; max = 2000, damp = 5, algorithm = "vanilla", out = "wls")
    @test round.(results.gbp.mean, digits=3) ≈ [1.0; 2.0; 1.0]

    H = [1.0 0.0 0.0; 25 -25 0; -50 -40 90; 0 1 0]
    z = [0; 1.795; 1.966; -0.066]
    v = [1e-2; 1e-2; 1e-2; 1e-6]     
    results, system = gbp(H, z, v; max = 30, variance = 1e60, algorithm = "kahan", out = ["display", "wls"]) 
    @test round.(results.gbp.mean, digits=5) ≈ [0.00579; -0.066; -0.00427]
end

@testset "DynamicGBP" begin
    H = [3.0 2.0 -1.0; -2.0 2.0 1.0; 1.0 1.0  1.0]
    z = [16.0; 13.0; 14.0] 
    v = [10.0; 20.0; 30.0]    
    d = [10 1 6.0 1.0; 200 3 4.0 1.0; 100 2 3.0 1.0]

    results, system = gbp(H, z, v, d; max = 2000, damp = 5, algorithm = "vanillaDynamic", out = ["wls", "display"])
    @test round.(results.gbp.mean[:, end], digits=3) ≈ [1.0; 2.0; 1.0]

    results, system = gbp(H, z, v, d; max = 2000, damp = 5, algorithm = "efficientDynamic")
    @test round.(results.gbp.mean[:, end], digits=3) ≈ [1.0; 2.0; 1.0]

    results, system = gbp(H, z, v, d; max = 2000, damp = 5, algorithm = "kahanDynamic")
    @test round.(results.gbp.mean[:, end], digits=3) ≈ [1.0; 2.0; 1.0]

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

@testset "AgeingGBP" begin
    H = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 3.0 2.0 -1.0]
    z = [1.0; 2.0; 3.0; 11.0] 
    v = [1.0; 1.0; 1.0; 1.0]    
    
    d = [100 4 6.0 1.0 1 500 1e57 2 1e60]
    results, system = gbp(H, z, v, d; max = 2000, damp = 5, algorithm = "vanillaAgeing", out = "display")     
    @test results.gbp.mean[:, end] ≈ [1.0; 2.0; 3.0]

    d = [100 4 6.0 1.0 2 500 1e57 0.00002 1e60]
    results, system = gbp(H, z, v, d; max = 2000, damp = 5, algorithm = "vanillaAgeing", out = "display")     
    @test results.gbp.mean[:, end] ≈ [1.0; 2.0; 3.0]

    d = [100 4 6.0 1.0 3 100 0.08 2 1e60]
    results, system = gbp(H, z, v, d; max = 2000, damp = 5, algorithm = "vanillaAgeing", out = "display")     
    @test results.gbp.mean[:, end] ≈ [1.0; 2.0; 3.0]
end