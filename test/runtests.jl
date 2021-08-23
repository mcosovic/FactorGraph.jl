using GaussBP
using SparseArrays
using Test

@testset "StaticGBP" begin
    ### Vanilla GBP with Damping and XLSX input
    gbp = graphicalModel("data33_14.xlsx"; variance = 1e60)
    for iteration = 1:50
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    for iteration = 51:1000
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    displayData(gbp)
    @test maximum(gbp.inference.mean ./ exact.estimate) < 1.0000000009

    ### Kahan GBP with Damping and HDF5 input
    gbp = graphicalModel("data33_14.h5"; variance = 1e60)
    for iteration = 1:50
        messageFactorVariableKahan(gbp)
        messageVariableFactorKahan(gbp)
    end
    for iteration = 51:1000
        messageDampFactorVariableKahan(gbp)
        messageVariableFactorKahan(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    displayData(gbp)
    @test maximum(gbp.inference.mean ./ exact.estimate) < 1.0000000009

    ### Efficient GBP with Damping and HDF5 input
    gbp = graphicalModel("data33_14.h5"; variance = 1e10)
    for iteration = 1:50
        messageFactorVariableEfficient(gbp)
        messageVariableFactorEfficient(gbp)
    end
    for iteration = 51:1000
        messageDampFactorVariableEfficient(gbp)
        messageVariableFactorEfficient(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    @test maximum(gbp.inference.mean ./ exact.estimate) < 1.00000000003

    ### Vanilla GBP with Damping with marginal computation in each iteration
    gbp = graphicalModel("data33_14.h5"; variance = 1e60)
    for iteration = 1:50
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
        marginal(gbp)
    end
    for iteration = 51:1000
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
        marginal(gbp)
    end
    exact = wls(gbp)
    displayData(gbp, exact)
    @test maximum(gbp.inference.mean ./ exact.estimate) < 1.0000000009

    ### Vanilla GBP with Damping with error metrics
    gbp = graphicalModel("data33_14.h5"; variance = 1e60)
    for iteration = 1:50
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    for iteration = 51:500
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    errors = errorMetric(gbp, exact)
    displayData(gbp, exact, errors)
    @test abs(errors.rmse[1] - exact.rmse[1]) < 1e-10
    @test abs(errors.mae[1] - exact.mae[1]) < 1e-10
    @test abs(errors.wrss[1] - exact.wrss[1]) < 1e-10
    @test errors.rmseGBPWLS[1] < 1e-10
    @test errors.maeGBPWLS[1] < 1e-10

    ### Vanilla GBP with Damping and passing arguments
    H = [3.0 2.0 -1.0; -2.0 2.0 1.0; 1.0 1.0 1.0]
    z = [6.0; 3.0; 4.0]
    v = [1.0; 1.0; 1.0]
    gbp = graphicalModel(H, z, v)
    for iteration = 1:5
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    for iteration = 6:1000
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    @test round.(gbp.inference.mean, digits = 3) ≈ [1.0; 2.0; 1.0]

    ### Kahan GBP and passing arguments
    H = [1.0 0.0 0.0; 25 -25 0; -50 -40 90; 0 1 0]
    z = [0; 1.795; 1.966; -0.066]
    v = [1e-2; 1e-2; 1e-2; 1e-6]

    gbp = graphicalModel(H, z, v, variance = 1e60)
    for iteration = 1:30
        messageFactorVariableKahan(gbp)
        messageVariableFactorKahan(gbp)
    end
    marginal(gbp)
    @test round.(gbp.inference.mean, digits = 5) ≈ [0.00579; -0.066; -0.00427]
end

@testset "DynamicGBP" begin
    H = [3.0 2.0 -1.0; -2.0 2.0 1.0; 1.0 1.0 1.0]
    z = [16.0; 13.0; 14.0]
    v = [10.0; 20.0; 30.0]

    ### Vanilla GBP with Damping and two dynamic updates
    gbp = graphicalModel(H, z, v)
    for iteration = 1:9
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    dynamic = [1 6.0 1.0; 3 4.0 1.0]
    dynamicInference(gbp, dynamic)
    for iteration = 10:99
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    dynamic = [2 3.0 1.0]
    dynamicInference(gbp, dynamic)
    for iteration = 100:1000
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    @test round.(gbp.inference.mean, digits = 3) ≈ [1.0; 2.0; 1.0]

    ### Vanilla GBP with Damping and one dynamic updates
    gbp = graphicalModel(H, z, v)
    for iteration = 1:9
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
        marginal(gbp)
    end
    dynamic = [1 6.0 1.0; 3 4.0 1.0; 2 3.0 1.0]
    dynamicInference(gbp, dynamic)
    for iteration = 10:1000
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
        marginal(gbp)
    end
    @test round.(gbp.inference.mean, digits = 3) ≈ [1.0; 2.0; 1.0]
end

@testset "AgeingGBP" begin
    H = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 3.0 2.0 -1.0]
    z = [1.0; 2.0; 3.0; 11.0]
    v = [1.0; 1.0; 1.0; 1.0]

    ### Vanilla GBP with linear ageing
    gbp = graphicalModel(H, z, v)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    dynamic = [4 6.0 1.0 1 1e57 2 1e60]
    dynamicInference(gbp, dynamic)
    for iteration = 101:1000
        ageingInference(gbp, dynamic)
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    @test round.(gbp.inference.mean, digits = 3) ≈ [1.0; 2.0; 3.0]

    ### Vanilla GBP with logarithmic ageing
    gbp = graphicalModel(H, z, v)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    dynamic = [4 6.0 1.0 2 1e57 0.00002 1e60]
    dynamicInference(gbp, dynamic)
    for iteration = 101:1000
        ageingInference(gbp, dynamic)
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    @test round.(gbp.inference.mean, digits = 3) ≈ [1.0; 2.0; 3.0]

    ### Vanilla GBP with exponential ageing
    gbp = graphicalModel(H, z, v)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    dynamic =  [4 6.0 1.0 3 0.08 2 1e60]
    dynamicInference(gbp, dynamic)
    for iteration = 101:1000
        ageingInference(gbp, dynamic)
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    @test round.(gbp.inference.mean, digits = 3) ≈ [1.0; 2.0; 3.0]
end