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
    dynamicFactor!(gbp; factor = 1, mean = 6, variance = 1)
    dynamicFactor!(gbp; factor = 3, mean = 4, variance = 1)
    for iteration = 10:99
        messageDampFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    dynamicFactor!(gbp; factor = 2, mean = 3, variance = 1)
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
    dynamicFactor!(gbp; factor = 1, mean = 6, variance = 1)
    dynamicFactor!(gbp; factor = 3, mean = 4, variance = 1)
    dynamicFactor!(gbp; factor = 2, mean = 3, variance = 1)
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
    for iteration = 1:900
        ageingVariance!(gbp; factor = 4, initial = 1, limit = 1e60, model = 1, a = 1e57, tau = iteration)
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
    for iteration = 1:1000
        ageingVariance!(gbp; factor = 4, initial = 1, limit = 1e60, model = 2, a = 1e57, b = 0.00002, tau = iteration)
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
    for iteration = 1:900
        ageingVariance!(gbp; factor = 4, initial = 1, limit = 1e60, model = 3, a = 0.08, b = 2, tau = iteration)
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    @test round.(gbp.inference.mean, digits = 3) ≈ [1.0; 2.0; 3.0]
end

@testset "GraphicalModel" begin
    ### Freeze factor node
    gbp = graphicalModel("data33_14.h5")
    for i = 1:2
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    T = gbp.inference
    beforeFreeze = [T.fromFactor T.toVariable T.meanFactorVariable T.varianceFactorVariable]
    idx = T.fromFactor .== 10
    beforeFreeze = beforeFreeze[idx, :]
    freezeFactor!(gbp; factor = 10)
    for i = 1:20
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    T = gbp.inference
    afterFreeze = [T.fromFactor T.toVariable T.meanFactorVariable T.varianceFactorVariable]
    afterFreeze = afterFreeze[idx, :]
    @test sum(beforeFreeze[:, 3] - afterFreeze[:, 3]) == 0
    @test sum(beforeFreeze[:, 4] - afterFreeze[:, 4]) == 0

    ### Defreeze factor node
    defreezeFactor!(gbp; factor = 10)
    for i = 1:20
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    T = gbp.inference
    afterDefreeze = [T.fromFactor T.toVariable T.meanFactorVariable T.varianceFactorVariable]
    afterDefreeze = afterDefreeze[idx, :]
    @test sum(afterDefreeze[:, 3] - afterFreeze[:, 3]) != 0
    @test sum(afterDefreeze[:, 4] - afterFreeze[:, 4]) != 0

    ### Freeze variable node
    gbp = graphicalModel("data33_14.h5")
    for i = 1:2
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    T = gbp.inference
    beforeFreeze = [T.fromVariable T.toFactor T.meanVariableFactor T.varianceVariableFactor]
    idx = T.fromVariable .== 9
    beforeFreeze = beforeFreeze[idx, :]
    freezeVariable!(gbp; variable = 9)
    for i = 1:20
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    T = gbp.inference
    afterFreeze = [T.fromVariable T.toFactor T.meanVariableFactor T.varianceVariableFactor]
    afterFreeze = afterFreeze[idx, :]
    @test sum(beforeFreeze[:, 3] - afterFreeze[:, 3]) == 0
    @test sum(beforeFreeze[:, 4] - afterFreeze[:, 4]) == 0

    ### Defreeze variable node
    defreezeVariable!(gbp; variable = 9)
    for i = 1:20
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    T = gbp.inference
    afterDefreeze = [T.fromFactor T.toVariable T.meanFactorVariable T.varianceFactorVariable]
    afterDefreeze = afterDefreeze[idx, :]
    @test sum(afterDefreeze[:, 3] - afterFreeze[:, 3]) != 0
    @test sum(afterDefreeze[:, 4] - afterFreeze[:, 4]) != 0

    ### Freeze and defreeze message variable node to factor node
    gbp = graphicalModel("data33_14.h5")
    for iteration = 1:5
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    beforeFreeze6 = copy(gbp.inference.meanVariableFactor[6])
    beforeFreeze57 = copy(gbp.inference.meanVariableFactor[57])
    freezeVariableFactor!(gbp; variable = 9, factor = 29)
    freezeVariableFactor!(gbp; variable = 2, factor = 10)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    afterFreeze6 = copy(gbp.inference.meanVariableFactor[6])
    afterFreeze57 = copy(gbp.inference.meanVariableFactor[57])
    @test beforeFreeze6 == afterFreeze6
    @test beforeFreeze57 == afterFreeze57
    defreezeVariableFactor!(gbp; variable = 9, factor = 29)
    defreezeVariableFactor!(gbp; variable = 2, factor = 10)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    @test afterFreeze6 != gbp.inference.meanVariableFactor[6]
    @test afterFreeze57 != gbp.inference.meanVariableFactor[57]

    ### Freeze and defreeze message factor node to variable node
    gbp = graphicalModel("data33_14.h5")
    for iteration = 1:5
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    beforeFreeze6 = copy(gbp.inference.meanFactorVariable[6])
    beforeFreeze57 = copy(gbp.inference.meanFactorVariable[57])
    freezeFactorVariable!(gbp; variable = 2, factor = 11)
    freezeFactorVariable!(gbp; variable = 14, factor = 25)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    afterFreeze6 = copy(gbp.inference.meanFactorVariable[6])
    afterFreeze57 = copy(gbp.inference.meanFactorVariable[57])
    @test beforeFreeze6 == afterFreeze6
    @test beforeFreeze57 == afterFreeze57
    defreezeFactorVariable!(gbp; variable = 2, factor = 11)
    defreezeFactorVariable!(gbp; variable = 14, factor = 25)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    @test afterFreeze6 != gbp.inference.meanFactorVariable[6]
    @test afterFreeze57 != gbp.inference.meanFactorVariable[57]

    ### Hide factor node
    H = [1.5 0.0 1.0;
         2.0 3.1 4.6;
         2.0 1.1 0.6;
         1.0 0.0 0.0;
         1.0 0.0 0.0]
    z = [0.5; 0.8; 4.1; 5.1; 10]
    v = [0.1; 1.0; 1.0; 0.0002;  0.0001]
    gbp = graphicalModel(H,z,v)
    hideFactor!(gbp; factor = 2)
    for iteration = 1:200
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    @test round.(gbp.inference.mean, digits = 4) ≈ round.(exact.estimate, digits = 4)

    H = [1.5 0.0 1.0;
         2.0 3.1 4.6;
         2.0 1.1 0.6;
         1.0 0.0 0.0;
         1.0 0.0 0.0]
    z = [0.5; 0.8; 4.1; 5.1; 10]
    v = [0.1; 1.0; 1.0; 0.0002;  0.0001]
    gbp = graphicalModel(H,z,v)
    hideFactor!(gbp; factor = 4)
    for iteration = 1:200
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    @test round.(gbp.inference.mean, digits = 4) ≈ round.(exact.estimate, digits = 4)

    ### Add factor node
    gbp = graphicalModel("data33_14.h5")
    for iteration = 1:15
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    mean = [2.2; 3.1; 0.5]
    variance = [0.001; 0.001; 0.001]
    jacobian = zeros(3, 14)
    jacobian[1, 5] = 2.0
    jacobian[2, 1] = 0.45; jacobian[2, 5] = 0.23; jacobian[2, 14] = 0.1
    jacobian[3, 3] = 0.25; jacobian[3, 6] = 0.8
    addFactors!(gbp; mean = mean, variance = variance, jacobian = jacobian)
    for iteration = 1:200
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    @test round.(gbp.inference.mean, digits = 4) ≈ round.(exact.estimate, digits = 4)

    ### Hide and add factor node
    gbp = graphicalModel("data33_14.h5")
    for iteration = 1:15
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    hideFactor!(gbp; factor = 2)
    mean = [2.2; 3.1; 0.5]
    variance = [0.001; 0.001; 0.001]
    jacobian = zeros(3, 14)
    jacobian[1, 5] = 2.0
    jacobian[2, 1] = 0.45; jacobian[2, 5] = 0.23; jacobian[2, 14] = 0.1
    jacobian[3, 3] = 0.25; jacobian[3, 6] = 0.8
    addFactors!(gbp; mean = mean, variance = variance, jacobian = jacobian)
    for iteration = 1:200
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    exact = wls(gbp)
    @test round.(gbp.inference.mean, digits = 4) ≈ round.(exact.estimate, digits = 4)
end


@testset "GraphicalModelTree" begin
    H = [2 3 0 0 0 0 0 0 0 0 0;
         0 4 6 8 2 3 0 0 0 0 0;
         0 0 0 5 0 0 0 0 0 0 0;
         0 0 0 0 2 0 0 0 0 0 0;
         0 0 0 0 0 8 0 2 0 1 0;
         0 0 0 0 0 4 0 0 0 0 0;
         0 0 0 0 0 3 7 0 2 0 0;
         0 0 0 0 0 0 0 5 0 0 0;
         0 0 0 0 0 0 0 3 0 0 0;
         2 0 0 0 0 0 0 0 0 0 0;
         2 0 0 0 0 0 0 0 0 0 0;
         1 0 0 0 0 0 0 0 0 0 2]
    z = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12]
    v = [3; 4; 2; 5; 1; 8; 9; 1; 2; 3; 8; 4]

    ### isTree
    gbp = graphicalModel(H, z, v)
    tree = isTree(gbp)
    @test tree == true

    gbp = graphicalModel("data33_14.h5")
    tree = isTree(gbp)
    @test tree == false

    # Forward-backward algorithm
    gbp = graphicalModel(H, z, v)
    for iteration = 1:100
        messageFactorVariableVanilla(gbp)
        messageVariableFactorVanilla(gbp)
    end
    marginal(gbp)
    mean = copy(gbp.inference.mean)

    gbp = graphicalModelTree(H, z, v; root = 1)
    while gbp.graph.forward
        forwardVariableFactor(gbp)
        forwardFactorVariable(gbp)
    end
    while gbp.graph.backward
        backwardVariableFactor(gbp)
        backwardFactorVariable(gbp)
    end
    marginal(gbp)
    @test sum(round.(mean - gbp.inference.mean, digits = 4)) == 0

    # Forward-backward algorithm - new root
    gbp = graphicalModelTree(H, z, v; root = 8)
    while gbp.graph.forward
        forwardVariableFactor(gbp)
        forwardFactorVariable(gbp)
    end
    while gbp.graph.backward
        backwardVariableFactor(gbp)
        backwardFactorVariable(gbp)
    end
    marginal(gbp)
    @test sum(round.(mean - gbp.inference.mean, digits = 4)) == 0

    # Forward-backward algorithm - new root
    gbp = graphicalModelTree(H, z, v; root = 2)
    while gbp.graph.forward
        forwardVariableFactor(gbp)
        forwardFactorVariable(gbp)
    end
    while gbp.graph.backward
        backwardVariableFactor(gbp)
        backwardFactorVariable(gbp)
    end
    marginal(gbp)
    @test sum(round.(mean - gbp.inference.mean, digits = 4)) == 0
end