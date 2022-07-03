struct WeightedLeastSquares
    estimate::Array{Float64,1}
    rmse::Float64
    mae::Float64
    wrss::Float64
end

struct ErrorMetric
    rmse::Float64
    mae::Float64
    wrss::Float64
end

struct ErrorMetricWiden
    rmse::Float64
    mae::Float64
    wrss::Float64
    rmseGBPWLS::Float64
    maeGBPWLS::Float64
end

########## Compute the GBP error metrics ##########
function errorMetric(gbp::Union{ContinuousModel, ContinuousTreeModel})
    observationGBP = gbp.system.coefficient * gbp.inference.mean
    NactiveRows = activeRows(gbp)

    rmse = 0.0; mae = 0.0; wrss = 0.0
    @inbounds for i = 1:length(gbp.system.observation)
        rmse += (gbp.system.observation[i] - observationGBP[i])^2
        mae += abs(gbp.system.observation[i] - observationGBP[i])
        wrss += (gbp.system.observation[i] - observationGBP[i])^2 / gbp.system.variance[i]
    end
    rmse = (rmse / NactiveRows)^(1/2)
    mae = mae / NactiveRows

    return ErrorMetric(rmse, mae, wrss)
end

########## Compute the GBP error metrics and error according to the WLS ##########
@inline function errorMetric(gbp::Union{ContinuousModel, ContinuousTreeModel}, wls::WeightedLeastSquares)
    observationGBP = gbp.system.coefficient * gbp.inference.mean
    NactiveRows = activeRows(gbp)

    rmse = 0.0; mae = 0.0; wrss = 0.0
    @inbounds for i = 1:length(gbp.system.observation)
        rmse += (gbp.system.observation[i] - observationGBP[i])^2
        mae += abs(gbp.system.observation[i] - observationGBP[i])
        wrss += (gbp.system.observation[i] - observationGBP[i])^2 / gbp.system.variance[i]
    end
    rmse = (rmse / NactiveRows)^(1/2)
    mae = mae / NactiveRows

    rmseGBPWLS = 0.0; maeGBPWLS = 0.0
    @inbounds for i = 1:gbp.graph.Nvariable
        rmseGBPWLS += (wls.estimate[i] - gbp.inference.mean[i])^2
        maeGBPWLS += abs(wls.estimate[i] - gbp.inference.mean[i])
    end
    rmseGBPWLS = (rmseGBPWLS / gbp.graph.Nvariable)^(1/2)
    maeGBPWLS = maeGBPWLS / gbp.graph.Nvariable

    return ErrorMetricWiden(rmse, mae, wrss, rmseGBPWLS, maeGBPWLS)
end

########## Compute WLS solution and error metrics ##########
function wls(gbp::Union{ContinuousModel, ContinuousTreeModel})
    W = spdiagm(0 =>  @. 1.0 / sqrt(gbp.system.variance))
    A = W * gbp.system.coefficient
    G = A' * A
    b = A' * W * gbp.system.observation
    x = G \ b

    observationWLS = gbp.system.coefficient * x
    NactiveRows = activeRows(gbp)

    rmse = 0.0; mae = 0.0; wrss = 0.0
    @inbounds for i = 1:length(gbp.system.observation)
        rmse += (gbp.system.observation[i] - observationWLS[i])^2
        mae += abs(gbp.system.observation[i] - observationWLS[i])
        wrss += (gbp.system.observation[i] - observationWLS[i])^2 / gbp.system.variance[i]
    end
    rmse = (rmse / NactiveRows)^(1/2)
    mae = mae / NactiveRows

    return WeightedLeastSquares(x, rmse, mae, wrss)
end

########## Find number of nonzeros rows ##########
@inline function activeRows(gbp)
    NactiveRows = length(gbp.system.observation)
    @inbounds for (k, i) in enumerate(gbp.system.observation)
        if i == 0.0
            if all(gbp.system.coefficientTranspose[:, k] .== 0.0)
                NactiveRows -= 1
            end
        end
    end

    return NactiveRows
end

########## Check keyword arguments ##########
function checkKeywords(prob, alpha, variance)
    #### Check convergence parameters
    if prob <= 0.0 || prob >= 1.0
        error("Invalid prob value.")
    end
    if alpha <= 0.0 || alpha >= 1.0
        error("Invalid alpha value.")
    end

    #### Check initial variance
    if variance < 0.0
        error("Invalid variance value.")
    end
end