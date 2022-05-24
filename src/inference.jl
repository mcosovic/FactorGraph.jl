########## Compute the GBP marginal vector ##########
function marginal(gbp::ContinuousModel)
    @inbounds for i in gbp.graph.iterateMarginal
        Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]

        for j in gbp.graph.colptrMarginal[i]
            Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
            Wcol += 1 / gbp.inference.varianceFactorVariable[j]
        end
        gbp.inference.variance[i] = 1 / Wcol
        gbp.inference.mean[i] = Mcol * gbp.inference.variance[i]
    end
end

########## Compute the GBP marginal vector for the tree factor graph ##########
function marginal(gbp::ContinuousTreeModel)
    @inbounds for i = 1:gbp.graph.Nvariable
        Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]

        for j in gbp.graph.incomingToVariable[i]
            Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
            Wcol += 1 / gbp.inference.varianceFactorVariable[j]
        end
        gbp.inference.variance[i] = 1 / Wcol
        gbp.inference.mean[i] = Mcol * gbp.inference.variance[i]
    end
end

########## Compute the BP marginal vector for the tree factor graph ##########
function marginal(bp::DiscreteTreeModel)
    @inbounds for i = 1:bp.graph.Nvariable
        bp.inference.marginal[i] = copy(bp.graph.messageDirect[i])
        for k in bp.graph.incomingToVariable[i]
            bp.inference.marginal[i] .*= bp.inference.messageFactorVariable[k]
        end
        normalizeConstant = sum(bp.inference.marginal[i])
        bp.inference.marginal[i] = bp.inference.marginal[i] ./ normalizeConstant
    end
end

########## Compute the BP unnormalized marginal vector for the tree factor graph ##########
function marginalUnnormalized(bp::DiscreteTreeModel)
    @inbounds for i = 1:bp.graph.Nvariable
        bp.inference.marginal[i] = copy(bp.graph.messageDirect[i])
        for k in bp.graph.incomingToVariable[i]
            bp.inference.marginal[i] .*= bp.inference.messageFactorVariable[k]
        end
    end
end

########## Dynamic the GBP update ##########
@inline function dynamicFactor!(gbp::ContinuousModel; factor = 0::Int64, mean = 0, variance = 0)
    if (gbp.system.jacobianTranspose.colptr[factor + 1] - gbp.system.jacobianTranspose.colptr[factor]) == 1
        idx = gbp.system.jacobianTranspose.colptr[factor]
        variable = trunc(Int, gbp.system.jacobianTranspose.rowval[idx])

        gbp.graph.meanDirect[variable] -= gbp.system.observation[factor] * gbp.system.jacobianTranspose[variable, factor] / gbp.system.variance[factor]
        gbp.graph.weightDirect[variable] -= gbp.system.jacobianTranspose[variable, factor]^2 / gbp.system.variance[factor]

        gbp.graph.meanDirect[variable] += mean * gbp.system.jacobianTranspose[variable, factor] / variance
        gbp.graph.weightDirect[variable] += gbp.system.jacobianTranspose[variable, factor]^2 / variance
    else
        k = gbp.graph.dynamic[factor]
        @inbounds for i in gbp.graph.rowptr[k]
            gbp.graph.meanIndirect[i] = mean
            gbp.graph.varianceIndirect[i] = variance
        end
    end
    gbp.system.observation[factor] = mean
    gbp.system.variance[factor] = variance
end

######### Ageing the GBP update ##########
@inline function ageingVariance!(gbp::ContinuousModel; factor = 0::Int64, initial = 0, limit = 0, model = 0::Int64, a = 0, b = 0, tau = 0)
    if gbp.system.variance[factor] < limit
        Nvariable = gbp.system.jacobianTranspose.colptr[factor + 1] - gbp.system.jacobianTranspose.colptr[factor]
        if Nvariable == 1
            idx = gbp.system.jacobianTranspose.colptr[factor]
            variable = trunc(Int, gbp.system.jacobianTranspose.rowval[idx])

            gbp.graph.meanDirect[variable] -= gbp.system.observation[factor] * gbp.system.jacobianTranspose[variable, factor] / gbp.system.variance[factor]
            gbp.graph.weightDirect[variable] -= gbp.system.jacobianTranspose[variable, factor]^2 / gbp.system.variance[factor]
        end

        if model == 1
            gbp.system.variance[factor] = a * tau + initial
        elseif model == 2
            gbp.system.variance[factor] = a * log10((tau + 1 + b) / (1 + b)) + initial
        elseif model == 3
            gbp.system.variance[factor] = initial * (1 + b)^(a * tau)
        end

        if gbp.system.variance[factor] > limit
            gbp.system.variance[factor] =  limit
        end

        if Nvariable == 1
            gbp.graph.meanDirect[variable] += gbp.system.observation[factor] * gbp.system.jacobianTranspose[variable, factor] / gbp.system.variance[factor]
            gbp.graph.weightDirect[variable] += gbp.system.jacobianTranspose[variable, factor]^2 / gbp.system.variance[factor]
        else
            k = gbp.graph.dynamic[factor]
            @inbounds for i in gbp.graph.rowptr[k]
                gbp.graph.varianceIndirect[i] = gbp.system.variance[factor]
            end
        end
    end
end

########## Set damping parameters ##########
function damping!(gbp::ContinuousModel; prob::Float64 = 0.6, alpha::Float64 = 0.4)
    gbp.graph.alphaNew = fill(1.0, gbp.graph.Nlink)
    gbp.graph.alphaOld = fill(0.0, gbp.graph.Nlink)
    bernoulliSample = randsubseq(collect(1:gbp.graph.Nlink), prob)
    @inbounds for i in bernoulliSample
        gbp.graph.alphaNew[i] = 1.0 - alpha
        gbp.graph.alphaOld[i] = alpha
    end
end