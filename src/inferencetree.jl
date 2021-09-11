########## Forward messages from variable to factor nodes ##########
function forwardVariableFactor(gbp::GraphicalModelTree)
    @inbounds for variable in gbp.graph.iterateVariable
        factor = gbp.graph.colForward[variable][1]

        gbp.graph.passVariableFactor += 1
        gbp.inference.fromVariable[gbp.graph.passVariableFactor] = variable
        gbp.inference.toFactor[gbp.graph.passVariableFactor] = factor

        Mcol = gbp.graph.meanDirect[variable]; Wcol = gbp.graph.weightDirect[variable]
        for j in gbp.graph.incomingToVariable[variable]
            Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
            Wcol += 1 / gbp.inference.varianceFactorVariable[j]
        end
        gbp.inference.varianceVariableFactor[gbp.graph.passVariableFactor] = 1 / Wcol
        gbp.inference.meanVariableFactor[gbp.graph.passVariableFactor] = Mcol / Wcol

        push!(gbp.graph.incomingToFactor[factor], gbp.graph.passVariableFactor)
        for (k, variables) in enumerate(gbp.graph.rowForward[factor])
            if variable == variables
                deleteat!(gbp.graph.rowForward[factor], k)
                deleteat!(gbp.graph.colForward[variable], 1)
            end
        end
        if length(gbp.graph.rowForward[factor]) == 1
            push!(gbp.graph.iterateFactor, factor)
        end
    end
    gbp.graph.iterateVariable = Int64[]
end

########## Forward messages from factor to variable nodes ##########
function forwardFactorVariable(gbp::GraphicalModelTree)
    @inbounds for factor in gbp.graph.iterateFactor
        variable = gbp.graph.rowForward[factor][1]

        gbp.graph.passFactorVariable += 1
        gbp.inference.fromFactor[gbp.graph.passFactorVariable] = factor
        gbp.inference.toVariable[gbp.graph.passFactorVariable] = variable

        Mrow = gbp.system.observation[factor]; Vrow = gbp.system.variance[factor]
        for j in gbp.graph.incomingToFactor[factor]
            coeff = gbp.system.jacobian[gbp.inference.toFactor[j], gbp.inference.fromVariable[j]]
            Mrow -= coeff * gbp.inference.meanVariableFactor[j]
            Vrow += coeff^2 * gbp.inference.varianceVariableFactor[j]
        end
        gbp.inference.meanFactorVariable[gbp.graph.passFactorVariable] = Mrow / gbp.system.jacobian[factor, variable]
        gbp.inference.varianceFactorVariable[gbp.graph.passFactorVariable] = Vrow / (gbp.system.jacobian[factor, variable]^2)

        push!(gbp.graph.incomingToVariable[variable], gbp.graph.passFactorVariable)
        for (k, factors) in enumerate(gbp.graph.colForward[variable])
            if factor == factors
                deleteat!(gbp.graph.colForward[variable], k)
                deleteat!(gbp.graph.rowForward[factor], 1)
            end
        end
        if length(gbp.graph.colForward[variable]) == 1 && variable != gbp.graph.root
            push!(gbp.graph.iterateVariable, gbp.graph.rowForward[factor][1])
        end
    end
    gbp.graph.iterateFactor = Int64[]

    if gbp.graph.passFactorVariable + gbp.graph.passVariableFactor == gbp.graph.Nlink
        gbp.graph.iterateVariable = [gbp.graph.root]
        gbp.graph.forward = false
    end
end

########## Backward messages from variable to factor nodes ##########
function backwardVariableFactor(gbp::GraphicalModelTree)
    @inbounds for variable in gbp.graph.iterateVariable
        for factor in gbp.graph.colBackward[variable]
            gbp.graph.passVariableFactor += 1
            gbp.inference.fromVariable[gbp.graph.passVariableFactor] = variable
            gbp.inference.toFactor[gbp.graph.passVariableFactor] = factor

            Mcol = gbp.graph.meanDirect[variable]; Wcol = gbp.graph.weightDirect[variable]
            for k in gbp.graph.incomingToVariable[variable]
                if factor != gbp.inference.fromFactor[k]
                    Mcol += gbp.inference.meanFactorVariable[k] / gbp.inference.varianceFactorVariable[k]
                    Wcol += 1 / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.varianceVariableFactor[gbp.graph.passVariableFactor] = 1 / Wcol
            gbp.inference.meanVariableFactor[gbp.graph.passVariableFactor] = Mcol / Wcol

            push!(gbp.graph.incomingToFactor[factor], gbp.graph.passVariableFactor)
            for (p, variables) in enumerate(gbp.graph.rowBackward[factor])
                if variable == variables
                    deleteat!(gbp.graph.rowBackward[factor], p)
                end
            end
            push!(gbp.graph.iterateFactor, factor)
        end
    end
    gbp.graph.iterateVariable = Int64[]
end

########## Backward messages from factor to variable nodes ##########
function backwardFactorVariable(gbp::GraphicalModelTree)
    @inbounds for factor in gbp.graph.iterateFactor
        for (kk, variable) in enumerate(gbp.graph.rowBackward[factor])
            gbp.graph.passFactorVariable += 1

            gbp.inference.fromFactor[gbp.graph.passFactorVariable] = factor
            gbp.inference.toVariable[gbp.graph.passFactorVariable] = variable

            Mrow = gbp.system.observation[factor]; Vrow = gbp.system.variance[factor]
            for j in gbp.graph.incomingToFactor[factor]
                if gbp.inference.fromVariable[j] != variable
                    coeff = gbp.system.jacobian[factor, gbp.inference.fromVariable[j]]
                    Mrow -= coeff * gbp.inference.meanVariableFactor[j]
                    Vrow += coeff^2 * gbp.inference.varianceVariableFactor[j]
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.passFactorVariable] = Mrow / gbp.system.jacobian[factor, variable]
            gbp.inference.varianceFactorVariable[gbp.graph.passFactorVariable] = Vrow / (gbp.system.jacobian[factor, variable]^2)
            push!(gbp.graph.incomingToVariable[variable], gbp.graph.passFactorVariable)
            for (k, factors) in enumerate(gbp.graph.colBackward[variable])
                if factor == factors
                    deleteat!(gbp.graph.colBackward[variable], k)
                end
            end
            if length(gbp.graph.colBackward[variable]) != 0
                push!(gbp.graph.iterateVariable, variable)
            end
        end
    end
    gbp.graph.iterateFactor = Int64[]

    if gbp.graph.passFactorVariable + gbp.graph.passVariableFactor == 2 * gbp.graph.Nlink
        gbp.graph.backward = false
    end
end