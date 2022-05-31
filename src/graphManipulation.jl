######### Freeze factor node ##########
function freezeFactor!(gbp::ContinuousModel; factor = 0::Int64)
    errorNodeIndex(factor, gbp.graph.Nfactor; name = "factor")

    factorLocal = gbp.graph.dynamic[factor]
    errorFactorLocal(factorLocal; state = "frozen")

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateFactor)
        if gbp.graph.iterateFactor[i] == factorLocal
            whereIs = i
            break
        end
    end
    errorWhereIs(whereIs; name = "factor node", state = "frozen")

    deleteat!(gbp.graph.iterateFactor, whereIs)
end

######### Defreeze factor node ##########
function defreezeFactor!(gbp::ContinuousModel; factor = 0::Int64)
    errorNodeIndex(factor, gbp.graph.Nfactor; name = "factor")

    factorLocal = gbp.graph.dynamic[factor]
    errorFactorLocal(factorLocal; state = "defrozen")

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateFactor)
        if gbp.graph.iterateFactor[i] > factorLocal
            whereIs = i
            break
        end
    end

    if whereIs != 0
        errorInside(whereIs, gbp.graph.iterateFactor, factorLocal; name = "factor node", state = "defrozen")
        insert!(gbp.graph.iterateFactor, whereIs, factorLocal)
    else
        errorLast(gbp.graph.iterateFactor, factorLocal; name = " factor node", state = "defrozen")
        push!(gbp.graph.iterateFactor, factorLocal)
    end
end

######### Freeze variable node ##########
function freezeVariable!(gbp::ContinuousModel; variable = 0::Int64)
    errorNodeIndex(variable, gbp.graph.Nvariable; name = "variable")

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateVariable)
        if gbp.graph.iterateVariable[i] == variable
            whereIs = i
            break
        end
    end
    errorWhereIs(whereIs; name = "variable node", state = "frozen")

    deleteat!(gbp.graph.iterateVariable, whereIs)
end

########## Defreeze variable node ##########
function defreezeVariable!(gbp::ContinuousModel; variable = 0::Int64)
    errorNodeIndex(variable, gbp.graph.Nvariable; name = "variable")

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateVariable)
        if gbp.graph.iterateVariable[i] > variable
            whereIs = i
            break
        end
    end

    if whereIs != 0
        errorInside(whereIs, gbp.graph.iterateVariable, variable; name = "variable node", state = "defrozen")
        insert!(gbp.graph.iterateVariable, whereIs, variable)
    else
        errorLast(gbp.graph.iterateVariable, variable; name = "variable node", state = "defrozen")
        push!(gbp.graph.iterateVariable, variable)
    end
end

######### Freeze Egdge: From variable to factor node ##########
function freezeVariableFactor!(gbp::ContinuousModel; variable = 0::Int64, factor = 0::Int64)
    errorNodeIndex(factor, gbp.graph.Nfactor; name = "factor")
    errorNodeIndex(variable, gbp.graph.Nvariable; name = "variable")
    errorFactorLocal(gbp.graph.dynamic[factor]; state = "frozen")
    errorEdge(gbp.system.jacobian, factor, variable)

    whereIs = 0
    @inbounds for (k, i) in enumerate(gbp.graph.colptr[variable])
        if gbp.inference.fromFactor[i] == factor
            whereIs = k
        end
    end
    errorWhereIs(whereIs; name = "edge", state = "frozen")

    deleteat!(gbp.graph.colptr[variable], whereIs)
    gbp.graph.toFactor[factor, variable] = 0
end

######### Defreeze Egdge: From variable to factor node ##########
function defreezeVariableFactor!(gbp::ContinuousModel; variable = 0::Int64, factor = 0::Int64)
    errorNodeIndex(factor, gbp.graph.Nfactor; name = "factor")
    errorNodeIndex(variable, gbp.graph.Nvariable; name = "variable")
    errorFactorLocal(gbp.graph.dynamic[factor]; state = "frozen")
    errorEdge(gbp.system.jacobian, factor, variable)

    @inbounds for i = 1:gbp.graph.Nlink
        if gbp.inference.fromFactor[i] == factor && gbp.inference.toVariable[i] == variable
            for j in gbp.graph.colptr[variable]
                if j == i
                    error("The edge is already defrozen.")
                end
            end
            push!(gbp.graph.colptr[variable], i)
            break
        end
    end

    @inbounds for i = 1:gbp.graph.Nlink
        if gbp.inference.fromVariable[i] == variable && gbp.inference.toFactor[i] == factor
            gbp.graph.toFactor[factor, variable] = i
            break
        end
    end
    sort!(gbp.graph.colptr[variable])
end

######### Freeze Egdge: From factor to variable node ##########
function freezeFactorVariable!(gbp::ContinuousModel; factor = 0::Int64, variable = 0::Int64)
    errorNodeIndex(factor, gbp.graph.Nfactor; name = "factor")
    errorNodeIndex(variable, gbp.graph.Nvariable; name = "variable")

    factorLocal = gbp.graph.dynamic[factor]
    errorFactorLocal(factorLocal; state = "frozen")
    errorEdge(gbp.system.jacobian, factor, variable)

    whereIs = 0
    @inbounds for (k, i) in enumerate(gbp.graph.rowptr[factorLocal])
        if gbp.inference.fromVariable[i] == variable
            whereIs = k
        end
    end
    errorWhereIs(whereIs; name = "edge", state = "frozen")

    deleteat!(gbp.graph.rowptr[factorLocal], whereIs)
    gbp.graph.toVariable[variable, factor] = 0
end

######### Defreeze Egdge: From factor to variable node ##########
function defreezeFactorVariable!(gbp::ContinuousModel; factor = 0::Int64, variable = 0::Int64)
    errorNodeIndex(factor, gbp.graph.Nfactor; name = "factor")
    errorNodeIndex(variable, gbp.graph.Nvariable; name = "variable")

    factorLocal = gbp.graph.dynamic[factor]
    errorFactorLocal(factorLocal; state = "frozen")
    errorEdge(gbp.system.jacobian, factor, variable)

    @inbounds for i = 1:gbp.graph.Nlink
        if gbp.inference.fromVariable[i] == variable && gbp.inference.toFactor[i] == factor
            for j in gbp.graph.rowptr[factorLocal]
                if j == i
                    error("The edge is already defrozen.")
                end
            end
            push!(gbp.graph.rowptr[factorLocal], i)
            break
        end
    end

    @inbounds for i = 1:gbp.graph.Nlink
        if gbp.inference.fromFactor[i] == factor && gbp.inference.toVariable[i] == variable
            gbp.graph.toVariable[variable, factor] = i
            break
        end
    end
    sort!(gbp.graph.rowptr[factorLocal])
end

######### Hide factor node ##########
function hideFactor!(gbp::ContinuousModel; factor = 0::Int64)
    errorNodeIndex(factor, gbp.graph.Nfactor; name = "factor")

    factorLocal = gbp.graph.dynamic[factor]
    if factorLocal != 0
        whereIs = 0
        @inbounds for i = 1:length(gbp.graph.iterateFactor)
            if gbp.graph.iterateFactor[i] == factorLocal
                whereIs = i
                break
            end
        end
        errorWhereIs(whereIs; name = "factor node", state = "hidden")
        deleteat!(gbp.graph.iterateFactor, whereIs)

        deleteLinksRow = gbp.graph.rowptr[factorLocal]
        variables = gbp.inference.fromVariable[deleteLinksRow]
        deleteLinksCol = Int[]
        @inbounds for i in variables
            links = gbp.graph.colptr[i]
            for (k, j) in enumerate(links)
                if gbp.inference.fromFactor[j] == factor
                    push!(deleteLinksCol, j)
                    deleteat!(gbp.graph.colptr[i], k)
                    deleteat!(gbp.graph.colptrMarginal[i], k)
                end
            end
        end

        @inbounds for i in deleteLinksRow
            gbp.inference.meanVariableFactor[i] = 0.0
            gbp.inference.varianceVariableFactor[i] = 0.0
        end
        @inbounds for i in deleteLinksCol
            gbp.inference.meanFactorVariable[i] = 0.0
            gbp.inference.varianceFactorVariable[i] = 0.0
        end
    else
        idx = gbp.system.jacobianTranspose.colptr[factor]
        variable = gbp.system.jacobianTranspose.rowval[idx]

        gbp.graph.meanDirect[variable] -= gbp.system.observation[factor] * gbp.system.jacobianTranspose[variable, factor] / gbp.system.variance[factor]
        gbp.graph.weightDirect[variable] -= gbp.system.jacobianTranspose[variable, factor]^2 / gbp.system.variance[factor]
        if gbp.graph.weightDirect[variable] == 0.0
            gbp.graph.weightDirect[variable] = 1 /  gbp.graph.virtualVariance
            gbp.graph.meanDirect[variable] =  gbp.graph.virtualMean /  gbp.graph.virtualVariance
        end
    end
    gbp.system.observation[factor] = 0.0
    @inbounds for i in gbp.system.jacobianTranspose.colptr[factor]:(gbp.system.jacobianTranspose.colptr[factor + 1] - 1)
        col = gbp.system.jacobianTranspose.rowval[i]
        gbp.system.jacobianTranspose.nzval[i] = 0.0
        gbp.system.jacobian[factor, col] = 0.0
        gbp.graph.toFactor[factor, col] = 0
    end
end

######### Add factor nodes ##########
function addFactors!(gbp::ContinuousModel; mean = 0.0, variance = 0.0, jacobian = [])
    ### Update system model
    append!(gbp.system.observation, mean)
    append!(gbp.system.variance, variance)
    newFactors = sparse(jacobian)
    newFactorsTranspose = copy(transpose(newFactors))
    gbp.system.jacobian = [gbp.system.jacobian; newFactors]
    gbp.system.jacobianTranspose = copy(transpose(gbp.system.jacobian))

    ### Add singly connected factor nodes and update dynamic vector
    Nfactor = length(mean)
    Nlink = 0; Nindirect = 0
    append!(gbp.graph.dynamic, fill(0, Nfactor))
    @inbounds for i = 1:Nfactor
        NvariableInRow = newFactorsTranspose.colptr[i + 1] - newFactorsTranspose.colptr[i]
        if NvariableInRow == 1
            variable = newFactorsTranspose.rowval[newFactorsTranspose.colptr[i]]
            if gbp.graph.weightDirect[variable] == 1 / gbp.graph.virtualVariance &&  gbp.graph.meanDirect[variable] == gbp.graph.virtualMean /  gbp.graph.virtualVariance
                gbp.graph.meanDirect[variable] = 0.0
                gbp.graph.weightDirect[variable] = 0.0
            end
            gbp.graph.meanDirect[variable] += mean[i] * newFactors[i, variable] / variance[i]
            gbp.graph.weightDirect[variable] += newFactors[i, variable]^2 / variance[i]
        else
            Nlink += NvariableInRow
            Nindirect += 1
            gbp.graph.dynamic[i + gbp.graph.Nfactor] = Nindirect + gbp.graph.Nindirect
        end
    end

    ### Update vectors related with indirect factor nodes
    toFactorLocal = fill(0, Nlink); temp = fill(0.0, Nlink)
    append!(gbp.graph.rowptr, [Int[] for i = 1:Nindirect])
    append!(gbp.graph.coefficient, temp)
    append!(gbp.graph.meanIndirect, temp)
    append!(gbp.graph.varianceIndirect, temp)
    append!(gbp.inference.toFactor, toFactorLocal)
    append!(gbp.inference.fromVariable, toFactorLocal)
    append!(gbp.inference.meanVariableFactor, temp)
    append!(gbp.inference.varianceVariableFactor, temp)
    append!(gbp.graph.alphaOld, temp)
    append!(gbp.graph.alphaNew, fill(1.0, Nlink))

    varianceInitial = fill(0.0, gbp.graph.Nvariable)
    meanInitial = fill(0.0, gbp.graph.Nvariable)
    idxr = gbp.graph.Nindirect + 1
    idxi = gbp.graph.Nlink + 1
    prev = 1
    @inbounds for col = 1:Nfactor
        if (newFactorsTranspose.colptr[col + 1] - newFactorsTranspose.colptr[col]) != 1
            for i = newFactorsTranspose.colptr[col]:(newFactorsTranspose.colptr[col + 1] - 1)
                row = newFactorsTranspose.rowval[i]

                gbp.graph.coefficient[idxi] = newFactorsTranspose[row, col]
                gbp.graph.meanIndirect[idxi] = mean[col]
                gbp.graph.varianceIndirect[idxi] = variance[col]
                gbp.inference.toFactor[idxi] = col + gbp.graph.Nfactor
                gbp.inference.fromVariable[idxi] = row

                toFactorLocal[idxi - gbp.graph.Nlink] = col
                push!(gbp.graph.rowptr[idxr], idxi)

                if meanInitial[row] == 0 && varianceInitial[row] == 0
                    Mcol = gbp.graph.meanDirect[row]; Wcol = gbp.graph.weightDirect[row]
                    for j in gbp.graph.colptrMarginal[row]
                        Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
                        Wcol += 1 / gbp.inference.varianceFactorVariable[j]
                    end
                    varianceInitial[row] = 1 / Wcol
                    meanInitial[row] = Mcol * varianceInitial[row]
                end

                idxi += 1
                rown = newFactorsTranspose.rowval[newFactorsTranspose.colptr[col + 1] - 1]
                if CartesianIndex(rown, col) == CartesianIndex(row, col)
                    idxr += 1
                end
            end
        end
    end

    @inbounds for i = gbp.graph.Nlink + 1:gbp.graph.Nlink + Nlink
        vari = gbp.inference.fromVariable[i]
        gbp.inference.meanVariableFactor[i] = meanInitial[vari]
        gbp.inference.varianceVariableFactor[i] = varianceInitial[vari]
    end

    start = gbp.graph.Nlink + 1
    last = length(gbp.inference.toFactor)
    links = collect(start:last)
    sendToFactor = sparse(toFactorLocal, gbp.inference.fromVariable[start:last], links, Nfactor, gbp.graph.Nvariable)
    gbp.graph.toFactor = [gbp.graph.toFactor; sendToFactor]

    ### Update vectors related with variable nodes
    append!(gbp.inference.toVariable, toFactorLocal)
    append!(gbp.inference.fromFactor, toFactorLocal)
    gbp.graph.toVariable = [gbp.graph.toVariable transpose(sendToFactor)]
    idxi = 1; prev = 1
    @inbounds for col = 1:gbp.graph.Nvariable
        for i = gbp.system.jacobian.colptr[col]:(gbp.system.jacobian.colptr[col + 1] - 1)
            row = gbp.system.jacobian.rowval[i]
            if gbp.system.jacobianTranspose.colptr[row + 1] - gbp.system.jacobianTranspose.colptr[row] != 1
                if prev == col
                    gbp.graph.colptr[col] = Int64[]
                    gbp.graph.colptrMarginal[col] = Int64[]
                end
                prev = col + 1

                if gbp.graph.toFactor[row, col] != 0
                    push!(gbp.graph.colptr[col], idxi)
                    push!(gbp.graph.colptrMarginal[col], idxi)
                end
                gbp.inference.toVariable[idxi] = col
                gbp.inference.fromFactor[idxi] = row

                if gbp.graph.toVariable[col, row] != 0
                    gbp.graph.toVariable[col, row] = idxi
                end
                idxi += 1
            end
        end
    end
    append!(gbp.inference.meanFactorVariable, temp)
    append!(gbp.inference.varianceFactorVariable, temp)
    append!(gbp.graph.iterateFactor, collect(gbp.graph.Nindirect + 1:gbp.graph.Nindirect + Nindirect))

    ### Update factor graf statistics
    gbp.graph.Nfactor += Nfactor
    gbp.graph.Nindirect += Nindirect
    gbp.graph.Nlink += Nlink
end


######### Error node ##########
@inline function errorNodeIndex(node, limit; name = "")
    if node < 0 || node > limit
        error("The $name node does not exist.")
    end
end

######### Error where is ##########
@inline function errorWhereIs(whereIs; name = "", state = "")
    if whereIs == 0
        error("The $name does not exist or it is already $state.")
    end
end

######### Error local factor ##########
@inline function errorFactorLocal(factorLocal; state = "")
    if factorLocal == 0
        error("The singly connected factor node cannot be $state.")
    end
end

######### Error inside ##########
@inline function errorInside(whereIs, iterate, node; name = "", state = "")
    if whereIs != 1 && iterate[whereIs - 1] == node || iterate[1] == node
        error("The $name is already $state.")
    end
end

######### Error last ##########
@inline function errorLast(iterate, node; name = "", state = "")
    if iterate[end] == node
        error("The $name is already $state.")
    end
end

######### Error egde ##########
@inline function errorEdge(jacobian, factor, variable)
    if jacobian[factor, variable] == 0.0
        error("The edge does not exist.")
    end
end
