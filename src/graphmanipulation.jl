######### Freeze factor node ##########
function freezeFactor!(gbp::GraphicalModel; factor = 0::Int64)
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
function defreezeFactor!(gbp::GraphicalModel; factor = 0::Int64)
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
        errorInside(whereIs, gbp.graph.iterateFactor, factorLocal; name = "factor node")
        insert!(gbp.graph.iterateFactor, whereIs, factorLocal)
    else
        errorLast(gbp.graph.iterateFactor, factorLocal; name = " factor node")
        push!(gbp.graph.iterateFactor, factorLocal)
    end
end

######### Freeze variable node ##########
function freezeVariable!(gbp::GraphicalModel; variable = 0::Int64)
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
function defreezeVariable!(gbp::GraphicalModel; variable = 0::Int64)
    errorNodeIndex(variable, gbp.graph.Nvariable; name = "variable")

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateVariable)
        if gbp.graph.iterateVariable[i] > variable
            whereIs = i
            break
        end
    end

    if whereIs != 0
        errorInside(whereIs, gbp.graph.iterateVariable, variable; name = "variable node")
        insert!(gbp.graph.iterateVariable, whereIs, variable)
    else
        errorLast(gbp.graph.iterateVariable, variable; name = "variable node")
        push!(gbp.graph.iterateVariable, variable)
    end
end

######### Freeze Egdge: From variable to factor node ##########
function freezeVariableFactor!(gbp::GraphicalModel; variable = 0::Int64, factor = 0::Int64)
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
end

######### Defreeze Egdge: From variable to factor node ##########
function defreezeVariableFactor!(gbp::GraphicalModel; variable = 0::Int64, factor = 0::Int64)
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
    sort!(gbp.graph.colptr[variable])
end

######### Freeze Egdge: From factor to variable node ##########
function freezeFactorVariable!(gbp::GraphicalModel; factor = 0::Int64, variable = 0::Int64)
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
end

######### Defreeze Egdge: From factor to variable node ##########
function defreezeFactorVariable!(gbp::GraphicalModel; factor = 0::Int64, variable = 0::Int64)
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
    sort!(gbp.graph.rowptr[factorLocal])
end

######### Hide factor node ##########
function hideFactor!(gbp::GraphicalModel; factor = 0::Int64)
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
    end
end

######### Add factor node ##########
function addFactors!(gbp::GraphicalModel; mean = 0.0, variance = 0.0, jacobian = [])
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
    idx = findall(!iszero, newFactorsTranspose)
    idxr = gbp.graph.Nindirect + 1
    idxi = gbp.graph.Nlink + 1
    prev = 1
    @inbounds for i in idx
        if (newFactorsTranspose.colptr[i[2] + 1] - newFactorsTranspose.colptr[i[2]]) != 1
            gbp.graph.coefficient[idxi] = newFactorsTranspose[i]
            gbp.graph.meanIndirect[idxi] = mean[i[2]]
            gbp.graph.varianceIndirect[idxi] = variance[i[2]]
            gbp.inference.toFactor[idxi] = i[2] + gbp.graph.Nfactor
            gbp.inference.fromVariable[idxi] = i[1]

            toFactorLocal[idxi - gbp.graph.Nlink] = i[2]
            push!(gbp.graph.rowptr[idxr], idxi)

            if meanInitial[i[1]] == 0 && varianceInitial[i[1]] == 0
                Mcol = gbp.graph.meanDirect[i[1]]; Wcol = gbp.graph.weightDirect[i[1]]
                for j in gbp.graph.colptrMarginal[i[1]]
                    Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
                    Wcol += 1 / gbp.inference.varianceFactorVariable[j]
                end
                varianceInitial[i[1]] = 1 / Wcol
                meanInitial[i[1]] = Mcol * varianceInitial[i[1]]
            end

            idxi += 1
            if idx[newFactorsTranspose.colptr[i[2] + 1] - 1] == i
                idxr += 1
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
    dropzeros!(gbp.system.jacobian)
    idx = findall(!iszero, gbp.system.jacobian); idxi = 1; prev = 1
    @inbounds for i in idx
        if gbp.system.jacobianTranspose.colptr[i[1] + 1] - gbp.system.jacobianTranspose.colptr[i[1]] != 1
            if prev == i[2]
                gbp.graph.colptr[i[2]] = Int64[]
                gbp.graph.colptrMarginal[i[2]] = Int64[]
            end
            prev = i[2] + 1

            push!(gbp.graph.colptr[i[2]], idxi)
            push!(gbp.graph.colptrMarginal[i[2]], idxi)
            gbp.inference.toVariable[idxi] = i[2]
            gbp.inference.fromFactor[idxi] = i[1]
            idxi += 1
        end
    end

    links = collect(1:idxi - 1)
    gbp.graph.toVariable = sparse(gbp.inference.toVariable, gbp.inference.fromFactor, links, gbp.graph.Nvariable, gbp.graph.Nfactor + Nfactor)

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
@inline function errorInside(whereIs, iterate, node; name = "")
    if whereIs != 1 && iterate[whereIs - 1] == node || iterate[1] == node
        error("The $name is already defrozen.")
    end
end

######### Error last ##########
@inline function errorLast(iterate, node; name = "")
    if iterate[end] == node
        error("The $name is already defrozen.")
    end
end

######### Error egde ##########
@inline function errorEdge(jacobian, factor, variable)
    if jacobian[factor, variable] == 0.0
        error("The edge does not exist.")
    end
end