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
function addFactor!(gbp::GraphicalModel; mean = 0.0:Float64, variance = 0.0::Float64, variable = [])
    Nvariables = size(variable, 1)
    for i = 1:Nvariables
        errorNodeIndex(variable[i, 1], gbp.graph.Nvariable; name = "variable")
    end
    gbp.graph.Nfactor += 1

    ## Update system model
    push!(gbp.system.observation, mean)
    push!(gbp.system.variance, variance)
    newRow = sparse(ones(Nvariables), variable[:, 1], variable[:, 2], 1, gbp.graph.Nvariable)
    gbp.system.jacobian = [gbp.system.jacobian; newRow]
    gbp.system.jacobianTranspose = copy(transpose(gbp.system.jacobian))

    ### Add singly connected factor node
    if Nvariables == 1
        vari = trunc(Int, variable[1, 1]); coef = variable[1, 2]
        if gbp.graph.weightDirect[vari] == 1 / gbp.graph.virtualVariance &&  gbp.graph.meanDirect[vari] == gbp.graph.virtualMean /  gbp.graph.virtualVariance
            gbp.graph.meanDirect[vari] = mean * coef / variance
            gbp.graph.weightDirect[vari] = coef^2 / variance
        else
            gbp.graph.meanDirect[vari] += mean * coef / variance
            gbp.graph.weightDirect[vari] += coef^2 / variance
        end
        push!(gbp.graph.dynamic, 0)
    end

    ### Add indirect factor node
    if Nvariables != 1
        # Compute target marginals for initial messages and set new links from variable to factor node
        idx = sortperm(variable[:, 1])
        vari = trunc.(Int, variable[idx, 1]); coeff = variable[idx, 2]

        append!(gbp.graph.coefficient, coeff)
        append!(gbp.inference.fromVariable, vari)
        append!(gbp.inference.toFactor, fill(gbp.graph.Nfactor, Nvariables))
        append!(gbp.graph.meanIndirect, fill(mean, Nvariables))
        append!(gbp.graph.varianceIndirect, fill(variance, Nvariables))
        append!(gbp.graph.alphaNew, fill(1.0, Nvariables))
        append!(gbp.graph.alphaOld, fill(0.0, Nvariables))

        @inbounds for (k, i) in enumerate(vari)
            Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]
            for j in gbp.graph.colptrMarginal[i]
                Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
                Wcol += 1 / gbp.inference.varianceFactorVariable[j]
            end
            varianceInitial = 1 / Wcol
            meanInitial = Mcol * gbp.inference.variance[i]

            push!(gbp.inference.meanVariableFactor, meanInitial)
            push!(gbp.inference.varianceVariableFactor, varianceInitial)
        end
        push!(gbp.graph.rowptr, collect(gbp.graph.Nlink + 1:(gbp.graph.Nlink + Nvariables)))

        # Set new links from factor to variable
        increment = 0; last = 0
        @inbounds for (k, i) in enumerate(gbp.inference.toVariable)
            if !isempty(vari) && (i - 1) == vari[1]
                insert!(gbp.inference.toVariable, k, vari[1])
                insert!(gbp.inference.fromFactor, k, gbp.graph.Nfactor)
                insert!(gbp.inference.meanFactorVariable, k, 0.0)
                insert!(gbp.inference.varianceFactorVariable, k, gbp.graph.virtualVariance)
                push!(gbp.graph.colptr[vari[1]], k)
                deleteat!(vari, 1)
                increment += 1
            end

            if last != i  gbp.graph.colptr[i] .+= increment end
            last = i
        end
        if !isempty(vari)
            push!(gbp.inference.toVariable, vari[1])
            push!(gbp.inference.fromFactor, gbp.graph.Nfactor)
            push!(gbp.inference.meanFactorVariable, 0.0)
            push!(gbp.inference.varianceFactorVariable, gbp.graph.virtualVariance)
            push!(gbp.graph.colptr[vari[1]], length(gbp.inference.toVariable))
        end
        gbp.graph.Nlink += Nvariables
        gbp.graph.colptrMarginal = deepcopy(gbp.graph.colptr)

        # Add other parameters
        gbp.graph.Nindirect += 1

        push!(gbp.graph.iterateFactor, gbp.graph.Nindirect)
        push!(gbp.graph.dynamic, gbp.graph.Nindirect)

        links = collect(1:gbp.graph.Nlink)
        vari = trunc.(Int, variable[idx, 1])
        gbp.graph.toVariable = sparse(gbp.inference.toVariable, gbp.inference.fromFactor, links, gbp.graph.Nvariable, gbp.graph.Nfactor)
        gbp.graph.toFactor = [gbp.graph.toFactor; sparse(fill(1, length(gbp.graph.rowptr[end])), vari, gbp.graph.rowptr[end], 1, gbp.graph.Nvariable)]
    end
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