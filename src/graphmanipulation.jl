######### Freeze factor node ##########
function freezeFactor!(gbp::GraphicalModel; factor = 0::Int64)
    if factor == 0
        error("The keyword factor is missing.")
    elseif factor > gbp.graph.Nfactor
        error("The factor node does not exist.")
    end

    factorLocal = gbp.graph.dynamic[factor]
    if factorLocal == 0
        error("The singly connected factor node cannot be frozen.")
    end

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateFactor)
        if gbp.graph.iterateFactor[i] == factorLocal
            whereIs = i
            break
        end
    end
    if whereIs == 0
        error("The factor node does not exist or it is already frozen.")
    end

    deleteat!(gbp.graph.iterateFactor, whereIs)
end

######### Defreeze factor node ##########
function defreezeFactor!(gbp::GraphicalModel; factor = 0::Int64)
    if factor == 0
        error("The keyword factor is missing.")
    elseif factor > gbp.graph.Nfactor
        error("The factor node does not exist.")
    end

    factorLocal = gbp.graph.dynamic[factor]
    if factorLocal == 0
        error("The singly connected factor node cannot be defrozen.")
    end

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateFactor)
        if gbp.graph.iterateFactor[i] > factorLocal
            whereIs = i
            break
        end
    end

    if whereIs != 0
        if whereIs != 1 && gbp.graph.iterateFactor[whereIs - 1] == factorLocal || gbp.graph.iterateFactor[1] == factorLocal
            error("The factor node is already defrozen.")
        end
        insert!(gbp.graph.iterateFactor, whereIs, factorLocal)
    else
        if gbp.graph.iterateFactor[end] == factorLocal
            error("The factor node is already defrozen.")
        end
        push!(gbp.graph.iterateFactor, factorLocal)
    end
end

######### Freeze variable node ##########
function freezeVariable!(gbp::GraphicalModel; variable = 0::Int64)
    if variable == 0
        error("The keyword variable is missing.")
    elseif variable > gbp.graph.Nvariable
        error("The variable node does not exist.")
    end

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateVariable)
        if gbp.graph.iterateVariable[i] == variable
            whereIs = i
            break
        end
    end
    if whereIs == 0
        error("The variable node does not exist or it is already frozen.")
    end

    deleteat!(gbp.graph.iterateVariable, whereIs)
end

########## Defreeze variable node ##########
function defreezeVariable!(gbp::GraphicalModel; variable = 0::Int64)
    if variable == 0
        error("The keyword variable is missing.")
    elseif variable > gbp.graph.Nvariable
        error("The variable node does not exist.")
    end

    whereIs = 0
    @inbounds for i = 1:length(gbp.graph.iterateVariable)
        if gbp.graph.iterateVariable[i] > variable
            whereIs = i
            break
        end
    end

    if whereIs != 0
        if whereIs != 1 && gbp.graph.iterateVariable[whereIs - 1] == variable || gbp.graph.iterateVariable[1] == variable
            error("The variable node is already defrozen.")
        end
        insert!(gbp.graph.iterateVariable, whereIs, variable)
    else
        if gbp.graph.iterateVariable[end] == variable
            error("The variable node is already defrozen.")
        end
        push!(gbp.graph.iterateVariable, variable)
    end
end

######### Freeze Egdge: From variable to factor node ##########
function freezeVariableFactor!(gbp::GraphicalModel; variable = 0::Int64, factor = 0::Int64)
    if variable == 0 || factor == 0
        error("The keyword variable or factor is missing.")
    elseif variable > gbp.graph.Nvariable || factor > gbp.graph.Nfactor
        error("The variable or factor node does not exist.")
    end

    value = gbp.system.jacobian[factor, variable]
    if value == 0
        error("The edge does not exist.")
    end

    whereIs = 0
    @inbounds for (k, i) in enumerate(gbp.graph.colptr[variable])
        if gbp.inference.fromFactor[i] == factor
            whereIs = k
        end
    end
    if whereIs == 0
        error("The edge node does not exist or it is already frozen.")
    end
    deleteat!(gbp.graph.colptr[variable], whereIs)
end

######### Freeze Egdge: From factor to variable node ##########
function freezeFactorVariable!(gbp::GraphicalModel; factor = 0::Int64, variable = 0::Int64)
    if variable == 0 || factor == 0
        error("The keyword variable or factor is missing.")
    elseif variable > gbp.graph.Nvariable || factor > gbp.graph.Nfactor
        error("The variable or factor node does not exist.")
    end

    value = gbp.system.jacobian[factor, variable]
    if value == 0
        error("The edge does not exist.")
    end

    factorLocal = gbp.graph.dynamic[factor]
    whereIs = 0
    @inbounds for (k, i) in enumerate(gbp.graph.rowptr[factorLocal])
        if gbp.inference.fromVariable[i] == variable
            whereIs = k
        end
    end
    if whereIs == 0
        error("The edge node does not exist or it is already frozen.")
    end
    deleteat!(gbp.graph.rowptr[factorLocal], whereIs)
end

