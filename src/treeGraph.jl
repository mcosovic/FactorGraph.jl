########## Continuous model: Forward messages from variable to factor nodes ##########
function forwardVariableFactor(gbp::ContinuousTreeModel)
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
            end
        end
        if length(gbp.graph.rowForward[factor]) == 1
            push!(gbp.graph.iterateFactor, factor)
        end
    end
    gbp.graph.iterateVariable = Int64[]
end

########## Continuous model: Forward messages from factor to variable nodes ##########
function forwardFactorVariable(gbp::ContinuousTreeModel)
    @inbounds for factor in gbp.graph.iterateFactor
        variable = gbp.graph.rowForward[factor][1]

        gbp.graph.passFactorVariable += 1
        gbp.inference.fromFactor[gbp.graph.passFactorVariable] = factor
        gbp.inference.toVariable[gbp.graph.passFactorVariable] = variable

        Mrow = gbp.system.observation[factor]; Vrow = gbp.system.variance[factor]
        for j in gbp.graph.incomingToFactor[factor]
            coeff = gbp.system.coefficient[gbp.inference.toFactor[j], gbp.inference.fromVariable[j]]
            Mrow -= coeff * gbp.inference.meanVariableFactor[j]
            Vrow += coeff^2 * gbp.inference.varianceVariableFactor[j]
        end
        gbp.inference.meanFactorVariable[gbp.graph.passFactorVariable] = Mrow / gbp.system.coefficient[factor, variable]
        gbp.inference.varianceFactorVariable[gbp.graph.passFactorVariable] = Vrow / (gbp.system.coefficient[factor, variable]^2)

        push!(gbp.graph.incomingToVariable[variable], gbp.graph.passFactorVariable)
        for (k, factors) in enumerate(gbp.graph.colForward[variable])
            if factor == factors
                deleteat!(gbp.graph.colForward[variable], k)
            end
        end
        if length(gbp.graph.colForward[variable]) == 1 && variable != gbp.graph.root
            push!(gbp.graph.iterateVariable, variable)
        end
    end
    gbp.graph.iterateFactor = Int64[]

    if gbp.graph.passFactorVariable + gbp.graph.passVariableFactor == gbp.graph.Nlink
        gbp.graph.iterateVariable = [gbp.graph.root]
        gbp.graph.forward = false
    end
end

########## Continuous model: Backward messages from variable to factor nodes ##########
function backwardVariableFactor(gbp::ContinuousTreeModel)
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
            for (k, variables) in enumerate(gbp.graph.rowBackward[factor])
                if variable == variables
                    deleteat!(gbp.graph.rowBackward[factor], k)
                end
            end
            push!(gbp.graph.iterateFactor, factor)
        end
    end
    gbp.graph.iterateVariable = Int64[]
end

########## Continuous model: Backward messages from factor to variable nodes ##########
function backwardFactorVariable(gbp::ContinuousTreeModel)
    @inbounds for factor in gbp.graph.iterateFactor
        for variable in gbp.graph.rowBackward[factor]
            gbp.graph.passFactorVariable += 1

            gbp.inference.fromFactor[gbp.graph.passFactorVariable] = factor
            gbp.inference.toVariable[gbp.graph.passFactorVariable] = variable

            Mrow = gbp.system.observation[factor]; Vrow = gbp.system.variance[factor]
            for j in gbp.graph.incomingToFactor[factor]
                if gbp.inference.fromVariable[j] != variable
                    coeff = gbp.system.coefficient[factor, gbp.inference.fromVariable[j]]
                    Mrow -= coeff * gbp.inference.meanVariableFactor[j]
                    Vrow += coeff^2 * gbp.inference.varianceVariableFactor[j]
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.passFactorVariable] = Mrow / gbp.system.coefficient[factor, variable]
            gbp.inference.varianceFactorVariable[gbp.graph.passFactorVariable] = Vrow / (gbp.system.coefficient[factor, variable]^2)
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

########## Discrete model: Forward messages from variable to factor nodes ##########
function forwardVariableFactor(bp::DiscreteTreeModel)
    @inbounds for variable in bp.graph.iterateVariable
        factor = bp.graph.colForward[variable][1]

        bp.graph.passVariableFactor += 1
        bp.inference.fromVariable[bp.graph.passVariableFactor] = variable
        bp.inference.toFactor[bp.graph.passVariableFactor] = factor

        bp.inference.messageVariableFactor[bp.graph.passVariableFactor] = copy(bp.graph.messageDirect[variable])
        for j in bp.graph.incomingToVariable[variable]
            bp.inference.messageVariableFactor[bp.graph.passVariableFactor] .*= bp.inference.messageFactorVariable[j]
        end

        for (k, variables) in enumerate(bp.system.probability[factor])
            if variables == variable
                bp.graph.incomingToFactor[factor][k] = bp.graph.passVariableFactor
            end
        end

        for (k, variables) in enumerate(bp.graph.rowForward[factor])
            if variable == variables
                deleteat!(bp.graph.rowForward[factor], k)
            end
        end

        if length(bp.graph.rowForward[factor]) == 1
            push!(bp.graph.iterateFactor, factor)
        end
    end
    bp.graph.iterateVariable = Int64[]
end

########## Discrete model: Forward messages from factor to variable nodes ##########
function forwardFactorVariable(bp::DiscreteTreeModel)
    @inbounds for factor in bp.graph.iterateFactor
        variable = bp.graph.rowForward[factor][1]

        bp.graph.passFactorVariable += 1
        bp.inference.fromFactor[bp.graph.passFactorVariable] = factor
        bp.inference.toVariable[bp.graph.passFactorVariable] = variable

        table = deepcopy(bp.system.table[factor])
        bp.inference.messageFactorVariable[bp.graph.passFactorVariable] = fill(0.0, bp.graph.state[variable])
        position = 0
        for idx in CartesianIndices(bp.system.table[factor])
            for (k, variables) in enumerate(bp.system.probability[factor])
                if variables != variable
                    whereIs = bp.graph.incomingToFactor[factor][k]
                    table[idx] *= bp.inference.messageVariableFactor[whereIs][idx[k]]
                else
                    position = k
                end
            end
            bp.inference.messageFactorVariable[bp.graph.passFactorVariable][idx[position]] += table[idx]
        end

        push!(bp.graph.incomingToVariable[variable], bp.graph.passFactorVariable)
        for (k, factors) in enumerate(bp.graph.colForward[variable])
            if factor == factors
                deleteat!(bp.graph.colForward[variable], k)
            end
        end
        if length(bp.graph.colForward[variable]) == 1 && variable != bp.graph.root
            push!(bp.graph.iterateVariable, variable)
        end
    end
    bp.graph.iterateFactor = Int64[]

    if bp.graph.passFactorVariable + bp.graph.passVariableFactor == bp.graph.Nlink
        bp.graph.iterateVariable = [bp.graph.root]
        bp.graph.forward = false
    end
end

########## Discrete model: Backward messages from variable to factor nodes ##########
function backwardVariableFactor(bp::DiscreteTreeModel)
    @inbounds for variable in bp.graph.iterateVariable
        for factor in bp.graph.colBackward[variable]
            bp.graph.passVariableFactor += 1
            bp.inference.fromVariable[bp.graph.passVariableFactor] = variable
            bp.inference.toFactor[bp.graph.passVariableFactor] = factor

            bp.inference.messageVariableFactor[bp.graph.passVariableFactor] = copy(bp.graph.messageDirect[variable])
            for k in bp.graph.incomingToVariable[variable]
                if factor != bp.inference.fromFactor[k]
                    bp.inference.messageVariableFactor[bp.graph.passVariableFactor] .*= bp.inference.messageFactorVariable[k]
                end
            end

            for (k, variables) in enumerate(bp.system.probability[factor])
                if variables == variable
                    bp.graph.incomingToFactor[factor][k] = bp.graph.passVariableFactor
                end
            end

            for (k, variables) in enumerate(bp.graph.rowBackward[factor])
                if variable == variables
                    deleteat!(bp.graph.rowBackward[factor], k)
                end
            end
            push!(bp.graph.iterateFactor, factor)
        end
    end
    bp.graph.iterateVariable = Int64[]
end

########## Discrete model: Backward messages from factor to variable nodes ##########
function backwardFactorVariable(bp::DiscreteTreeModel)
    @inbounds for factor in bp.graph.iterateFactor
        for variable in bp.graph.rowBackward[factor]
            bp.graph.passFactorVariable += 1

            bp.inference.fromFactor[bp.graph.passFactorVariable] = factor
            bp.inference.toVariable[bp.graph.passFactorVariable] = variable

            table = deepcopy(bp.system.table[factor])
            bp.inference.messageFactorVariable[bp.graph.passFactorVariable] = fill(0.0, bp.graph.state[variable])
            position = 0
            for idx in CartesianIndices(bp.system.table[factor])
                for (k, variables) in enumerate(bp.system.probability[factor])
                    if variables != variable
                        whereIs = bp.graph.incomingToFactor[factor][k]
                        table[idx] *= bp.inference.messageVariableFactor[whereIs][idx[k]]
                    else
                        position = k
                    end
                end
               bp.inference.messageFactorVariable[bp.graph.passFactorVariable][idx[position]] += table[idx]
            end

            push!(bp.graph.incomingToVariable[variable], bp.graph.passFactorVariable)
            for (k, factors) in enumerate(bp.graph.colBackward[variable])
                if factor == factors
                    deleteat!(bp.graph.colBackward[variable], k)
                end
            end
            if length(bp.graph.colBackward[variable]) != 0
                push!(bp.graph.iterateVariable, variable)
            end
        end
    end
    bp.graph.iterateFactor = Int64[]

    if bp.graph.passFactorVariable + bp.graph.passVariableFactor == 2 * bp.graph.Nlink
        bp.graph.backward = false
    end
end

########## Continuous model: Check that the factor graph has a tree structure ##########
function isTree(gbp::Union{ContinuousModel, ContinuousTreeModel})
    ## Pass through the rows
    Nlink = copy(gbp.graph.Nlink)
    rowptr = [Int[] for i = 1:gbp.graph.Nfactor]
    iterateFactor = Int64[]
    @inbounds for i = 1:gbp.graph.Nfactor
        if gbp.system.coefficientTranspose.colptr[i + 1] - gbp.system.coefficientTranspose.colptr[i] != 1
            for j = gbp.system.coefficientTranspose.colptr[i]:(gbp.system.coefficientTranspose.colptr[i + 1] - 1)
                row = gbp.system.coefficientTranspose.rowval[j]
                push!(rowptr[i], row)
            end
        end
    end

    ## Pass through the columns
    iterateVariable = fill(0, gbp.graph.Nvariable)
    colptr = [Int[] for col = 1:gbp.graph.Nvariable]
    counter = 0
    @inbounds for col = 1:gbp.graph.Nvariable
        for i = gbp.system.coefficient.colptr[col]:(gbp.system.coefficient.colptr[col + 1] - 1)
            row = gbp.system.coefficient.rowval[i]
            if gbp.system.coefficientTranspose.colptr[row + 1] - gbp.system.coefficientTranspose.colptr[row] != 1
                push!(colptr[col], row)
            end
        end
        if length(colptr[col]) == 1
            counter += 1
            iterateVariable[counter] = col
        end
    end
    resize!(iterateVariable, counter)

    ## Pilling the factor graph
    hasSingleVariable = true; hasSingleFactor = true;
    while hasSingleFactor || hasSingleVariable
        iterateFactor = Int64[]
        for variable in iterateVariable
            factor = colptr[variable][1]
            for (k, variables) in enumerate(rowptr[factor])
                if variable == variables
                    deleteat!(rowptr[factor], k)
                    Nlink -= 1
                end
            end
            colptr[variable] = []
            if length(rowptr[factor]) == 1
                push!(iterateFactor, factor)
            end
        end
        if Nlink == 0 break end
        if isempty(iterateFactor)
            hasSingleFactor = false
        end

        iterateVariable = Int64[]
        for factor in iterateFactor
            variable = rowptr[factor][1]
            for (k, factors) in enumerate(colptr[variable])
                if factor == factors
                    deleteat!(colptr[variable], k)
                    Nlink -= 1
                end
            end
            rowptr[factor] = []
            if length(colptr[variable]) == 1
                push!(iterateVariable, variable)
            end
        end
        if Nlink == 0 break end
        if isempty(iterateVariable)
            hasSingleVariable = false
        end
    end

    isTree = true
    @inbounds for i in colptr
        if length(i) != 0
            isTree = false
            break
        end
    end

    return isTree
end

########## Discrete model: Check that the factor graph has a tree structure ##########
function isTree(gbp::DiscreteTreeModel)
    ## Pass through the rows
    Nlink = copy(gbp.graph.Nlink)
    rowptr = [Int[] for i = 1:gbp.graph.Nfactor]
    iterateFactor = Int64[]
    @inbounds for i = 1:gbp.graph.Nfactor
        if gbp.system.incidenceTranspose.colptr[i + 1] - gbp.system.incidenceTranspose.colptr[i] != 1
            for j = gbp.system.incidenceTranspose.colptr[i]:(gbp.system.incidenceTranspose.colptr[i + 1] - 1)
                row = gbp.system.incidenceTranspose.rowval[j]
                push!(rowptr[i], row)
            end
        end
    end

    ## Pass through the columns
    iterateVariable = fill(0, gbp.graph.Nvariable)
    colptr = [Int[] for col = 1:gbp.graph.Nvariable]
    counter = 0
    @inbounds for col = 1:gbp.graph.Nvariable
        for i = gbp.system.incidence.colptr[col]:(gbp.system.incidence.colptr[col + 1] - 1)
            row = gbp.system.incidence.rowval[i]
            if gbp.system.incidenceTranspose.colptr[row + 1] - gbp.system.incidenceTranspose.colptr[row] != 1
                push!(colptr[col], row)
            end
        end
        if length(colptr[col]) == 1
            counter += 1
            iterateVariable[counter] = col
        end
    end
    resize!(iterateVariable, counter)

    ## Pilling the factor graph
    hasSingleVariable = true; hasSingleFactor = true;
    while hasSingleFactor || hasSingleVariable
        iterateFactor = Int64[]
        for variable in iterateVariable
            factor = colptr[variable][1]
            for (k, variables) in enumerate(rowptr[factor])
                if variable == variables
                    deleteat!(rowptr[factor], k)
                    Nlink -= 1
                end
            end
            colptr[variable] = []
            if length(rowptr[factor]) == 1
                push!(iterateFactor, factor)
            end
        end
        if Nlink == 0 break end
        if isempty(iterateFactor)
            hasSingleFactor = false
        end

        iterateVariable = Int64[]
        for factor in iterateFactor
            variable = rowptr[factor][1]
            for (k, factors) in enumerate(colptr[variable])
                if factor == factors
                    deleteat!(colptr[variable], k)
                    Nlink -= 1
                end
            end
            rowptr[factor] = []
            if length(colptr[variable]) == 1
                push!(iterateVariable, variable)
            end
        end
        if Nlink == 0 break end
        if isempty(iterateVariable)
            hasSingleVariable = false
        end
    end

    isTree = true
    @inbounds for i in colptr
        if length(i) != 0
            isTree = false
            break
        end
    end

    return isTree
end