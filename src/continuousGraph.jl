mutable struct ContinuousSystem
    coefficient::SparseMatrixCSC{Float64,Int64}
    coefficientTranspose::SparseMatrixCSC{Float64,Int64}
    observation::Array{Float64,1}
    variance::Array{Float64,1}
end

mutable struct ContinuousGraph
    Nvariable::Int64
    Nfactor::Int64
    Nindirect::Int64
    Nlink::Int64
    virtualMean::Float64
    virtualVariance::Float64
    meanDirect::Array{Float64,1}
    weightDirect::Array{Float64,1}
    meanIndirect::Array{Float64,1}
    varianceIndirect::Array{Float64,1}
    coefficient::Array{Float64,1}
    toVariable::SparseMatrixCSC{Int64, Int64}
    toFactor::SparseMatrixCSC{Int64, Int64}
    rowptr::Vector{Vector{Int64}}
    colptr::Vector{Vector{Int64}}
    colptrMarginal::Vector{Vector{Int64}}
    alphaNew::Array{Float64,1}
    alphaOld::Array{Float64,1}
    iterateFactor::Array{Int64,1}
    iterateVariable::Array{Int64,1}
    iterateMarginal::Array{Int64,1}
    dynamic::Array{Int64,1}
end

mutable struct ContinuousTreeGraph
    Nvariable::Int64
    Nfactor::Int64
    Nlink::Int64
    root::Int64
    virtualMean::Float64
    virtualVariance::Float64
    meanDirect::Array{Float64,1}
    weightDirect::Array{Float64,1}
    rowForward::Vector{Vector{Int64}}
    rowBackward::Vector{Vector{Int64}}
    colForward::Vector{Vector{Int64}}
    colBackward::Vector{Vector{Int64}}
    incomingToFactor::Vector{Vector{Int64}}
    incomingToVariable::Vector{Vector{Int64}}
    iterateFactor::Array{Int64,1}
    iterateVariable::Array{Int64,1}
    passFactorVariable::Int64
    passVariableFactor::Int64
    forward::Bool
    backward::Bool
end

mutable struct ContinuousInference
    fromFactor::Array{Int64,1}
    toVariable::Array{Int64,1}
    meanFactorVariable::Array{Float64,1}
    varianceFactorVariable::Array{Float64,1}
    fromVariable::Array{Int64,1}
    toFactor::Array{Int64,1}
    meanVariableFactor::Array{Float64,1}
    varianceVariableFactor::Array{Float64,1}
    mean::Array{Float64,1}
    variance::Array{Float64,1}
end

struct ContinuousModel
    graph::ContinuousGraph
    inference::ContinuousInference
    system::ContinuousSystem
end

struct ContinuousTreeModel
    graph::ContinuousTreeGraph
    inference::ContinuousInference
    system::ContinuousSystem
end

########## Form a graph with continuous variables and initialize messages ##########
function continuousModel(
    args...;
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    mean::Float64 = 0.0,
    variance::Float64 = 1e10)

    checkKeywords(prob, alpha, variance)
    system = readContinuousArguments(args)
    graph, inference = makeContinuousGraph(system, mean, variance, prob, alpha)

    return ContinuousModel(graph, inference, system)
end

########## Form a tree factor graph and initialize messages ##########
function continuousTreeModel(
    args...;
    mean::Float64 = 0.0,
    variance::Float64 = 1e10,
    root::Int64 = 1)

    system = readContinuousArguments(args)
    graph, inference = makeContinuousTreeGraph(system, mean, variance, root)

    return ContinuousTreeModel(graph, inference, system)
end


########## Read system model ##########
function readContinuousArguments(args)
    if typeof(args[1]) == Array{Float64, 2}
        coefficient = sparse(args[1])
    else
        coefficient = args[1]
    end

    coefficientTranspose = copy(transpose(coefficient))
    observation = args[2]
    variance = args[3]

    return ContinuousSystem(coefficient, coefficientTranspose, observation, variance)
end

########## Produce the graphical model ##########
function makeContinuousGraph(system, meanVirtual, varianceVirtual, dampProbability, dampAlpha)
    ### Number of factor and variable nodes
    Nfactor, Nvariable = size(system.coefficient)

    ### Find graph numbers, set the direct mean and variance and internal factor numeration
    Ndirect = 0; Nlink = 0; Nindirect = 0
    meanDirect = fill(meanVirtual / varianceVirtual, Nvariable)
    weightDirect = fill(1 / varianceVirtual, Nvariable)
    dynamic = fill(0, Nfactor)

    @inbounds for i = 1:Nfactor
        NvariableInRow = system.coefficientTranspose.colptr[i + 1] - system.coefficientTranspose.colptr[i]
        if NvariableInRow == 1
            Ndirect += NvariableInRow
            variable = system.coefficientTranspose.rowval[system.coefficientTranspose.colptr[i]]
            meanDirect[variable] = 0.0
            weightDirect[variable] = 0.0
        else
            Nlink += NvariableInRow
            Nindirect += 1
            dynamic[i] = Nindirect
        end
    end

    ### Pass through the columns
    colptr = [Int[] for i = 1:Nvariable]
    colptrMarginal = [Int[] for i = 1:Nvariable]
    toVariable = fill(0, Nlink)
    fromFactor = similar(toVariable)
    idxi = 1
    @inbounds for col = 1:Nvariable
        for i = system.coefficient.colptr[col]:(system.coefficient.colptr[col + 1] - 1)
            row = system.coefficient.rowval[i]
            NvariableInRow = system.coefficientTranspose.colptr[row + 1] - system.coefficientTranspose.colptr[row]
            if NvariableInRow == 1
                meanDirect[col] += system.observation[row] * system.coefficient[row, col] / system.variance[row]
                weightDirect[col] += system.coefficient[row, col]^2 / system.variance[row]
            else
                push!(colptr[col], idxi)
                push!(colptrMarginal[col], idxi)
                toVariable[idxi] = col
                fromFactor[idxi] = row
                idxi += 1
            end
        end
    end

    ### Pass through the rows and send messages from singly-connected factor nodes to all indirect links
    coefficient = fill(0.0, Nlink)
    meanIndirect = similar(coefficient)
    varianceIndirect = similar(coefficient)
    meanFactorVariable = similar(coefficient)
    varianceFactorVariable = similar(coefficient)
    meanVariableFactor = similar(coefficient)
    varianceVariableFactor = similar(coefficient)

    fromVariable = fill(0, Nlink)
    toFactor = similar(fromVariable)
    rowptr = [Int[] for i = 1:Nindirect]
    idxi = 1
    @inbounds for (col, val) in enumerate(dynamic)
        if val != 0
            for i = system.coefficientTranspose.colptr[col]:(system.coefficientTranspose.colptr[col + 1] - 1)
                row = system.coefficientTranspose.rowval[i]

                coefficient[idxi] = system.coefficientTranspose[row, col]
                meanIndirect[idxi] = system.observation[col]
                varianceIndirect[idxi] = system.variance[col]

                toFactor[idxi] = col
                fromVariable[idxi] = row
                push!(rowptr[dynamic[col]], idxi)

                varianceVariableFactor[idxi] = 1 / weightDirect[row]
                meanVariableFactor[idxi] = meanDirect[row] * varianceVariableFactor[idxi]
                idxi += 1
            end
        end
    end

    ### Message send indices
    links = collect(1:idxi - 1)
    sendToFactor = sparse(toFactor, fromVariable, links, Nfactor, Nvariable)
    sendToVariable = sparse(toVariable, fromFactor, links, Nvariable, Nfactor)

    ### Set damping parameters
    alphaNew = fill(1.0, Nlink)
    alphaOld = fill(0.0, Nlink)
    bernoulliSample = randsubseq(collect(1:Nlink), dampProbability)
    @inbounds for i in bernoulliSample
        alphaNew[i] = 1.0 - dampAlpha
        alphaOld[i] = dampAlpha
    end

    ### Initialize marginal mean and variance vectors
    mean = fill(0.0, Nvariable)
    variance = fill(0.0, Nvariable)

    ### Iteration counters
    iterateFactor = collect(1:Nindirect)
    iterateVariable = collect(1:Nvariable)
    iterateMarginal = copy(iterateVariable)

    return ContinuousGraph(Nvariable, Nfactor, Nindirect, Nlink,
            meanVirtual, varianceVirtual, meanDirect, weightDirect, meanIndirect, varianceIndirect, coefficient,
            sendToVariable, sendToFactor, rowptr, colptr, colptrMarginal, alphaNew, alphaOld,
            iterateFactor, iterateVariable, iterateMarginal, dynamic),
           ContinuousInference(fromFactor, toVariable, meanFactorVariable, varianceFactorVariable,
            fromVariable, toFactor, meanVariableFactor, varianceVariableFactor, mean, variance)
end

########## Produce the tree graphical model ##########
function makeContinuousTreeGraph(system, virtualMean, virtualVariance, root)
    ### Number of factor and variable nodes
    Nfactor, Nvariable = size(system.coefficient)

    ### Find graph numbers, set the direct mean and variance, set factor numeration and pass through the rows
    Nlink = 0; Nindirect = 0
    meanDirect = fill(virtualMean / virtualVariance, Nvariable)
    weightDirect = fill(1 / virtualVariance, Nvariable)

    rowForward = [Int[] for i = 1:Nfactor]; rowBackward = [Int[] for i = 1:Nfactor]
    incomingToVariable = [Int[] for i = 1:Nvariable]
    @inbounds for i = 1:Nfactor
        NvariableInRow = system.coefficientTranspose.colptr[i + 1] - system.coefficientTranspose.colptr[i]
        if NvariableInRow == 1
            variable = system.coefficientTranspose.rowval[system.coefficientTranspose.colptr[i]]
            meanDirect[variable] = 0.0
            weightDirect[variable] = 0.0
        else
            Nlink += NvariableInRow
            Nindirect += 1
            for j = system.coefficientTranspose.colptr[i]:(system.coefficientTranspose.colptr[i + 1] - 1)
                row = system.coefficientTranspose.rowval[j]
                push!(rowForward[i], row)
                push!(rowBackward[i], row)
            end
        end
    end

    ## Pass through the columns
    iterateVariable = fill(0, Nvariable)
    counter = 0
    colForward = [Int[] for i = 1:Nvariable]; colBackward = [Int[] for i = 1:Nvariable]
    incomingToFactor = [Int[] for i = 1:Nfactor]
    @inbounds for col = 1:Nvariable
        for i = system.coefficient.colptr[col]:(system.coefficient.colptr[col + 1] - 1)
            row = system.coefficient.rowval[i]
            NvariableInRow = system.coefficientTranspose.colptr[row + 1] - system.coefficientTranspose.colptr[row]
            if NvariableInRow == 1
                meanDirect[col] += system.observation[row] * system.coefficient[row, col] / system.variance[row]
                weightDirect[col] += system.coefficient[row, col]^2 / system.variance[row]
            else
                push!(colForward[col], row)
                push!(colBackward[col], row)
            end
        end
        if length(colForward[col]) == 1 && col != root
            counter += 1
            iterateVariable[counter] = col
        end
    end
    resize!(iterateVariable, counter)

    ### Initialize data
    fromVariable = fill(0, Nlink)
    toFactor = fill(0, Nlink)
    meanVariableFactor = fill(0.0, Nlink)
    varianceVariableFactor = fill(0.0, Nlink)

    fromFactor = fill(0, Nlink)
    toVariable = fill(0, Nlink)
    meanFactorVariable = fill(0.0, Nlink)
    varianceFactorVariable = fill(0.0, Nlink)

    mean = fill(0, Nvariable)
    variance = fill(0, Nvariable)

    passFactorVariable = 0
    passVariableFactor = 0
    iterateFactor = Int64[]

    return ContinuousTreeGraph(Nvariable, Nfactor, Nlink, root, virtualMean, virtualVariance, meanDirect, weightDirect,
            rowForward, rowBackward, colForward, colBackward, incomingToFactor, incomingToVariable,
            iterateFactor, iterateVariable, passFactorVariable, passVariableFactor, true, true),
           ContinuousInference(fromFactor, toVariable, meanFactorVariable, varianceFactorVariable,
            fromVariable, toFactor, meanVariableFactor, varianceVariableFactor, mean, variance)
end

