struct SystemModel
    jacobian::SparseMatrixCSC{Float64,Int64}
    jacobianTranspose::SparseMatrixCSC{Float64,Int64}
    observation::Array{Float64,1}
    variance::Array{Float64,1}
    data::String
end

mutable struct FactorGraph
    Nvariable::Int64
    Nfactor::Int64
    Nindirect::Int64
    Nlink::Int64
    meanDirect::Array{Float64,1}
    weightDirect::Array{Float64,1}
    meanIndirect::Array{Float64,1}
    varianceIndirect::Array{Float64,1}
    coefficient::Array{Float64,1}
    sendToVariable::Array{Int64,1}
    sendToFactor::Array{Int64,1}
    rowptr::Array{Int64,1}
    colptr::Array{Int64,1}
    alphaNew::Array{Float64,1}
    alphaOld::Array{Float64,1}
    iterateFactor::Array{Int64,1}
    iterateVariable::Array{Int64,1}
    dynamic::Array{Int64,1}
    ageing::Array{Int64,1}
end

mutable struct Inference
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

struct GraphicalModel
    graph::FactorGraph
    inference::Inference
    system::SystemModel
end

########## Form a graph and initialize messages ##########
function graphicalModel(
    args...;
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    mean::Float64 = 0.0,
    variance::Float64 = 1e10)

    checkKeywords(prob, alpha, variance)
    if checkSystemInputs(args)
        system = juliaOut(args)
    else
        system = juliaIn(args)
    end

    graph, inference = makeGraph(system, mean, variance, prob, alpha)

    return GraphicalModel(graph, inference, system)
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

    #### Check initial messages
    if variance < 0.0
        error("Invalid variance value.")
    end
end

########## Type of the system data ##########
function checkSystemInputs(args)
    fromfile = false
    @inbounds for i in args
        if typeof(i) == String
            fromfile = true
            break
        end
    end

    return fromfile
end

########## Load from HDF5 or XLSX files ##########
function juliaOut(args)
    #### Check the package is installed
    pathtoGaussBP = Base.find_package("GaussBP")
    if isnothing(pathtoGaussBP)
        throw(ErrorException("GaussBP not found in install packages"))
    end
    packagepath = abspath(joinpath(dirname(pathtoGaussBP), ".."))

    extension = ".h5"; path = ""; dataname = ""; fullpath = ""

    @inbounds for i = 1:length(args)
        try
            extension = string(match(r"\.[A-Za-z0-9]+$", args[i]).match)
        catch
            extension = ""
        end
        if extension == ".h5" || extension == ".xlsx"
            fullpath = args[i]
            path = dirname(args[i])
            dataname = basename(args[i])
            break
        end
    end

    if isempty(extension)
        throw(ErrorException("the input DATA extension is not found"))
    elseif extension != ".h5" && extension != ".xlsx"
        throw(DomainError(extension, "the input DATA extension is not supported"))
    end

    if path == ""
        path = joinpath(packagepath, "src/example/")
        fullpath = joinpath(packagepath, "src/example/", dataname)
    end

    if !(dataname in cd(readdir, path))
        throw(DomainError(dataname, "the input DATA is not found"))
    end

    #### Read from HDF5 or XLSX file
    if extension == ".h5"
        list = h5read(fullpath, "/jacobian")::Array{Float64,2}
        jacobian = sparse(list[:,1], list[:,2], list[:,3])::SparseMatrixCSC{Float64,Int64}
        jacobianTranspose = sparse(list[:,2], list[:,1], list[:,3])::SparseMatrixCSC{Float64,Int64}

        observation = h5read(fullpath, "/observation")::Array{Float64,1}
        variance = h5read(fullpath, "/variance")::Array{Float64,1}
    elseif extension == ".xlsx"
        xf = XLSX.openxlsx(fullpath, mode = "r")
        if "jacobian" in XLSX.sheetnames(xf)
            start = startxlsx(xf["jacobian"])
            list = xf["jacobian"][:][start:end, :]
            jacobian = sparse(list[:, 1], list[:, 2], list[:, 3])
            jacobianTranspose = sparse(list[:, 2], list[:, 1], list[:, 3])
        else
            throw(ErrorException("error opening sheet jacobian"))
        end
        if "measurement" in XLSX.sheetnames(xf)
            start = startxlsx(xf["measurement"])
            list = xf["measurement"][:][start:end, :]
            observation = list[:, 1]
            variance = list[:, 2]
        else
            throw(ErrorException("error opening sheet measurement"))
        end
    else
        error("the input data is not a valid format")
    end

    return SystemModel(jacobian, jacobianTranspose, observation, variance, dataname)
end

########## Read in-Julia system model ##########
function juliaIn(args)
    if typeof(args[1]) == Array{Float64, 2}
        jacobian = sparse(args[1])
    else
        jacobian = args[1]
    end

    jacobianTranspose = copy(transpose(jacobian))
    observation = args[2]
    variance = args[3]

    return SystemModel(jacobian, jacobianTranspose, observation, variance, "noname")
end

########## Start row in xlsx-file ##########
function startxlsx(xf)
    start = 1
    @inbounds for r in XLSX.eachrow(xf)
        if !isa(r[1], String)
            start = XLSX.row_number(r)
            break
        end
    end
    return start
end

########## Produce the factor graph and define singly-connected, virtual and indirect factor nodes ##########
function makeGraph(system, meanVirtual, varianceVirtual, dampProbability, dampAlpha)
    Nfactor, Nvariable = size(system.jacobian)
    Ndirect = 0
    Nlink = 0
    NvariableInColumn = 0
    direct = fill(0, Nvariable)

    ### Find number of links
    @inbounds for i = 1:Nfactor
        NvariableInColumn = system.jacobianTranspose.colptr[i + 1] - system.jacobianTranspose.colptr[i]

        if NvariableInColumn == 1
            Ndirect += NvariableInColumn
            row = system.jacobianTranspose.rowval[nzrange(system.jacobianTranspose, i)][1]
            direct[row] = 1
        else
            Nlink += NvariableInColumn
        end
    end

    ### Define singly-connected and indirect factor nodes arrays and dynamic data
    Nindirect = Nfactor - Ndirect
    meanDirect = fill(0.0, Nvariable)
    weightDirect = fill(0.0, Nvariable)

    factorIndex = fill(0, Nlink)
    variableIndex = similar(factorIndex)
    coefficient = fill(0.0, Nlink)
    meanIndirect = fill(0.0, Nlink)
    varianceIndirect = similar(meanIndirect)
    dynamic = fill(0, Nfactor)

    rowptr = fill(0, Nindirect)
    colptr = fill(0, Nvariable)
    idxi = 1; idxr = 1
    idxT = findall(!iszero, system.jacobianTranspose)
    @inbounds for i in idxT
        if (system.jacobianTranspose.colptr[i[2] + 1] - system.jacobianTranspose.colptr[i[2]]) == 1
            meanDirect[i[1]] += system.observation[i[2]] * system.jacobianTranspose[i] / system.variance[i[2]]
            weightDirect[i[1]] += system.jacobianTranspose[i]^2 / system.variance[i[2]]
        else
            factorIndex[idxi] = idxr
            variableIndex[idxi] = i[1]

            coefficient[idxi] = system.jacobianTranspose[i]

            meanIndirect[idxi] = system.observation[i[2]]
            varianceIndirect[idxi] = system.variance[i[2]]

            idxi += 1
            colptr[i[1]] += 1
            rowptr[idxr] = idxi

            if idxT[system.jacobianTranspose.colptr[i[2] + 1] - 1] == i
                dynamic[i[2]] = idxr
                idxr += 1
            end
        end
        if direct[i[1]] == 0
            meanDirect[i[1]] = meanVirtual / varianceVirtual
            weightDirect[i[1]] = 1 / varianceVirtual
            direct[i[1]] = 1
        end
    end
    pushfirst!(colptr, 1)
    pushfirst!(rowptr, 1)
    colptr = cumsum(colptr)

    ### Indices for sending messages
    temp = collect(1:Nlink)
    sortVariableIndex = sortperm(variableIndex)
    sendToFactor = temp[sortVariableIndex]

    newFactorIndex = factorIndex[sortVariableIndex]
    sortFactorIndex = sortperm(newFactorIndex)
    sendToVariable = temp[sortFactorIndex]

    ### Indices for messages inference output
    lookup = Dict{Int, Int}()
    for (k, i) in enumerate(dynamic)
        if i != 0
            lookup[i] = k
        end
    end

    factorIdx = similar(factorIndex)
    @inbounds for i = 1:Nindirect
        for j in rowptr[i]:(rowptr[i + 1] - 1)
            factorIdx[sendToVariable[j]] = lookup[i]
            factorIndex[j] = lookup[i]
        end
    end
    toVariable = variableIndex[sendToFactor]


    ### Initialize arrays for messages
    meanFactorVariable = similar(coefficient)
    varianceFactorVariable = similar(coefficient)
    meanVariableFactor = similar(coefficient)
    varianceVariableFactor = similar(coefficient)

    ### Set damping parameters
    alphaNew = fill(1.0, Nlink)
    alphaOld = fill(0.0, Nlink)
    bernoulliSample = randsubseq(collect(1:Nlink), dampProbability)
    @inbounds for i in bernoulliSample
        alphaNew[i] = 1.0 - dampAlpha
        alphaOld[i] = dampAlpha
    end

    ### Pass messages from singly-connected factor nodes to all indirect links
    @inbounds for i = 1:Nvariable
        for j in colptr[i]:(colptr[i + 1] - 1)
            k = sendToFactor[j]
            varianceVariableFactor[k] = 1 / weightDirect[i]
            meanVariableFactor[k] = meanDirect[i] * varianceVariableFactor[k]
        end
    end

    ### Initialize marginal mean and variance vectors
    mean = fill(0.0, Nvariable)
    variance = fill(0.0, Nvariable)

    ### Iteration counters
    ageing = [0]
    iterateFactor = collect(1:Nindirect)
    iterateVariable = collect(1:Nvariable)


    return FactorGraph(Nvariable, Nfactor, Nindirect, Nlink, meanDirect, weightDirect, meanIndirect, varianceIndirect, coefficient, sendToVariable, sendToFactor, rowptr, colptr, alphaNew, alphaOld, iterateFactor, iterateVariable, dynamic, ageing),
           Inference(factorIdx, toVariable, meanFactorVariable, varianceFactorVariable, variableIndex, factorIndex, meanVariableFactor, varianceVariableFactor, mean, variance)
end

########## Set damping parameters ##########
function damping!(gbp::GraphicalModel; prob::Float64 = 0.6, alpha::Float64 = 0.4)
    bernoulliSample = randsubseq(collect(1:gbp.graph.Nlink), prob)
    @inbounds for i in bernoulliSample
        gbp.graph.alphaNew[i] = 1.0 - alpha
        gbp.graph.alphaOld[i] = alpha
    end
end

########## Dynamic the GBP update ##########
@inline function dynamicInference!(gbp::GraphicalModel, dynamic)
    factor = trunc.(Int, dynamic[:, 1])
    @inbounds for (k, i) in enumerate(factor)
        if (gbp.system.jacobianTranspose.colptr[i + 1] - gbp.system.jacobianTranspose.colptr[i]) == 1
            idx = gbp.system.jacobianTranspose.colptr[i]
            variable = trunc(Int, gbp.system.jacobianTranspose.rowval[idx])

            gbp.graph.meanDirect[variable] -= gbp.system.observation[i] * gbp.system.jacobianTranspose[variable, i] / gbp.system.variance[i]
            gbp.graph.weightDirect[variable] -= gbp.system.jacobianTranspose[variable, i]^2 / gbp.system.variance[i]

            gbp.graph.meanDirect[variable] += dynamic[k, 2] * gbp.system.jacobianTranspose[variable, i] / dynamic[k, 3]
            gbp.graph.weightDirect[variable] += gbp.system.jacobianTranspose[variable, i]^2 / dynamic[k, 3]
        else
            j = gbp.graph.dynamic[i]
            @inbounds for n in gbp.graph.rowptr[j]:(gbp.graph.rowptr[j + 1] - 1)
                gbp.graph.meanIndirect[n] = dynamic[k, 2]
                gbp.graph.varianceIndirect[n] = dynamic[k, 3]
            end
        end
        gbp.system.observation[i] = dynamic[k, 2]
        gbp.system.variance[i] = dynamic[k, 3]
    end
    gbp.graph.ageing[1] = 0
end

######### Ageing the GBP update ##########
@inline function ageingInference!(gbp::GraphicalModel, dynamic)
    gbp.graph.ageing[1] += 1
    factor = trunc.(Int, dynamic[:, 1])
    @inbounds for (k, i) in enumerate(factor)
        if gbp.system.variance[i] < dynamic[k, 7]
            if (gbp.system.jacobianTranspose.colptr[i + 1] - gbp.system.jacobianTranspose.colptr[i]) == 1
                idx = gbp.system.jacobianTranspose.colptr[i]
                variable = trunc(Int, gbp.system.jacobianTranspose.rowval[idx])

                gbp.graph.meanDirect[variable] -= gbp.system.observation[i] * gbp.system.jacobianTranspose[variable, i] / gbp.system.variance[i]
                gbp.graph.weightDirect[variable] -= gbp.system.jacobianTranspose[variable, i]^2 / gbp.system.variance[i]
            end
            if dynamic[k, 4] == 1
                gbp.system.variance[i] = dynamic[k, 5] * gbp.graph.ageing[1] + dynamic[k, 3]
            elseif dynamic[k, 4] == 2
                d = 1 + dynamic[k, 6]
                gbp.system.variance[i] = dynamic[k, 5] * log10((gbp.graph.ageing[1] + d) / d) + dynamic[k, 3]
            elseif dynamic[k, 4] == 3
                d = 1 + dynamic[k, 6]
                gbp.system.variance[i] = dynamic[k, 3] * d^(dynamic[k, 5] * gbp.graph.ageing[1])
            end
            if gbp.system.variance[i] > dynamic[k, 7]
                gbp.system.variance[i] =  dynamic[k, 7]
            end

            if (gbp.system.jacobianTranspose.colptr[i + 1] - gbp.system.jacobianTranspose.colptr[i]) == 1
                gbp.graph.meanDirect[variable] += gbp.system.observation[i] * gbp.system.jacobianTranspose[variable, i] / gbp.system.variance[i]
                gbp.graph.weightDirect[variable] += gbp.system.jacobianTranspose[variable, i]^2 / gbp.system.variance[i]
            else
                j = gbp.graph.dynamic[i]
                @inbounds for n in gbp.graph.rowptr[j]:(gbp.graph.rowptr[j + 1] - 1)
                    gbp.graph.varianceIndirect[n] = gbp.system.variance[i]
                end
            end
        end
    end
end

######### Freeze factor node ##########
function freezeFactor!(gbp; factor = 0)
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
    for i = 1:length(gbp.graph.iterateFactor)
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
function defreezeFactor!(gbp; factor = 0)
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
    for i = 1:length(gbp.graph.iterateFactor)
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
function freezeVariable!(gbp; variable = 0)
    if variable == 0
        error("The keyword variable is missing.")
    elseif variable > gbp.graph.Nvariable
        error("The variable node does not exist.")
    end

    whereIs = 0
    for i = 1:length(gbp.graph.iterateVariable)
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

######### Defreeze variable node ##########
function defreezeVariable!(gbp; variable = 0)
    if variable == 0
        error("The keyword variable is missing.")
    elseif variable > gbp.graph.Nvariable
        error("The variable node does not exist.")
    end

    whereIs = 0
    for i = 1:length(gbp.graph.iterateVariable)
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