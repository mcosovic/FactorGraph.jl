mutable struct SystemModel
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

mutable struct FactorGraphTree
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

struct GraphicalModelTree
    graph::FactorGraphTree
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

########## Produce the graphical model ##########
function makeGraph(system, meanVirtual, varianceVirtual, dampProbability, dampAlpha)
    ### Number of factor and variable nodes
    Nfactor, Nvariable = size(system.jacobian)

    ### Find graph numbers, set the direct mean and variance and internal factor numeration
    Ndirect = 0; Nlink = 0; Nindirect = 0
    meanDirect = fill(meanVirtual / varianceVirtual, Nvariable)
    weightDirect = fill(1 / varianceVirtual, Nvariable)
    dynamic = fill(0, Nfactor)

    @inbounds for i = 1:Nfactor
        NvariableInRow = system.jacobianTranspose.colptr[i + 1] - system.jacobianTranspose.colptr[i]
        if NvariableInRow == 1
            Ndirect += NvariableInRow
            variable = system.jacobianTranspose.rowval[system.jacobianTranspose.colptr[i]]
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
        for i = system.jacobian.colptr[col]:(system.jacobian.colptr[col + 1] - 1)
            row = system.jacobian.rowval[i]
            NvariableInRow = system.jacobianTranspose.colptr[row + 1] - system.jacobianTranspose.colptr[row]
            if NvariableInRow == 1
                meanDirect[col] += system.observation[row] * system.jacobian[row, col] / system.variance[row]
                weightDirect[col] += system.jacobian[row, col]^2 / system.variance[row]
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
            for i = system.jacobianTranspose.colptr[col]:(system.jacobianTranspose.colptr[col + 1] - 1)
                row = system.jacobianTranspose.rowval[i]

                coefficient[idxi] = system.jacobianTranspose[row, col]
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

    return FactorGraph(Nvariable, Nfactor, Nindirect, Nlink,
                        meanVirtual, varianceVirtual, meanDirect, weightDirect, meanIndirect, varianceIndirect, coefficient,
                        sendToVariable, sendToFactor, rowptr, colptr, colptrMarginal, alphaNew, alphaOld,
                        iterateFactor, iterateVariable, iterateMarginal, dynamic),
           Inference(fromFactor, toVariable, meanFactorVariable, varianceFactorVariable,
                     fromVariable, toFactor, meanVariableFactor, varianceVariableFactor, mean, variance)
end

########## Set damping parameters ##########
function damping!(gbp::GraphicalModel; prob::Float64 = 0.6, alpha::Float64 = 0.4)
    gbp.graph.alphaNew = fill(1.0, gbp.graph.Nlink)
    gbp.graph.alphaOld = fill(0.0, gbp.graph.Nlink)
    bernoulliSample = randsubseq(collect(1:gbp.graph.Nlink), prob)
    @inbounds for i in bernoulliSample
        gbp.graph.alphaNew[i] = 1.0 - alpha
        gbp.graph.alphaOld[i] = alpha
    end
end

########## Form a tree factor graph and initialize messages ##########
function graphicalModelTree(
    args...;
    mean::Float64 = 0.0,
    variance::Float64 = 1e10,
    root::Int64 = 1)

    if checkSystemInputs(args)
        system = juliaOut(args)
    else
        system = juliaIn(args)
    end

    graph, inference = graphicalModelTree(system, mean, variance, root)

    return GraphicalModelTree(graph, inference, system)
end

########## Produce the graphical model ##########
function graphicalModelTree(system, virtualMean, virtualVariance, root)
    ### Number of factor and variable nodes
    Nfactor, Nvariable = size(system.jacobian)

    ### Find graph numbers, set the direct mean and variance, set factor numeration and pass through the rows
    Nlink = 0; Nindirect = 0
    meanDirect = fill(virtualMean / virtualVariance, Nvariable)
    weightDirect = fill(1 / virtualVariance, Nvariable)

    rowForward = [Int[] for i = 1:Nfactor]; rowBackward = [Int[] for i = 1:Nfactor]
    incomingToVariable = [Int[] for i = 1:Nvariable]
    @inbounds for i = 1:Nfactor
        NvariableInRow = system.jacobianTranspose.colptr[i + 1] - system.jacobianTranspose.colptr[i]
        if NvariableInRow == 1
            variable = system.jacobianTranspose.rowval[system.jacobianTranspose.colptr[i]]
            meanDirect[variable] = 0.0
            weightDirect[variable] = 0.0
        else
            Nlink += NvariableInRow
            Nindirect += 1
            for j = system.jacobianTranspose.colptr[i]:(system.jacobianTranspose.colptr[i + 1] - 1)
                row = system.jacobianTranspose.rowval[j]
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
        for i = system.jacobian.colptr[col]:(system.jacobian.colptr[col + 1] - 1)
            row = system.jacobian.rowval[i]
            NvariableInRow = system.jacobianTranspose.colptr[row + 1] - system.jacobianTranspose.colptr[row]
            if NvariableInRow == 1
                meanDirect[col] += z[row] * system.jacobian[row, col] / v[row]
                weightDirect[col] += system.jacobian[row, col]^2 / v[row]
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

    return FactorGraphTree(Nvariable, Nfactor, Nlink, root, virtualMean, virtualVariance, meanDirect, weightDirect,
            rowForward, rowBackward, colForward, colBackward, incomingToFactor, incomingToVariable,
            iterateFactor, iterateVariable, passFactorVariable, passVariableFactor, true, true),
        Inference(fromFactor, toVariable, meanFactorVariable, varianceFactorVariable,
            fromVariable, toFactor, meanVariableFactor, varianceVariableFactor, mean, variance)
end

########## Check that the factor graph has a tree structure ##########
function isTree(gbp::Union{GraphicalModel, GraphicalModelTree})
    ## Pass through the rows
    rowptr = [Int[] for i = 1:gbp.graph.Nfactor]
    iterateFactor = Int64[]
    @inbounds for i = 1:gbp.graph.Nfactor
        if gbp.system.jacobianTranspose.colptr[i + 1] - gbp.system.jacobianTranspose.colptr[i] != 1
            for j = gbp.system.jacobianTranspose.colptr[i]:(gbp.system.jacobianTranspose.colptr[i + 1] - 1)
                row = gbp.system.jacobianTranspose.rowval[j]
                push!(rowptr[i], row)
            end
        end
    end

    ## Pass through the columns
    iterateVariable = fill(0, gbp.graph.Nvariable)
    colptr = [Int[] for col = 1:gbp.graph.Nvariable]
    counter = 0
    @inbounds for col = 1:gbp.graph.Nvariable
        for i = gbp.system.jacobian.colptr[col]:(gbp.system.jacobian.colptr[col + 1] - 1)
            row = gbp.system.jacobian.rowval[i]
            if gbp.system.jacobianTranspose.colptr[row + 1] - gbp.system.jacobianTranspose.colptr[row] != 1
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
    @inbounds while hasSingleVariable || hasSingleFactor
        for i in iterateVariable
            for (k, j) in enumerate(rowptr[colptr[i][1]])
                if i == j
                    deleteat!(rowptr[colptr[i][1]], k)
                    deleteat!(colptr[i], 1)
                end
            end
            if length(rowptr[colptr[i][1]]) == 1
                push!(iterateFactor, colptr[i][1])
            end
        end
        if isempty(iterateFactor)
            hasSingleFactor = false
        end
        iterateVariable = Int64[]

        for i in iterateFactor
            for (k, j) in enumerate(colptr[rowptr[i][1]])
                if i == j
                    deleteat!(colptr[rowptr[i][1]], k)
                    deleteat!(rowptr[i], 1)
                end
            end
            if length(colptr[rowptr[i][1]]) == 1
                push!(iterateVariable, rowptr[i][1])
            end
        end
        if isempty(iterateVariable)
            hasSingleVariable = false
        end
        iterateFactor = Int64[]
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