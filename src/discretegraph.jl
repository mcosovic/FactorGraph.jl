mutable struct DiscreteSystem
    probability::Vector{Vector{Int64}}
    table::Vector{Array{Float64, N} where N}
    jacobian::SparseMatrixCSC{Float64,Int64}
    jacobianTranspose::SparseMatrixCSC{Float64,Int64}
    data::String
end

mutable struct DiscreteTreeGraph
    Nvariable::Int64
    Nfactor::Int64
    Nlink::Int64
    root::Int64
    virtualMessage::Float64
    state::Array{Int64,1}
    messageDirect::Vector{Vector{Float64}}
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

mutable struct DiscreteInference
    fromFactor::Array{Int64,1}
    toVariable::Array{Int64,1}
    messageFactorVariable::Vector{Vector{Float64}}
    fromVariable::Array{Int64,1}
    toFactor::Array{Int64,1}
    messageVariableFactor::Vector{Vector{Float64}}
    marginal::Vector{Vector{Float64}}
end

struct DiscreteTreeModel
    graph::DiscreteTreeGraph
    inference::DiscreteInference
    system::DiscreteSystem
end

########## Form a graph with continuous variables and initialize messages ##########
function discreteTreeModel(
    args...;
    virtual::Float64 = 1.0,
    root::Int64 = 1)

    if checkFileOrArguments(args)
        system = readDiscreteFile(args)
    else
        system = readDiscreteArguments(args)
    end

    graph, inference = makeDiscreteTreeGraph(system, virtual, root)

    return DiscreteTreeModel(graph, inference, system)
end

########## Read in-Julia discrete system model ##########
function readDiscreteArguments(args)
    if isa(args[1], Dict)
        probability = []; table = []
        @inbounds for (key, value) in args[1]
            if key == "probability"
                for k in value
                    push!(probability, k)
                end
            end
            if key == "table"
                for k in value
                    push!(table, k)
                end
            end
        end
    else
        probability = args[1]
        table = args[2]
    end

    jacobian, jacobianTranspose = incidenceMatrix(probability)

    return DiscreteSystem(probability, table, jacobian, jacobianTranspose, "noname")
end

########## Load from HDF5 or XLSX files ##########
function readDiscreteFile(args)
    fullpath, extension, dataname = checkImportFile(args)

    #### Read from HDF5 or XLSX file
    if extension == ".h5"
        list = h5open(fullpath, "r")
        probability = []; table = Any[]
        probabilityH5 = h5read(fullpath, "/probability")
        tableH5 = h5read(fullpath, "/table")
        for key in sort(collect(keys(probabilityH5)))
            push!(probability, vec(probabilityH5[key]))
            if size(tableH5[key], 2) == 1
                push!(table, vec(tableH5[key]))
            else
                push!(table, tableH5[key])
            end
        end
    elseif extension == ".xlsx"
        xf = XLSX.openxlsx(fullpath, mode = "r")
        if "graphical model" in XLSX.sheetnames(xf)
            list = xf["graphical model"][:][1:end, :]

            probability = Any[]; table = Any[]; cpt = []
            Nlist = size(list, 1);  flag = false; last = 0
            @inbounds for i = 1:Nlist - 1
                if !ismissing(list[i, 1])
                    if list[i, 1] == "conditional independence"
                        for (q, k) in enumerate(list[i, 2:end])
                            if ismissing(k)
                                last = q
                                break
                            end
                        end
                        push!(probability, list[i, 2:last])
                    end

                    if list[i, 1] == "discrete variable states"
                        cpt = zeros(tuple(list[i, 2:last]...))
                    end
                    if list[i, 1] == "conditional probability table"
                        flag = true
                    end
                end
                if flag && !ismissing(list[i, 2])
                    idx = CartesianIndex(tuple(list[i, 2:last]...))
                    cpt[idx] = list[i, last + 1]
                end
                if !ismissing(list[i + 1, 1]) && list[i + 1, 1] == "conditional independence"
                    push!(table, cpt)
                    flag = false
                end
                if i + 1 == Nlist && !ismissing(list[i + 1, 2])
                    idx = CartesianIndex(tuple(list[i + 1, 2:last]...))
                    cpt[idx] = list[i + 1, last + 1]
                    push!(table, cpt)
                end
            end
        else
            throw(ErrorException("Error opening sheet: graphical model."))
        end
    else
        error("The input data is not a valid format.")
    end

    jacobian, jacobianTranspose = incidenceMatrix(probability)

    return DiscreteSystem(probability, table, jacobian, jacobianTranspose, dataname)
end

########## Form incidence matrix ##########
function incidenceMatrix(probabilities)
    row = Int64[]; col = Int64[]

    @inbounds for (k, probability) in enumerate(probabilities)
        for variable in probability
            push!(row, k)
            push!(col, variable)
        end
    end
    jacobian = sparse(row, col, 1)
    jacobianTranspose = sparse(col, row, 1)

    return jacobian, jacobianTranspose
end

########## Produce the graphical model ##########
function makeDiscreteTreeGraph(system, virtualMessage, root)
    ### Number of factor and variable nodes
    Nfactor, Nvariable = size(system.jacobian)

    ### Pass through probabilities
    Nlink = 0
    states = fill(0, Nvariable)
    messageDirect = [Float64[] for i = 1:Nvariable]
    rowForward =  [Int64[] for i = 1:Nfactor]
    colForward = [Int64[] for i = 1:Nvariable]
    @inbounds for (k, probability) in enumerate(system.probability)
        for (kk, variable) in enumerate(probability)
            if states[variable] == 0
                states[variable] = size(system.table[k], kk)
                messageDirect[variable] = fill(virtualMessage, states[variable])
            elseif states[variable] != size(system.table[k], kk)
                error("The possible states of the variable $variable do not agree through tables.")
            end
            if length(probability) != 1
                push!(rowForward[k], variable)
                push!(colForward[variable], k)
                Nlink += 1
            else
                messageDirect[variable] .*= system.table[k]
            end
        end
    end
    rowBackward = deepcopy(rowForward)
    colBackward = deepcopy(colForward)

    ### Initialize the first forward step
    iterateVariable = fill(0, Nvariable)
    counter = 0
    @inbounds for col = 1:Nvariable
        if length(colForward[col]) == 1 && col != root
            counter += 1
            iterateVariable[counter] = col
        end
    end
    resize!(iterateVariable, counter)

    ### Initialize data
    fromVariable = fill(0, Nlink)
    toFactor = fill(0, Nlink)
    messageVariableFactor = [Float64[] for i = 1:Nlink]

    fromFactor = fill(0, Nlink)
    toVariable = fill(0, Nlink)
    messageFactorVariable = [Float64[] for i = 1:Nlink]

    marginal = [Float64[] for i = 1:Nvariable]
    incomingToVariable = [Int64[] for i = 1:Nvariable]
    incomingToFactor = [zeros(length(system.probability[i])) for i = 1:Nfactor]
    passFactorVariable = 0
    passVariableFactor = 0
    iterateFactor = Int64[]

    return DiscreteTreeGraph(Nvariable, Nfactor, Nlink, root, virtualMessage, states, messageDirect,
        rowForward, rowBackward, colForward, colBackward, incomingToFactor, incomingToVariable,
        iterateFactor, iterateVariable, passFactorVariable, passVariableFactor, true, true),
           DiscreteInference(fromFactor, toVariable, messageFactorVariable, fromVariable, toFactor, messageVariableFactor, marginal)
end