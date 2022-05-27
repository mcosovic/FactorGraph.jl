mutable struct DiscreteSystem
    probability::Vector{Vector{Int64}}
    table::Vector{Array{Float64, N} where N}
    incidence::SparseMatrixCSC{Float64,Int64}
    incidenceTranspose::SparseMatrixCSC{Float64,Int64}
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


    system = readDiscreteArguments(args)
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

    incidence, incidenceTranspose = incidenceMatrix(probability)

    return DiscreteSystem(probability, table, incidence, incidenceTranspose)
end

########## Incidence matrix ##########
function incidenceMatrix(probabilities)
    row = Int64[]; col = Int64[]

    @inbounds for (k, probability) in enumerate(probabilities)
        for variable in probability
            push!(row, k)
            push!(col, variable)
        end
    end
    incidence = sparse(row, col, 1)
    incidenceTranspose = sparse(col, row, 1)

    return incidence, incidenceTranspose
end

########## Produce the graphical model ##########
function makeDiscreteTreeGraph(system, virtualMessage, root)
    ### Number of factor and variable nodes
    Nfactor, Nvariable = size(system.incidence)

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