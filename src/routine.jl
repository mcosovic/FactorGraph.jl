struct FactorGraph
    row::Array{Int64,1}
    col::Array{Int64,1}
    rowptr::Array{Int64,1}
    colptr::Array{Int64,1}
    Mind::Array{Float64,1}
    Vind::Array{Float64,1}
    coeff::Array{Float64,1}
    coeffInv::Array{Float64,1}
    Mdir::Array{Float64,1}
    Wdir::Array{Float64,1}
    Nvar::Int64
    Nfac::Int64
    Ndir::Int64
    Nind::Int64
    Nlink::Int64
    Ndyn::Int64
    dirDyn::Dict{Int64,Array{Int64,1}}
    facDyn::Dict{Int64,Array{Int64,1}}
    observationDyn::Array{Float64,1}
    varianceDyn::Array{Float64,1}
end

struct BeliefPropagation
    iterNative::Int64
    iterDamp::Int64
    iterBump::Int64
    iterDampBump::Int64
    iterCount::Array{Int64,1}
    flagDyn::Array{Int64,1}
    updateCount::Array{Int64,1}
    Mfac_var::Array{Float64,1}
    Wfac_var::Array{Float64,1}
    Mvar_fac::Array{Float64,1}
    Vvar_fac::Array{Float64,1}
    alpha1::Array{Float64,1}
    alpha2::Array{Float64,1}
    to_fac::Array{Int64,1}
    to_var::Array{Int64,1}
end


########## Produce the factor graph and define singly-connected, virtual and indirect factor nodes ##########
function factors(system, settings)
    Nfac, Nvar = size(system.jacobian)
    Ndir = 0
    Nlink = 0
    var_in_column = 0

    ### Find number of links
    @inbounds for i = 1:Nfac
        var_in_column = system.jacobianTranspose.colptr[i + 1] - system.jacobianTranspose.colptr[i]

        if var_in_column == 1
            Ndir += var_in_column
        else
            Nlink += var_in_column
        end
    end

    ### Define singly-connected and indirect factor nodes arrays and dynamic data
    Nind = Nfac - Ndir
    Mdir = fill(0.0, Nvar)
    Wdir = fill(0.0, Nvar)

    row = fill(0, Nlink)
    col = similar(row)
    coeff = fill(0.0, Nlink)
    coeffInv = similar(coeff)
    Mind = fill(0.0, Nlink)
    Vind = similar(Mind)

    rowptr = fill(0, Nind)
    colptr = fill(0, Nvar)
    idxi = 1
    idxr = 1

    if settings.dynamic
        Ndyn, ~ = size(system.dynamic)
        dirDyn = Dict(node => Int64[] for node in collect(1:Nvar))
        facDyn = Dict(node => Int64[] for node in collect(1:Nfac))  
        observationDyn = copy(system.observation)
        varianceDyn = copy(system.variance) 
    else
        Ndyn = 0
        dirDyn = Dict(1 => Int64[])
        facDyn = Dict(1 => Int64[])    
        observationDyn = [0.0]
        varianceDyn = [0.0]
    end
 
    idxT = findall(!iszero, system.jacobianTranspose)
    @inbounds for i in idxT
        if (system.jacobianTranspose.colptr[i[2] + 1] - system.jacobianTranspose.colptr[i[2]]) == 1
            Mdir[i[1]] += system.observation[i[2]] * system.jacobianTranspose[i] / system.variance[i[2]]
            Wdir[i[1]] += (system.jacobianTranspose[i]^2) / system.variance[i[2]]

            if settings.dynamic
                push!(dirDyn[i[1]], i[2])
                facDyn[i[2]] = [1; i[1]]
            end
        else
            row[idxi] = idxr
            col[idxi] = i[1]

            coeff[idxi] = system.jacobianTranspose[i]
            coeffInv[idxi] = 1 / system.jacobianTranspose[i]

            Mind[idxi] = system.observation[i[2]]
            Vind[idxi] = system.variance[i[2]]

            idxi += 1

            colptr[i[1]] += 1
            rowptr[idxr] = idxi

            if idxT[system.jacobianTranspose.colptr[i[2] + 1] - 1] == i
                if settings.dynamic
                    facDyn[i[2]] = [2; idxr]
                end
                idxr += 1
            end
        end
    end
    pushfirst!(colptr, 1)
    pushfirst!(rowptr, 1)
    colptr = cumsum(colptr)

    return FactorGraph(row, col, rowptr, colptr, Mind, Vind, coeff, coeffInv, Mdir, Wdir, Nvar, Nfac, Ndir, Nind, Nlink, Ndyn, dirDyn, facDyn, observationDyn, varianceDyn)
end


########## Initialize the GBP algorithm ##########
function initialize(graph, settings)
    ### Desing iteration scheme 
    dabu = minimum([settings.bump, settings.damp])
    iterNative = 0
    if settings.bump == settings.damp == settings.maxIter
        iterNative = settings.maxIter
    else
        iterNative = dabu - 1
    end

    iterDamp = 0
    if settings.bump == settings.maxIter
        iterDamp = settings.maxIter - iterNative
    elseif settings.bump > settings.damp
        iterDamp = settings.bump - settings.damp
    end

    iterBump = 0
    if settings.damp == settings.maxIter
        iterBump = settings.maxIter - iterNative
    elseif settings.bump < settings.damp
        iterBump = settings.damp - settings.bump
    end

    iterDampBump = 0
    if settings.damp < settings.maxIter && settings.bump < settings.maxIter
       iterDampBump = settings.maxIter - iterDamp - iterBump - iterNative
    end

    if !settings.outIterate && !settings.dynamic
        iterCount = [settings.maxIter]
    else
        iterCount = [0]
    end
    flagDyn = [1]   
    updateCount = [1]

    ### Initialize arrays for messages
    Mfac_var = similar(graph.coeff)
    Wfac_var = similar(graph.coeff)
    Mvar_fac = fill(settings.mean, graph.Nlink)
    Vvar_fac = fill(settings.variance, graph.Nlink)

    ### Set damping parameters
    bernoulli_sample = randsubseq(collect(1:graph.Nlink), settings.prob)
    alpha1 = fill(1.0, graph.Nlink)
    alpha2 = fill(0.0, graph.Nlink)

    @inbounds for i in bernoulli_sample
        alpha1[i] = 1.0 - settings.alpha
        alpha2[i] = settings.alpha
    end

    ### Indices for sending messages
    temp = collect(1:graph.Nlink)
    sort_col_idx = sortperm(graph.col)
    to_fac = temp[sort_col_idx]

    new_row = graph.row[sort_col_idx]
    sort_row_idx = sortperm(new_row)
    to_var = temp[sort_row_idx]

    ### Pass messages from singly-connected factor nodes to all corresponding indirect links
    idx = findall(!iszero, graph.Wdir)
    @inbounds for i in idx
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            k = to_fac[j]
            Vvar_fac[k] = 1 / graph.Wdir[i]
            Mvar_fac[k] = graph.Mdir[i] * Vvar_fac[k]
        end
    end

    return BeliefPropagation(iterNative, iterDamp, iterBump, iterDampBump, iterCount, flagDyn, updateCount, Mfac_var, Wfac_var, Mvar_fac, Vvar_fac, alpha1, alpha2, to_fac, to_var)
end

########## Dynamic the GBP update ##########
@inline function graph_dynamic(settings, system, graph, bp)
    bp.iterCount[1] += 1 

    while graph.Ndyn >= bp.flagDyn[1] && trunc(Int, system.dynamic[bp.flagDyn[1], 1]) == bp.iterCount[1]  
        u = bp.flagDyn[1]
        factor = trunc(Int, system.dynamic[u, 2])
        mean = system.dynamic[u, 3] 
        variance = system.dynamic[u, 4] 
        
        if settings.outDisplay
            println("\n Update Factor Nodes")
            A = [system.dynamic[u, 1:2]; system.dynamic[u, 3:4]; graph.observationDyn[factor]; graph.varianceDyn[factor]]
            pretty_table(A', header = ["Iteration", "Factor Node", "New Mean", "New Variance", "Old Mean", "Old Variance"],
            alignment=[:r, :r, :r, :r, :r, :r], formatters = ft_printf(["%3.0f", "%3.0f","%3.3f", "%3.3e", "%3.3f", "%3.3e"], [1, 2, 3, 4, 5, 6]))
        end

        graph.observationDyn[factor] = mean
        graph.varianceDyn[factor] = variance

        if graph.facDyn[factor][1] == 2
            i = graph.facDyn[factor][2]
            @inbounds for k in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
                graph.Mind[k] = mean
                graph.Vind[k] = variance
            end
        else
            variable = graph.facDyn[factor][2]
            update_factors = graph.dirDyn[variable]

            graph.Mdir[variable] = mean * system.jacobianTranspose[variable, factor] / variance
            graph.Wdir[variable] = system.jacobianTranspose[variable, factor]^2 / variance

            @inbounds for i in update_factors
                if i != factor
                    graph.Mdir[variable] += graph.observationDyn[i] * system.jacobianTranspose[variable, i] / graph.varianceDyn[i]
                    graph.Wdir[variable] += (system.jacobianTranspose[variable, i]^2) / graph.varianceDyn[i]
                end
            end
        end
        bp.flagDyn[1] += 1
    end
end

########## Start row in xlsx-file ##########
function startxlsx(xf)
    start = 1
    for r in XLSX.eachrow(xf)
        if !isa(r[1], String)
            start = XLSX.row_number(r)
            break
        end
    end
    return start
end