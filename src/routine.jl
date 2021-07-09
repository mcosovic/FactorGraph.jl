### Struct variables
struct FactorGraph
    row::Array{Int64,1}
    rowptr::Array{Int64,1}
    col::Array{Int64,1}
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
end

struct BeliefPropagation
    Mfac_var::Array{Float64,1}
    Wfac_var::Array{Float64,1}
    Mvar_fac::Array{Float64,1}
    Vvar_fac::Array{Float64,1}
    alpha1::Array{Float64,1}
    alpha2::Array{Float64,1}
    to_fac::Array{Int64,1}
    to_var::Array{Int64,1}
    iterCount::Array{Int64,1}
end

struct Results
    iterations::Union{Int64, Array{Int64,1}}
    means::Union{Array{Float64,2}, Array{Float64,1}} 
    variances::Union{Array{Float64,2}, Array{Float64,1}} 
end

struct ResultsError
    iterations::Union{Int64, Array{Int64,1}}
    means::Union{Array{Float64,2}, Array{Float64,1}} 
    variances::Union{Array{Float64,2}, Array{Float64,1}} 
    rmse::Array{Float64,1}
    mae::Array{Float64,1}
    wrss::Array{Float64,1}
end

struct ResultsWLS
    iterations::Union{Int64, Array{Int64,1}}
    means::Union{Array{Float64,2}, Array{Float64,1}} 
    variances::Union{Array{Float64,2}, Array{Float64,1}} 
    meansWLS::Array{Float64,1}
    rmseWLS::Array{Float64,1}
    maeWLS::Array{Float64,1}
    wrssWLS::Array{Float64,1}
    rmseBPWLS::Array{Float64,1}
    maeBPWLS::Array{Float64,1}  
end

struct ResultsErrorWLS
    iterations::Union{Int64, Array{Int64,1}}
    means::Union{Array{Float64,2}, Array{Float64,1}} 
    variances::Union{Array{Float64,2}, Array{Float64,1}} 
    rmse::Array{Float64,1} 
    mae::Array{Float64,1}
    wrss::Array{Float64,1}
    meansWLS::Array{Float64,1}
    rmseWLS::Array{Float64,1}
    maeWLS::Array{Float64,1}
    wrssWLS::Array{Float64,1}
    rmseBPWLS::Array{Float64,1}
    maeBPWLS::Array{Float64,1}  
end

struct Recursion
    Mcol::Array{Float64,1}
    Wcol::Array{Float64,1}
end

### Produce the factor graph and define singly-connected, virtual and indirect factor nodes
function factors(system, settings)
    Nfac, Nvar = size(system.J)
    Ndir = 0
    Nlink = 0
    var_in_column = 0
    dir_position = fill(0, Nfac)

    ########## Find number of links ##########
    @inbounds for i = 1:Nfac
        var_in_column = system.Jt.colptr[i + 1] - system.Jt.colptr[i]

        if var_in_column == 1
            Ndir += var_in_column
            dir_position[i] = system.Jt.rowval[system.Jt.colptr[i]]
        else
            Nlink += var_in_column
        end
    end

    ########## Define position of virtual factor nodes ##########
    idx = findall(!iszero, dir_position)
    virtual = fill(1, Nvar)

    @inbounds for i in idx
        virtual[dir_position[i]] = 0
    end

    ########## Define singly-connected, virtual and indirect factor nodes arrays ##########
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

    idxT = findall(!iszero, system.Jt)
    @inbounds for i in idxT
        if (system.Jt.colptr[i[2] + 1] - system.Jt.colptr[i[2]]) == 1
            Mdir[i[1]] += system.z[i[2]] * system.Jt[i] / system.v[i[2]]
            Wdir[i[1]] += (system.Jt[i]^2) / system.v[i[2]]
        else
            row[idxi] = idxr
            col[idxi] = i[1]

            coeff[idxi] = system.Jt[i]
            coeffInv[idxi] = 1 / system.Jt[i]

            Mind[idxi] = system.z[i[2]]
            Vind[idxi] = system.v[i[2]]

            idxi += 1

            colptr[i[1]] += 1
            rowptr[idxr] = idxi

            if idxT[system.Jt.colptr[i[2] + 1] - 1] == i
                idxr += 1
            end
        end
        if virtual[i[1]] !== 0
            Mdir[i[1]] = settings.mean / settings.variance
            Wdir[i[1]] = 1 / settings.variance
        end
    end
    pushfirst!(colptr, 1)
    pushfirst!(rowptr, 1)
    colptr = cumsum(colptr)

    return FactorGraph(row, rowptr, col, colptr, Mind, Vind, coeff, coeffInv, Mdir, Wdir, Nvar, Nfac, Ndir, Nind, Nlink)
end

### Initialize the algorithm
function initialize(graph, settings)
    ########## Initialize arrays for messages ##########
    Mfac_var = similar(graph.coeff)
    Wfac_var = similar(graph.coeff)
    Mvar_fac = similar(graph.coeff)
    Vvar_fac = similar(graph.coeff)

    ########## Set damping parameters ##########
    bernoulli_sample = randsubseq(collect(1:graph.Nlink), settings.prob)
    alpha1 = fill(1.0, graph.Nlink)
    alpha2 = fill(0.0, graph.Nlink)

    @inbounds for i in bernoulli_sample
        alpha1[i] = 1.0 - settings.alpha
        alpha2[i] = settings.alpha
    end

    ########## Indices for sending messages ##########
    temp = collect(1:graph.Nlink)
    sort_col_idx = sortperm(graph.col)
    to_fac = temp[sort_col_idx]

    new_row = graph.row[sort_col_idx]
    sort_row_idx = sortperm(new_row)
    to_var = temp[sort_row_idx]

    ########## Pass messages from singly-connected factor nodes to all indirect links ##########
    @inbounds for i = 1:graph.Nvar
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            k = to_fac[j]
            Vvar_fac[k] = 1 / graph.Wdir[i]
            Mvar_fac[k] = graph.Mdir[i] * Vvar_fac[k]
        end
    end

    ########## Iteration counter ##########
    iterCount = [0]
    
    return BeliefPropagation(Mfac_var, Wfac_var, Mvar_fac, Vvar_fac, alpha1, alpha2, to_fac, to_var, iterCount)
end

### Initialize results data
function resultsdata(system, settings, Nvar)
    if settings.iterate
        iterations = collect(1:settings.maxIter)
        means = zeros(Nvar, settings.maxIter)
        variances = zeros(Nvar, settings.maxIter)
        if settings.loss
            rmse = fill(0.0, settings.maxIter)
            mae = fill(0.0, settings.maxIter)
            wrss = fill(0.0, settings.maxIter)
        end
    else
        means = fill(0.0, Nvar)
        variances = fill(0.0, Nvar)
        iterations = settings.maxIter
        if settings.loss
            rmse = [0.0]
            mae = [0.0]
            wrss = [0.0]
        end
    end

    if settings.wls
        N, ~ = size(system.J)
        W = spdiagm(0 =>  @. 1.0 / sqrt(system.v))
        H = W * system.J
        G = H' * H
        b = H' * W * system.z
        meansWLS = G \ b

        z_wls = system.J * meansWLS
        rmseWLS = [((sum(system.z - z_wls).^2) / N)^(1/2)]
        maeWLS = [(sum(abs.(system.z - z_wls))) / N]
        wrssWLS = [sum(((system.z - z_wls).^2) ./ system.v)]
        if settings.iterate
            rmseBPWLS = fill(0.0, settings.maxIter)
            maeBPWLS = fill(0.0, settings.maxIter)
        else
            rmseBPWLS = [0.0]
            maeBPWLS = [0.0]
        end
    end

    return if settings.loss && !settings.wls
                ResultsError(iterations, means, variances, rmse, mae, wrss) 
            elseif settings.loss && settings.wls
                ResultsErrorWLS(iterations, means, variances, rmse, mae, wrss, meansWLS, rmseWLS, maeWLS, wrssWLS, rmseBPWLS, maeBPWLS) 
            elseif !settings.loss && settings.wls
                ResultsWLS(iterations, means, variances, meansWLS, rmseWLS, maeWLS, wrssWLS, rmseBPWLS, maeBPWLS) 
            else
                Results(iterations, means, variances) 
            end 
end

### Initialize arrays for recursion
function recursionini(graph)
    Mcol = fill(0.0, graph.Nvar)
    Wcol = fill(0.0, graph.Nvar)

    return Recursion(Mcol, Wcol)
end

### Error metrics
function stats(system, results, settings)
    N, Nv = size(system.J)
    if settings.loss 
        z_bp = system.J * results.means
        if settings.iterate 
            for i = 1:settings.maxIter
                results.rmse[i] = ((sum(system.z - z_bp[:, i]).^2) / N)^(1/2)
                results.mae[i] = (sum(abs.(system.z - z_bp[:, i]))) / N
                results.wrss[i] = sum(((system.z - z_bp[:, i]).^2) ./ system.v)
            end
        else
            results.rmse[1] = ((sum(system.z - z_bp).^2) / N)^(1/2)
            results.mae[1] = (sum(abs.(system.z - z_bp))) / N
            results.wrss[1] = sum(((system.z - z_bp).^2) ./ system.v)
        end
    end

    if settings.wls
        if settings.iterate
            for i = 1:settings.maxIter
            results.rmseBPWLS[i] = ((sum(results.meansWLS - results.means[:, i]).^2) / Nv)^(1/2)
            results.maeBPWLS[i] = (sum(abs.(results.meansWLS - results.means[:, i]))) / Nv
            end 
        else
            results.rmseBPWLS[1] = ((sum(results.meansWLS - results.means).^2) / Nv)^(1/2)
            results.maeBPWLS[1] = (sum(abs.(results.meansWLS - results.means))) / Nv
        end
    end
end

### Display
function displaydata(graph, bp, results, settings, prep, infe)
    ########## Factor graph data and iterate scheme ##########
    data = ["Factor graph data and Gaussian belief propagation iterations scheme" "" ""
            "Number of Factor Nodes" "."^10 graph.Nfac;
            "Number of Variable Nodes" "."^10 graph.Nvar;
            "Number of Direct Links" "."^10 graph.Ndir;
            "Number of Indirect Links" "."^10 graph.Nlink;
            "" "" "";
            "Number of iterations" "."^10 settings.maxIter;
            "Number of iterations where means and variances are computed" "."^10 settings.IterNative;
            "Number of iterations where means and variances are computed using damping" "."^10 settings.IterDamp;
            "Number of iterations where only means are computed" "."^10 settings.IterBump;
            "Number of iterations where only means are computed using damping" "."^10 settings.IterDampBump]

        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 12],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.0f", 3),
        highlighters = (hl_cell( [(1,1)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [1], body_hlines_format = Tuple('─' for _ = 1:4))    
     
    ########## BP running time ##########
    data = ["" "" ""
            "Gaussian belief propagation time in milliseconds" "" "";
            "Preprocesing" "."^10 1000 * prep;
            "Inference" "."^10 1000 * infe;
            "Total" "."^10 1000 * infe + 1000 * prep]
    pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 12],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4f", 3),
        highlighters = (hl_cell( [(2,1)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))
   
    ########## BP and WLS error metrics ##########
    if settings.loss && !settings.wls
        data = ["" "" ""
        "Error Metrics" "" ""
        "Root mean square error of the Gaussian belief propagation" "."^10 results.rmse[end];
        "Mean absolute error of the Gaussian belief propagation" "."^10 results.mae[end];
        "Weighted residual sum of squares error of the Gaussian belief propagation" "."^10 results.wrss[end]]
        
        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 12],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4e", 3),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))
    end
    if settings.wls
        data = ["" "" ""
                "Error Metrics" "" ""
                "Root mean square error of the Gaussian belief propagation" "."^10 results.rmse[end];
                "Root mean square error of the weighted least-squares method" "."^10 results.rmseWLS
                "Root mean square error ratio" "."^10 results.rmse[end] / results.rmseWLS
                "" "" "";
                "Mean absolute error of the Gaussian belief propagation" "."^10 results.mae[end];
                "Mean absolute error of the weighted least-squares method" "."^10 results.maeWLS
                "Mean absolute error ratio" "."^10 results.mae[end] / results.maeWLS
                "" "" "";
                "Weighted residual sum of squares of the belief propagation" "."^10 results.wrss[end];
                "Weighted residual sum of squares of the weighted least-squares method" "."^10 results.wrssWLS;
                "Weighted residual sum of squares ratio" "."^10 results.wrss[end] / results.wrssWLS]

        high = Highlighter((data,i,j) -> i == 5 && j == 3 && data[5,3] > 1.1, Crayon(bold = true, background = :red))
        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 12],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4e", 3),
        highlighters = (high, hl_cell([(2,1), (8,3), (5,3)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))
    end

    if settings.wls
        ########## BP and WLS state variables ##########
        A = [collect(1:graph.Nvar) results.means[:, end] results.meansWLS results.means[:, end] ./ results.meansWLS]
        A = A[reverse(sortperm(A[:, 4])),  :]

        println("\n State Variables Display")
        pretty_table(A, header = ["State Variable", "Belief Propagation Estimate", "Weighted Least-squares Estimate", "Maximum to Minimum Ratio"],
            alignment=[:r, :r, :r, :r], formatters = ft_printf(["%3.0f", "%3.6f","%3.6f","%3.6e"], [1, 2, 3, 4]))
    else
        ########## BP state variables ##########
        println("\n State Variables Display")
        pretty_table([collect(1:graph.Nvar) results.means[:, end]], header = ["State Variable", "Estimate"],
            alignment=[:r, :r], formatters = ft_printf(["%3.0f", "%3.6f"], [1, 2]))
    end
end
