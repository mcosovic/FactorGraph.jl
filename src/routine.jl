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
    IterNative::Int64
    IterDamp::Int64
    IterBump::Int64
    IterDampBump::Int64
end

struct Recursion
    Mcol::Array{Float64,1}
    Wcol::Array{Float64,1}
end


### Produce the factor graph and define singly-connected, virtual and indirect factor nodes
function factors(system, meanVirtual, variVirtual)
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
            Mdir[i[1]] = meanVirtual / variVirtual
            Wdir[i[1]] = 1 / variVirtual
        end
    end
    pushfirst!(colptr, 1)
    pushfirst!(rowptr, 1)
    colptr = cumsum(colptr)

    return FactorGraph(row, rowptr, col, colptr, Mind, Vind, coeff, coeffInv, Mdir, Wdir, Nvar, Nfac, Ndir, Nind, Nlink)
end


### Initialize the algorithm
function initialize(graph, maxIter, damp, bump, prob, alpha)
    ########## Initialize arrays for messages ##########
    Mfac_var = similar(graph.coeff)
    Wfac_var = similar(graph.coeff)
    Mvar_fac = similar(graph.coeff)
    Vvar_fac = similar(graph.coeff)

    ########## Set damping parameters ##########
    bernoulli_sample = randsubseq(collect(1:graph.Nlink), prob)
    alpha1 = fill(1.0, graph.Nlink)
    alpha2 = fill(0.0, graph.Nlink)

    @inbounds for i in bernoulli_sample
        alpha1[i] = 1.0 - alpha
        alpha2[i] = alpha
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

    ########## Iteration schemes ##########
    dabu = minimum([bump, damp])
    IterNative = 0
    if bump == damp == maxIter
        IterNative = maxIter
    else
        IterNative = dabu - 1
    end

    IterDamp = 0
    if bump == maxIter
        IterDamp = maxIter - IterNative
    elseif bump > damp
        IterDamp = bump - damp
    end

    IterBump = 0
    if damp == maxIter
        IterBump = maxIter - IterNative
    elseif bump < damp
        IterBump = damp - bump
    end

    IterDampBump = 0
    if damp < maxIter && bump < maxIter
       IterDampBump = maxIter - IterDamp - IterBump - IterNative
    end

    return BeliefPropagation(Mfac_var, Wfac_var, Mvar_fac, Vvar_fac, alpha1, alpha2, to_fac, to_var, IterNative, IterDamp, IterBump, IterDampBump)
end

########## Initialize arrays for recursion ##########
function recursionini(graph)
    Mcol = fill(0.0, graph.Nvar)
    Wcol = fill(0.0, graph.Nvar)

    return Recursion(Mcol, Wcol)
end


### Results
function results(system, graph, bp, Xbp, wls, prep, infe)
    ########## Graph data and iterate scheme ##########

    data = ["Factor graph data and belief propagation iterations scheme" "" ""
            "Number of Factor Nodes" "."^10 graph.Nfac;
            "Number of Variable Nodes" "."^10 graph.Nvar;
            "Number of Direct Links" "."^10 graph.Ndir;
            "Number of Indirect Links" "."^10 graph.Nlink;
            "Number of iterations where means and variances are computed" "."^10 bp.IterNative;
            "Number of iterations where means and variances are computed using damping" "."^10 bp.IterDamp;
            "Number of iterations where only means are computed" "."^10 bp.IterBump;
            "Number of iterations where only means are computed using damping" "."^10 bp.IterDampBump]

        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 12],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.0f", 3),
        highlighters = (hl_cell( [(1,1)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [1], body_hlines_format = Tuple('─' for _ = 1:4))    

            
    ########## BP running time ##########
    data = ["" "" ""
            "Belief Propagation Time in Milliseconds" "" "";
            "Preprocesing" "."^10 1000 * prep;
            "Inference" "."^10 1000 * infe;
            "Total" "."^10 1000 * infe + 1000 * prep]
    pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 12],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4f", 3),
        highlighters = (hl_cell( [(2,1)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))

    ########## BP and WLS error metrics ##########
    N, ~ = size(system.J)
    z_bp = system.J * Xbp
    wrss_bp = sum(((system.z - z_bp).^2) ./ system.v)
    rmse_bp = ((sum(system.z - z_bp).^2) / N)^(1/2)

    if wls == "builtin" || wls == "lu"
        W = spdiagm(0 =>  @. 1.0 / sqrt(system.v))
        H = W * system.J
        G = H' * H
        b = H' * W * system.z
        if wls == "builtin"
            Xwls = G \ b
        end
        if wls == "lu"
            F = lu(G)
            L, U, p, q, Rs = F.:(:)
            Xwls = U \  (L \ ((Rs .* b)[p]))
            Xwls = Xwls[sortperm(q)]
        end
        z_wls = system.J * Xwls
        wrss_wls = sum(((system.z - z_wls).^2) ./ system.v)
        rmse_wls = ((sum(system.z - z_wls).^2) / N)^(1/2)

        A = [collect(1:length(Xbp)) Xbp Xwls Xbp ./ Xwls]
        A = A[reverse(sortperm(A[:, 4])),  :]

        data = ["" "" ""
                "Error Metrics" "" ""
                "Weighted residual sum of squares of the belief propagation" "."^10 wrss_bp;
                "Weighted residual sum of squares of the weighted least-squares method" "."^10 wrss_wls;
                "Weighted residual sum of squares ratio" "."^10 wrss_bp / wrss_wls;
                "Root mean square of the belief propagation" "."^10 rmse_bp;
                "Root mean square of the weighted least-squares method" "."^10 rmse_wls
                "Root mean square ratio" "."^10 rmse_bp / rmse_wls]
        high = Highlighter((data,i,j) -> i == 5 && j == 3 && data[5,3] > 1.1, Crayon(bold = true, background = :red))
    else
        data = ["" "" ""
                "Error Metrics" "" ""
                "Weighted residual sum of squares of the belief propagation" "."^10 wrss_bp;
                "Root mean square of the belief propagation" "."^10 rmse_bp]
        high = Highlighter(hl_col(2, crayon"dark_gray"))

    end
    pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 12],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4e", 3),
        highlighters = (high, hl_cell([(2,1), (8,3), (5,3)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))

    if wls == "builtin" || wls == "lu"
        println("\n State Variables Display")
        pretty_table(A, header = ["State Variable", "Belief Propagation Estimate", "Weighted Least-squares Estimate", "Maximum to Minimum Ratio"],
            alignment=[:r, :r, :r, :r], formatters = ft_printf(["%3.0f", "%3.6f","%3.6f","%3.6e"], [1, 2, 3, 4]))
    else
        println("\n State Variables Display")
        pretty_table([collect(1:length(Xbp)) Xbp], header = ["State Variable", "Belief Propagation Estimate"],
            alignment=[:r, :r], formatters = ft_printf(["%3.0f", "%3.6f"], [1, 2]))
    end
end
