struct ResultsGBP
    iteration::Union{Int64, Array{Int64,1}}
    mean::Union{Array{Float64,2}, Array{Float64,1}} 
    variance::Union{Array{Float64,2}, Array{Float64,1}} 
    rmse::Union{Array{Float64,1}, Nothing} 
    mae::Union{Array{Float64,1}, Nothing} 
    wrss::Union{Array{Float64,1}, Nothing} 
end

struct ResultsWLS
    mean::Union{Array{Float64,2}, Array{Float64,1}, Nothing}
    rmse::Union{Array{Float64,1}, Nothing}
    mae::Union{Array{Float64,1}, Nothing}
    wrss::Union{Array{Float64,1}, Nothing}
    rmseGBPWLS::Union{Array{Float64,1}, Nothing}
    maeGBPWLS::Union{Array{Float64,1}, Nothing} 
end

struct Results
    gbp::ResultsGBP
    wls::ResultsWLS
end

########## Initialize results data ##########
function initialize_results(settings, system, graph)
    #### Dynamic data
    if settings.dynamic
        dynamic = trunc.(Int, unique(system.dynamic[:, 1]))
        whichDyn = dynamic .<= settings.maxIter
        dynamic = dynamic[whichDyn]
        Ndynamic = length(dynamic) + 1
    end

    #### Initialize the GBP results
    if settings.outIterate
        iteration = collect(1:settings.maxIter)
        mean = zeros(graph.Nvar, settings.maxIter)
        variance = zeros(graph.Nvar, settings.maxIter)
    elseif settings.dynamic
        iteration = [dynamic .- 1; settings.maxIter]
        mean = zeros(graph.Nvar, Ndynamic)
        variance = zeros(graph.Nvar, Ndynamic)
    else
        iteration = settings.maxIter
        mean = fill(0.0, graph.Nvar)
        variance = fill(0.0, graph.Nvar)
    end

    #### Initialize the GBP error metrics
    if settings.outEvaluation
        if settings.outIterate
            rmse = fill(0.0, settings.maxIter)
            mae = fill(0.0, settings.maxIter)
            wrss = fill(0.0, settings.maxIter)
        elseif settings.dynamic
            rmse = fill(0.0, Ndynamic)
            mae = fill(0.0, Ndynamic)
            wrss = fill(0.0, Ndynamic)
        else
            rmse = [0.0]
            mae = [0.0]
            wrss = [0.0]
        end
    else
        rmse = nothing
        mae = nothing
        wrss = nothing
    end

    if settings.outWls
        if settings.outIterate
            rmseGBPWLS = fill(0.0, settings.maxIter)
            maeGBPWLS = fill(0.0, settings.maxIter)
        elseif settings.dynamic
            rmseGBPWLS = fill(0.0, Ndynamic)
            maeGBPWLS = fill(0.0, Ndynamic)
        else
            rmseGBPWLS = [0.0]
            maeGBPWLS = [0.0]
        end

        if !settings.dynamic
            meanWLS = fill(0.0, graph.Nvar)
            rmseWLS = [0.0]      
            maeWLS = [0.0]
            wrssWLS = [0.0]
        else 
            meanWLS = zeros(graph.Nvar, Ndynamic)
            rmseWLS = fill(0.0, Ndynamic)       
            maeWLS = fill(0.0, Ndynamic)
            wrssWLS = fill(0.0, Ndynamic)
        end
    else
        meanWLS = nothing
        rmseWLS = nothing
        maeWLS = nothing
        wrssWLS = nothing
        rmseGBPWLS = nothing
        maeGBPWLS = nothing 
    end

    return Results(ResultsGBP(iteration, mean, variance, rmse, mae, wrss), ResultsWLS(meanWLS, rmseWLS, maeWLS, wrssWLS, rmseGBPWLS, maeGBPWLS))    
end

########## Compute the GBP marginals and error metrics ##########
@inline function marginal(settings, system, graph, bp, results)
    if bp.iterCount[1] == results.gbp.iteration[bp.updateCount[1]]
        k = bp.updateCount[1]
        @inbounds for i = 1:graph.Nvar
            Mcol = 0.0; Wcol = 0.0
    
            for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
                Mcol += bp.Mfac_var[j] * bp.Wfac_var[j]
                Wcol += bp.Wfac_var[j]
            end
            results.gbp.variance[i, k] = 1 / (Wcol + graph.Wdir[i])
            results.gbp.mean[i, k] = (Mcol + graph.Mdir[i]) * results.gbp.variance[i, k]        
        end 
        
        if settings.outEvaluation
            observationGBP = system.jacobian * results.gbp.mean[:, k]
            
            if settings.dynamic
                z = graph.observationDyn
                v = graph.varianceDyn
            else
                z = system.observation
                v = system.variance
            end

            @inbounds for i = 1:graph.Nfac
                results.gbp.rmse[k] += (z[i] - observationGBP[i])^2 
                results.gbp.mae[k] += abs(z[i] - observationGBP[i]) 
                results.gbp.wrss[k] += (z[i] - observationGBP[i]) / v[i]
            end
            results.gbp.rmse[k] = (results.gbp.rmse[k] / graph.Nfac)^(1/2)
            results.gbp.mae[k] = results.gbp.mae[k] / graph.Nfac
        end 
        bp.updateCount[1] += 1    
    end 
end

########## Compute WLS solution and error metrics ##########
function wls_solution(H, z, v, graph, results, idx)
    W = spdiagm(0 =>  @. 1.0 / sqrt(v))
    A = W * H
    G = A' * A
    b = A' * W * z
    results.wls.mean[:, idx] = G \ b

    observationWLS = H * results.wls.mean[:, idx]

    @inbounds for i = 1:graph.Nfac
        results.wls.rmse[idx] += (z[i] - observationWLS[i])^2 
        results.wls.mae[idx] += abs(z[i] - observationWLS[i]) 
        results.wls.wrss[idx] += (z[i] - observationWLS[i]) / v[i]
    end
    results.wls.rmse[idx] = (results.wls.rmse[idx] / graph.Nfac)^(1/2)
    results.wls.mae[idx] = results.wls.mae[idx] / graph.Nfac
end

########## Compute WLS solution and error metrics depending on modules ##########
function wls_metrics(system, settings, graph, results)
    if settings.dynamic
        #### Initialize observation and variance values
        graph.observationDyn[:] = system.observation
        graph.varianceDyn[:] = system.variance

        #### Save WLS solution before update occurs
        dynamic = trunc.(Int, unique(system.dynamic[:, 1]))
        dynamic = dynamic[dynamic .<= settings.maxIter]
        saveWLS = [dynamic .- 1; settings.maxIter]

        flagDyn = 1; flagSave = 1
        
        @inbounds for i = 1:settings.maxIter
            while graph.Ndyn >= flagDyn && system.dynamic[flagDyn, 1] == i 
                factor = trunc(Int, system.dynamic[flagDyn, 2])

                graph.observationDyn[factor] = system.dynamic[flagDyn, 3] 
                graph.varianceDyn[factor] = system.dynamic[flagDyn, 4] 

                flagDyn += 1
            end
            
            if saveWLS[flagSave] == i
                wls_solution(system.jacobian, graph.observationDyn, graph.varianceDyn, graph, results, flagSave)

                if !settings.outIterate
                    results.wls.rmseGBPWLS[flagSave] = ((sum(results.wls.mean[:, flagSave] - results.gbp.mean[:, flagSave]).^2) / graph.Nvar)^(1/2)
                    results.wls.maeGBPWLS[flagSave] = (sum(abs.(results.wls.mean[:, flagSave] - results.gbp.mean[:, flagSave]))) / graph.Nvar
                end
            end

            if settings.outIterate
                results.wls.rmseGBPWLS[i] = ((sum(results.wls.mean[:, flagSave] - results.gbp.mean[:, i]).^2) / graph.Nvar)^(1/2)
                results.wls.maeGBPWLS[i] = (sum(abs.(results.wls.mean[:, flagSave] - results.gbp.mean[:, i]))) / graph.Nvar
            end

            if saveWLS[flagSave] == i
                flagSave += 1
            end
        end
    else
        wls_solution(system.jacobian, system.observation, system.variance, graph, results, 1)

        if settings.outIterate
            @inbounds for i = 1:settings.maxIter
                results.wls.rmseGBPWLS[i] = ((sum(results.wls.mean - results.gbp.mean[:, i]).^2) / graph.Nvar)^(1/2)
                results.wls.maeGBPWLS[i] = (sum(abs.(results.wls.mean - results.gbp.mean[:, i]))) / graph.Nvar
            end
        else
            results.wls.rmseGBPWLS[1] = ((sum(results.wls.mean - results.gbp.mean).^2) / graph.Nvar)^(1/2)
            results.wls.maeGBPWLS[1] = (sum(abs.(results.wls.mean - results.gbp.mean))) / graph.Nvar
        end
    end
end

########## Display: Factor graph data and iterate scheme ##########
function display_stat(graph, bp, settings, system, algorithm)
    if settings.outDisplay
        data = ["The Gaussian belief propagation and factor graph data" "" "";
                "Input data" "."^10 system.data;
                "Algorithm type" "."^10 algorithm;
                "" "" ""]

        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 14],
            noheader = true, alignment = [:l, :l, :r],
            highlighters = (hl_cell( [(1,1)], crayon"bold"), hl_col(2, crayon"dark_gray")),
            body_hlines = [1], body_hlines_format = Tuple('─' for _ = 1:4))  

        data = ["Number of factor nodes" "."^10 graph.Nfac;
                "Number of variable nodes" "."^10 graph.Nvar;
                "Number of direct links" "."^10 graph.Ndir;
                "Number of indirect links" "."^10 graph.Nlink;
                "" "" "";
                "Number of iterations" "."^10 settings.maxIter;
                "Number of iterations where means and variances are computed" "."^10 bp.iterNative;
                "Number of iterations where means and variances are computed using damping" "."^10 bp.iterDamp;
                "Number of iterations where only means are computed" "."^10 bp.iterBump;
                "Number of iterations where only means are computed using damping" "."^10 bp.iterDampBump]

        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 14],
            noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.0f", 3),
            highlighters = (hl_col(2, crayon"dark_gray")),
            body_hlines_format = Tuple('─' for _ = 1:4))  
    end  
end

########## Display ##########
function display_data(graph, bp, results, settings)
    #### GBP and WLS error metrics
    if settings.outEvaluation && !settings.outWls
        data = ["" "" ""
                "Error Metrics" "" ""
                "Root mean square error of the Gaussian belief propagation" "."^10 results.gbp.rmse[end];
                "Mean absolute error of the Gaussian belief propagation" "."^10 results.gbp.mae[end];
                "Weighted residual sum of squares error of the Gaussian belief propagation" "."^10 results.gbp.wrss[end]]
        
        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 14],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4e", 3),
        highlighters = (hl_cell( [(2,1)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))
    end
    if settings.outEvaluation && settings.outWls
        data = ["" "" ""
                "Error Metrics" "" ""
                "Root mean square error of the Gaussian belief propagation" "."^10 results.gbp.rmse[end];
                "Root mean square error of the weighted least-squares method" "."^10 results.wls.rmse[end];
                "Root mean square error ratio" "."^10 results.gbp.rmse[end] / results.wls.rmse[end];
                "" "" "";
                "Mean absolute error of the Gaussian belief propagation" "."^10 results.gbp.mae[end];
                "Mean absolute error of the weighted least-squares method" "."^10 results.wls.mae[end];
                "Mean absolute error ratio" "."^10 results.gbp.mae[end] / results.wls.mae[end];
                "" "" "";
                "Weighted residual sum of squares of the belief propagation" "."^10 results.gbp.wrss[end];
                "Weighted residual sum of squares of the weighted least-squares method" "."^10 results.wls.wrss[end];
                "Weighted residual sum of squares ratio" "."^10 results.gbp.wrss[end] / results.wls.wrss[end]]

        high = Highlighter((data,i,j) -> i == 5 && j == 3 && data[5,3] > 1.1, Crayon(bold = true, background = :red))
        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 14],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4e", 3),
        highlighters = (high, hl_cell([(2,1), (8,3), (5,3)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))
    end

    if settings.outWls
        #### GBP and WLS state variables
        A = [collect(1:graph.Nvar) results.gbp.mean[:, end] results.wls.mean[:, end] results.gbp.mean[:, end] ./ results.wls.mean[:, end]]
        A = A[reverse(sortperm(A[:, 4])),  :]

        println("\n State Variables Display")
        pretty_table(A, header = ["State Variable", "Belief Propagation Estimate", "Weighted Least-squares Estimate", "Estimates Ratio"],
        columns_width = [15, 28, 32, 20], alignment=[:r, :r, :r, :r], formatters = ft_printf(["%3.0f", "%3.6f", "%3.6f", "%3.6e"], [1, 2, 3, 4]))
    else
        #### GBP state variables
        println("\n State Variables Display")
        pretty_table([collect(1:graph.Nvar) results.gbp.mean[:, end]], header = ["State Variable", "Estimate"],
            alignment=[:r, :r], formatters = ft_printf(["%3.0f", "%3.6f"], [1, 2]))
    end
end