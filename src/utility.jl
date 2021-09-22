struct WeightedLeastSquares
    estimate::Array{Float64,1}
    rmse::Float64
    mae::Float64
    wrss::Float64
end

struct ErrorMetric
    rmse::Float64
    mae::Float64
    wrss::Float64
end

struct ErrorMetricWiden
    rmse::Float64
    mae::Float64
    wrss::Float64
    rmseGBPWLS::Float64
    maeGBPWLS::Float64
end

########## Compute the GBP error metrics ##########
function errorMetric(gbp::Union{ContinuousModel, ContinuousTreeModel})
    observationGBP = gbp.system.jacobian * gbp.inference.mean
    NactiveRows = activeRows(gbp)

    rmse = 0.0; mae = 0.0; wrss = 0.0
    @inbounds for i = 1:length(gbp.system.observation)
        rmse += (gbp.system.observation[i] - observationGBP[i])^2
        mae += abs(gbp.system.observation[i] - observationGBP[i])
        wrss += (gbp.system.observation[i] - observationGBP[i]) / gbp.system.variance[i]
    end
    rmse = (rmse / NactiveRows)^(1/2)
    mae = mae / NactiveRows

    return ErrorMetric(rmse, mae, wrss)
end

########## Compute the GBP error metrics and error according to the WLS ##########
@inline function errorMetric(gbp::Union{ContinuousModel, ContinuousTreeModel}, wls::WeightedLeastSquares)
    observationGBP = gbp.system.jacobian * gbp.inference.mean
    NactiveRows = activeRows(gbp)

    rmse = 0.0; mae = 0.0; wrss = 0.0
    @inbounds for i = 1:length(gbp.system.observation)
        rmse += (gbp.system.observation[i] - observationGBP[i])^2
        mae += abs(gbp.system.observation[i] - observationGBP[i])
        wrss += (gbp.system.observation[i] - observationGBP[i]) / gbp.system.variance[i]
    end
    rmse = (rmse / NactiveRows)^(1/2)
    mae = mae / NactiveRows

    rmseGBPWLS = 0.0; maeGBPWLS = 0.0
    @inbounds for i = 1:gbp.graph.Nvariable
        rmseGBPWLS += (wls.estimate[i] - gbp.inference.mean[i])^2
        maeGBPWLS += abs(wls.estimate[i] - gbp.inference.mean[i])
    end
    rmseGBPWLS = (rmseGBPWLS / gbp.graph.Nvariable)^(1/2)
    maeGBPWLS = maeGBPWLS / gbp.graph.Nvariable

    return ErrorMetricWiden(rmse, mae, wrss, rmseGBPWLS, maeGBPWLS)
end

########## Compute WLS solution and error metrics ##########
function wls(gbp::Union{ContinuousModel, ContinuousTreeModel})
    W = spdiagm(0 =>  @. 1.0 / sqrt(gbp.system.variance))
    A = W * gbp.system.jacobian
    G = A' * A
    b = A' * W * gbp.system.observation
    x = G \ b

    observationWLS = gbp.system.jacobian * x
    NactiveRows = activeRows(gbp)

    rmse = 0.0; mae = 0.0; wrss = 0.0
    @inbounds for i = 1:length(gbp.system.observation)
        rmse += (gbp.system.observation[i] - observationWLS[i])^2
        mae += abs(gbp.system.observation[i] - observationWLS[i])
        wrss += (gbp.system.observation[i] - observationWLS[i]) / gbp.system.variance[i]
    end
    rmse = (rmse / NactiveRows)^(1/2)
    mae = mae / NactiveRows

    return WeightedLeastSquares(x, rmse, mae, wrss)
end

########## Display: GBP and WLS ##########
function displayData(args...)
    gbpIdx = 0; wlsIdx = 0; errorIdx = 0; widenIdx = 0
    @inbounds for (k, i) in enumerate(args)
        T = split.(string(typeof(i)), ".")[end]
        if T == "ContinuousModel" || T == "ContinuousTreeModel"
            gbpIdx = k
        end
        if T == "WeightedLeastSquares"
            wlsIdx = k
        end
        if T == "ErrorMetric" || T == "ErrorMetricWiden"
            errorIdx = k
        end
    end
    if gbpIdx == 0
       error("FactorGraph model is missing.")
    end

    if errorIdx != 0 && wlsIdx == 0
        data = ["" "" ""
                "Error Metrics" "" ""
                "Root mean square error of the Gaussian belief propagation" "."^10 getfield(args[errorIdx], :rmse);
                "Mean absolute error of the Gaussian belief propagation" "."^10 getfield(args[errorIdx], :mae);
                "Weighted residual sum of squares error of the Gaussian belief propagation" "."^10 getfield(args[errorIdx], :wrss)]

        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 14],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4e", 3),
        highlighters = (hl_cell( [(2,1)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))
    end

    if errorIdx != 0 && wlsIdx != 0
        data = ["" "" ""
                "Error Metrics" "" ""
                "Root mean square error of the Gaussian belief propagation" "."^10 getfield(args[errorIdx], :rmse);
                "Root mean square error of the weighted least-squares method" "."^10  getfield(args[wlsIdx], :rmse);
                "Root mean square error ratio" "."^10 getfield(args[errorIdx], :rmse) / getfield(args[wlsIdx], :rmse);
                "" "" "";
                "Mean absolute error of the Gaussian belief propagation" "."^10 getfield(args[errorIdx], :mae);
                "Mean absolute error of the weighted least-squares method" "."^10 getfield(args[wlsIdx], :mae);
                "Mean absolute error ratio" "."^10 getfield(args[errorIdx], :mae) / getfield(args[wlsIdx], :mae);
                "" "" "";
                "Weighted residual sum of squares of the belief propagation" "."^10 getfield(args[errorIdx], :wrss);
                "Weighted residual sum of squares of the weighted least-squares method" "."^10 getfield(args[wlsIdx], :wrss);
                "Weighted residual sum of squares ratio" "."^10 getfield(args[errorIdx], :wrss) / getfield(args[wlsIdx], :wrss)]

        high = Highlighter((data,i,j) -> i == 5 && j == 3 && data[5,3] > 1.1, Crayon(bold = true, background = :red))
        pretty_table(data, tf = tf_borderless, columns_width = [75, 11, 14],
        noheader = true, alignment = [:l, :l, :r], formatters = ft_printf("%1.4e", 3),
        highlighters = (high, hl_cell([(2,1), (8,3), (5,3)], crayon"bold"), hl_col(2, crayon"dark_gray")),
        body_hlines = [2], body_hlines_format = Tuple('─' for _ = 1:4))
    end

    if wlsIdx == 0
        println("\n State Variables Display")
        pretty_table([collect(1:getfield(args[gbpIdx].graph, :Nvariable)) getfield(args[gbpIdx].inference, :mean)], header = ["State Variable", "Estimate"],
        alignment=[:r, :r], formatters = ft_printf(["%3.0f", "%3.6f"], [1, 2]))
    end

    if wlsIdx != 0
        A = [collect(1:getfield(args[gbpIdx].graph, :Nvariable)) getfield(args[gbpIdx].inference, :mean) getfield(args[wlsIdx], :estimate) getfield(args[gbpIdx].inference, :mean) ./ getfield(args[wlsIdx], :estimate)]
        A = A[reverse(sortperm(A[:, 4])),  :]

        println("\n State Variables Display")
        pretty_table(A, header = ["State Variable", "Belief Propagation Estimate", "Weighted Least-squares Estimate", "Estimates Ratio"],
        columns_width = [15, 28, 32, 20], alignment=[:r, :r, :r, :r], formatters = ft_printf(["%3.0f", "%3.6f", "%3.6f", "%3.6e"], [1, 2, 3, 4]))
    end
end

########## Find number of nonzeros rows ##########
@inline function activeRows(gbp)
    NactiveRows = length(gbp.system.observation)
    @inbounds for (k, i) in enumerate(gbp.system.observation)
        if i == 0.0
            if all(gbp.system.jacobianTranspose[:, k] .== 0.0)
                NactiveRows -= 1
            end
        end
    end

    return NactiveRows
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

########## Check data for import ##########
function checkImportFile(args)
    #### Check the package is installed
    pathtoFactorGraph = Base.find_package("FactorGraph")
    if isnothing(pathtoFactorGraph)
        throw(ErrorException("FactorGraph not found in install packages"))
    end
    packagepath = abspath(joinpath(dirname(pathtoFactorGraph), ".."))

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

    return fullpath, extension, dataname
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

    #### Check initial variance
    if variance < 0.0
        error("Invalid variance value.")
    end
end

########## Type of the system data ##########
function checkFileOrArguments(args)
    fromfile = false
    @inbounds for i in args
        if typeof(i) == String
            fromfile = true
            break
        end
    end

    return fromfile
end