struct SystemModel
    J::SparseMatrixCSC{Float64,Int64}
    Jt::SparseMatrixCSC{Float64,Int64}
    z::Array{Float64,1}
    v::Array{Float64,1}
    data::String
end

struct Settings
    maxIter::Int64
    IterNative::Int64
    IterDamp::Int64
    IterBump::Int64
    IterDampBump::Int64
    prob::Float64
    alpha::Float64
    mean::Float64
    variance::Float64
    algorithm::String
    iterate::Bool
    loss::Bool
    wls::Bool
    displayshow::Bool
end


### Type of the system data
function checkInputs(args)
    fromfile = false
    for i in args
        if typeof(i) == String
            fromfile = true
            break
        end
    end

    return fromfile
end

### Load from path
function outJulia(args)
    pathtoGaussBP = Base.find_package("GaussBP")
    if pathtoGaussBP .== nothing
        throw(ErrorException("GaussBP not found in install packages"))
    end
    packagepath = abspath(joinpath(dirname(pathtoGaussBP), ".."))
    extension = ".h5"; path = ""; dataname = ""; fullpath = ""

    for i = 1:length(args)
        try
            extension = string(match(r"\.[A-Za-z0-9]+$", args[i]).match)
        catch
            extension = ""
        end
        if extension == ".h5" || extension == ".csv" 
            fullpath = args[i]
            path = dirname(args[i])
            dataname = basename(args[i])
            break
        end
    end

    if isempty(extension)
        throw(ErrorException("the input DATA extension is not found"))
    elseif extension != ".h5" && extension != ".csv"
        throw(DomainError(extension, "the input DATA extension is not supported"))
    end

    if path == ""
        path = joinpath(packagepath, "src/example/")
        fullpath = joinpath(packagepath, "src/example/", dataname)
    end

    if dataname in cd(readdir, path)
        println("The input system: $dataname")
    else
        throw(DomainError(dataname, "the input DATA is not found"))
    end

    if extension == ".h5"
        list = h5read(fullpath, "/H")
        jacobian = sparse(list[:,1], list[:,2], list[:,3])
        jacobianT = sparse(list[:,2], list[:,1], list[:,3])

        observation = h5read(fullpath, "/z")
        noise = h5read(fullpath, "/v")
    elseif extension == ".csv"
        data = DataFrame(load(fullpath))

        list = dropmissing!(data[:,[1,2,3]])
        jacobian = sparse(list.row, list.column, list.coefficient)
        jacobianT = sparse(list.column, list.row, list.coefficient)

        observation = dropmissing!(data[:,[4]]).observation
        noise = dropmissing!(data[:,[5]]).variance
    else
        error("The input data is not a valid format.")
    end

    return SystemModel(jacobian, jacobianT, observation, noise, dataname)
end

#### Read system model from in-Julia
function inJulia(args)
    if typeof(args[1]) == Array{Float64, 2}
        jacobian = sparse(args[1])
    else
        jacobian = args[1]
    end

    jacobianT = copy(transpose(jacobian))
    observation = args[2]
    noise = args[3]
    
    return SystemModel(jacobian, jacobianT, observation, noise, "noname")
end


#### Check keyword arguments 
function checkKeywords(maxIter, damp, bump, prob, alpha, mean, variance, algorithm, out)
    ########## Check iteration scheme ##########
    if maxIter < 0
        error("The upper limit on GBP iterations 'max' has invalid value.")
    end
    if damp > maxIter || damp < 0
        error("Applied randomized damping parameter damp has invalid value.")
    end
    if bump > maxIter || bump <= 0
        error("Applied variance computation parameter bump has invalid value.")
    end

    ########## Desing iteration scheme ##########
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

    ########## Convergence Parameters ##########
    if prob < 0.0 || prob > 1.0
        error("Invalid prob value.")
    end

    if alpha < 0.0 || alpha > 1.0
        error("Invalid alpha value.")
    end

    ########## Virtual Factor Nodes ##########
    if variance < 0.0
        error("Invalid variance value.")
    end

    ########## GBP algorithms ##########
    if !(algorithm in ["gbp", "efficient", "kahan"])
        error("Invalid algorithm key.")
    end

    ########## Return and display data ##########
    if !isa(out, Array)
        out = [out]
    end
    out = vec(out)
    iterate = false; loss = false; wls = false; displayshow = false
    for i in out
        if i == "iterate"
            iterate = true 
        end
        if i == "error"
            loss = true
        end
        if i == "wls"
            wls = true
        end
        if i == "display"
            displayshow = true
        end

    end

    return Settings(maxIter, IterNative, IterDamp, IterBump, IterDampBump, prob, alpha, mean, variance, algorithm, iterate, loss, wls, displayshow)
end