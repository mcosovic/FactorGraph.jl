struct Settings
    maxIter::Int64
    damp::Int64
    bump::Int64
    prob::Float64 
    alpha::Float64 
    mean::Float64 
    variance::Float64
    dynamic::Bool
    outIterate::Bool 
    outEvaluation::Bool
    outWls::Bool
    outDisplay::Bool
    packagepath::String
end

struct SystemModel
    jacobian::SparseMatrixCSC{Float64,Int64}
    jacobianTranspose::SparseMatrixCSC{Float64,Int64}
    observation::Array{Float64,1}
    variance::Array{Float64,1}
    dynamic::Array{Float64,2}
    data::String
end

@enum AlgorithmType vanilla efficient kahan vanillaDynamic efficientDynamic kahanDynamic
@enum OutputControl iteration evaluation wls terminal


########## Check keyword arguments ##########
function check_keywords(maxIter, damp, bump, prob, alpha, mean, variance, algorithm, out)
    #### Check the package is installed 
    pathtoGaussBP = Base.find_package("GaussBP")
    if isnothing(pathtoGaussBP)
        throw(ErrorException("GaussBP not found in install packages"))
    end
    packagepath = abspath(joinpath(dirname(pathtoGaussBP), ".."))

    #### Check iteration scheme 
    if maxIter < 0
        error("The upper limit on GBP iterations 'max' has invalid value.")
    end
    if damp > maxIter || damp < 0
        error("Applied randomized damping parameter damp has invalid value.")
    end
    if bump > maxIter || bump < 0
        error("Applied variance computation parameter bump has invalid value.")
    end

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
   
    #### Check dynamic model
    dynamic = false
    if algorithm in [vanillaDynamic, efficientDynamic, kahanDynamic]
        dynamic = true
    end

    #### Check output data
    if !isa(out, Array)
        out = [out]
    end
    out = vec(out)
        
    outIterate = false; outEvaluation = false; outWls = false; outDisplay = false
    for i in out
        if i == iteration
            outIterate = true 
        end
        if i == evaluation
            outEvaluation = true
        end
        if i == wls
            outWls = true
        end
        if i == terminal
            outDisplay = true
        end
    end

    return Settings(maxIter, damp, bump, prob, alpha, mean, variance, dynamic, outIterate, outEvaluation, outWls, outDisplay, packagepath)
end

########## Type of the system data ########## 
function check_inputs(args)
    fromfile = false
    for i in args
        if typeof(i) == String
            fromfile = true
            break
        end
    end

    return fromfile
end

########## Load from HDF5 or XLSX files ########## 
function julia_out(args, settings)
    extension = ".h5"; path = ""; dataname = ""; fullpath = ""

    for i = 1:length(args)
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
        path = joinpath(settings.packagepath, "src/example/")
        fullpath = joinpath(settings.packagepath, "src/example/", dataname)
    end

    if !(dataname in cd(readdir, path))
        throw(DomainError(dataname, "the input DATA is not found"))
    end

    #### Read from HDF5 or XLSX file
    dynamic = Array{Float64}(undef, 0, 0)
    if extension == ".h5"
        list = h5read(fullpath, "/jacobian")::Array{Float64,2}
        jacobian = sparse(list[:,1], list[:,2], list[:,3])::SparseMatrixCSC{Float64,Int64}
        jacobianTranspose = sparse(list[:,2], list[:,1], list[:,3])::SparseMatrixCSC{Float64,Int64}

        observation = h5read(fullpath, "/observation")::Array{Float64,1}
        variance = h5read(fullpath, "/variance")::Array{Float64,1}

        if settings.dynamic
            dynamic = sortslices(h5read(fullpath, "/dynamic"), dims = 1, by = x->x[1], rev = false)::Array{Float64,2}
            dynamic_start(dynamic[1, 1])
        end
    elseif extension == ".xlsx"
        xf = XLSX.openxlsx(fullpath, mode = "r")
        if "jacobian" in XLSX.sheetnames(xf)
            start = startxlsx(xf["jacobian"])
            list = xf["jacobian"][:][start:end, :]::Array{Float64,2}
            jacobian = sparse(list[:, 1], list[:, 2], list[:, 3])::SparseMatrixCSC{Float64,Int64}
            jacobianTranspose = sparse(list[:, 2], list[:, 1], list[:, 3])::SparseMatrixCSC{Float64,Int64}
        else
            throw(ErrorException("error opening sheet jacobian"))
        end
        if "measurement" in XLSX.sheetnames(xf)
            start = startxlsx(xf["measurement"])
            list = xf["measurement"][:][start:end, :]::Array{Float64,2}
            observation = list[:, 1]::Array{Float64,1}
            variance = list[:, 2]::Array{Float64,1}
        else
            throw(ErrorException("error opening sheet measurement"))
        end
        if settings.dynamic
            if "dynamic" in XLSX.sheetnames(xf)
                start = startxlsx(xf["dynamic"])
                dynamic = xf["dynamic"][:][start:end, :]::Array{Float64,2}
                dynamic = sortslices(dynamic, dims = 1, by = x->x[1], rev = false)::Array{Float64,2}
                dynamic_start(dynamic[1, 1])
            else
                throw(ErrorException("error opening sheet dynamic"))
            end
        end
    else
        error("the input data is not a valid format")
    end

    return SystemModel(jacobian, jacobianTranspose, observation, variance, dynamic, dataname)
end

########## Read in-Julia system model ##########
function julia_in(args, settings)
    if typeof(args[1]) == Array{Float64, 2}
        jacobian = sparse(args[1])
    else
        jacobian = args[1]
    end

    jacobianTranspose = copy(transpose(jacobian))
    observation = args[2]
    variance = args[3]
 
    dynamic = Array{Float64}(undef, 0, 0)
    if settings.dynamic
        try
            dynamic = sortslices(args[4], dims = 1, by = x->x[1], rev = false)::Array{Float64,2}
        catch
            throw(DomainError("the dynamic update data is not found"))
        end
        dynamic_start(dynamic[1, 1])
    end
    
    return SystemModel(jacobian, jacobianTranspose, observation, variance, dynamic, "noname") 
end

########## Check Dynamic start update point ##########
function dynamic_start(start)
    if start <= 1
        throw(DomainError("the dynamic update point has invalid value"))
    end
end