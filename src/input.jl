struct SystemModel
    J::SparseMatrixCSC{Float64,Int64}
    Jt::SparseMatrixCSC{Float64,Int64}
    z::Array{Float64,1}
    v::Array{Float64,1}
    data::String
end


#### Read system model
function model(dataname, maxIter, damp, bump, method, algorithm, path)
    pathtoGaussBP = Base.find_package("GaussBP")
    if pathtoGaussBP == nothing
        throw(ErrorException("GaussBP not found in install packages"))
    end

    extension = ""
    try
        extension = match(r"\.[A-Za-z0-9]+$", dataname).match
    catch
        extension = ""
    end
    if isempty(extension)
        throw(ErrorException("the input DATA extension is not found"))
    elseif extension != ".h5" && extension != ".csv"
        throw(DomainError(extension, "the input DATA extension is not supported"))
    end

    if maxIter < 0
        error("The upper limit on BP iterations max has invalid value.")
    end
    if damp > maxIter || damp < 0
        error("Applied randomized damping parameter damp has invalid value.")
    end
    if bump > maxIter || bump <= 0
        error("Applied variance computation parameter bump has invalid value.")
    end

    if !(method in ["passing", "recursion"])
        error("Invalid method key.")
    end

    if !(algorithm in ["sum", "kahan"])
        error("Invalid algorithm key.")
    end

    if path == "from_package"
        package_dir = abspath(joinpath(dirname(Base.find_package("GaussBP")), ".."))
        path = joinpath(package_dir, "src/example/")
    end
    system = string(path, dataname)

    if extension == ".h5"
        list = h5read(system, "/H")
        jacobian = sparse(list[:,1], list[:,2], list[:,3])
        jacobianT = sparse(list[:,2], list[:,1], list[:,3])

        observation = h5read(system, "/b")
        noise = h5read(system, "/v")
    elseif extension == ".csv"
        data = DataFrame(load(system))

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
