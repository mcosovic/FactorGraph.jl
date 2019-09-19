#######################
#  Read system model  #
#######################


#-------------------------------------------------------------------------------
function get_extension(DATA)
    try
        match(r"\.[A-Za-z0-9]+$", DATA).match
    catch
        error("The input DATA extansion is missing.")
    end
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
function model(DATA, PATH)
    if PATH == "from_package"
        package_dir = abspath(joinpath(dirname(Base.find_package("GaussianBP")), ".."))
        PATH = joinpath(package_dir, "src/example/")
    end

    extension = get_extension(DATA)

    system = string(PATH, DATA)

    if extension == ".h5"
        list = h5read(system, "/H")
        jacobian = sparse(list[:,1], list[:,2], list[:,3])

        observation = h5read(system, "/b")
        noise = h5read(system, "/v")

    elseif extension == ".csv"
        data = DataFrame(load(system))

        list = dropmissing!(data[:,[1,2,3]])
        jacobian = sparse(list.row, list.column, list.coefficient)

        observation = dropmissing!(data[:,[4]]).observation
        noise = dropmissing!(data[:,[5]]).variance
    else
        error("The input data is not a valid format.")
    end

    return jacobian, observation, noise
end
#-------------------------------------------------------------------------------
