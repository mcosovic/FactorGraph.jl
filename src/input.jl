#######################
#  Read system model  #
#######################


#-------------------------------------------------------------------------------
function model(DATA, PATH)
    if PATH == "from_package"
        package_dir = abspath(joinpath(dirname(Base.find_package("GaussianBP")), ".."))
        PATH = joinpath(package_dir, "src/data/")
    end

    system = string(PATH, DATA, ".h5")

    list = h5read(system, "/H")
    jacobian = sparse(list[:,1], list[:,2], list[:,3])

    observation = h5read(system, "/b")
    noise = h5read(system, "/v")

    return jacobian, observation, noise
end
#-------------------------------------------------------------------------------
