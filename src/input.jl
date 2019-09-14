#######################
#  Read system model  #
#######################


#------------------------------------------------------------------------
function model(DATA, PATH)
    system = string(PATH, DATA, ".h5")

    list = h5read(system, "/H")
    jacobian = sparse(list[:,1], list[:,2], list[:,3])

    observation = h5read(system, "/b")
    noise = h5read(system, "/v")

    return jacobian, observation, noise
end
#------------------------------------------------------------------------
