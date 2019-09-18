module Instability

using SparseArrays
using HDF5
using Test
using Random
using PrettyTables


#################################################
#  The belief propagation type instability test #
#################################################


##################################
#  Read data from SimpleTest.h5  #
##################################
function model(DATA, PATH)
    system = string(PATH, DATA, ".h5")

    list = h5read(system, "/H")
    jacobian = sparse(list[:,1], list[:,2], list[:,3])

    observation = h5read(system, "/b")
    noise = h5read(system, "/v")

    return jacobian, observation, noise
end

jacobian, observation, noise = model("SimpleTest", "test/")


################################################
#  Check for type stability of main functions  #
################################################
include("../src/factorgraph.jl")
include("../src/auxiliary.jl")
include("../src/initialize.jl")
include("../src/inference.jl")
include("../src/evaluation.jl")
include("../src/bp_simple_passing.jl")
include("../src/bp_kahan_passing.jl")
include("../src/bp_simple_recursion.jl")

@inferred bp_simple_passing(jacobian, observation, noise, 10, 5, 10, 0.6, 0.5, 0.0, 1e-3, "on", "on", "on")
@inferred bp_kahan_passing(jacobian, observation, noise, 10, 5, 10, 0.6, 0.5, 0.0, 1e-3, "on", "on", "on")
@inferred bp_simple_recursion(jacobian, observation, noise, 10, 5, 10, 0.6, 0.5, 0.0, 1e-3, "on", "on", "on")

end # Instability
