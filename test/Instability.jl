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

    Hlist = h5read(system, "/H")
    H = sparse(Hlist[:,1], Hlist[:,2], Hlist[:,3])

    b = h5read(system, "/b")
    v = h5read(system, "/v")

    return H, b, v
end

H, b, v = model("SimpleTest", "test/")


################################################
#  Check for type stability of main functions  #
################################################
include("../src/factorgraph.jl")
include("../src/auxiliary.jl")
include("../src/initialize.jl")
include("../src/summation.jl")
include("../src/inference.jl")
include("../src/evaluation.jl")
include("../src/simplybp.jl")
include("../src/neumaierbp.jl")
@inferred bps(H, b, v, 10, 5, 10, 0.6, 0.5, 0.0, 1e-3, "on")
@inferred bpn(H, b, v, 10, 5, 10, 0.6, 0.5, 0.0, 1e-3, "on")


end # Instability
