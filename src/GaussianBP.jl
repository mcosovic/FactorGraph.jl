module GaussianBP

export bp, bpn

using SparseArrays
using HDF5
using Random
using Printf

############
# Includes #
############

include("input.jl")
include("factorgraph.jl")
include("auxiliary.jl")
include("initialize.jl")
include("inference.jl")
include("summation.jl")
include("printdata.jl")
include("evaluation.jl")

###############
# Gaussian BP #
###############

include("simplybp.jl")
include("neumaierbp.jl")

end # SimplyGBP
