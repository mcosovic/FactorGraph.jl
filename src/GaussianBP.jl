module GaussianBP

export bp, bpn

using SparseArrays
using HDF5
using Random


######################
# Default Parameters #
######################

# const data = "data_33_14"

# const MAXI = 15
# const DAMP = 10
#
# const PROB = 0.6
# const ALPH = 0.4
#
# const MEAN = 0.0
# const VARI = 1e3


############
# Includes #
############

include("input.jl")
include("factorgraph.jl")
include("auxiliary.jl")
include("initialize.jl")
include("inference.jl")
include("summation.jl")


###############
# Gaussian BP #
###############

include("simplybp.jl")
include("neumaierbp.jl")

end # SimplyGBP
