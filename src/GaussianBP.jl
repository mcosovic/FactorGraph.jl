module GaussianBP

export bp, bpn, bptime

using SparseArrays
using HDF5
using Random
using Printf


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

###########
# BP Time #
###########

include("printdata.jl")

function bptime(data::String="data_33_14", MAXI::Int64=15, DAMP::Int64=10, PROB::Float64=0.6, ALPH::Float64=0.4, MEAN::Float64=0.0, VARI::Float64=1e3)
    xbp, fgraph, init, infe, solu = bp(data, MAXI, DAMP, PROB, ALPH, MEAN, VARI)
    bp_time(fgraph, init, infe, solu)

    return xbp
end

end # SimplyGBP
