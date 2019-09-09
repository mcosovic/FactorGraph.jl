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

###########
# BP Time #
###########


# function bptime(data::String="data_33_14", MAXI::Int64=15, DAMP::Int64=10, PROB::Float64=0.6, ALPH::Float64=0.4, MEAN::Float64=0.0, VARI::Float64=1e3)
#     xbp, fgraph, init, infe, solu = bp(data, MAXI, DAMP, PROB, ALPH, MEAN, VARI)
#     bp_time(fgraph, init, infe, solu)
#
#     return xbp
# end
#
# function bpntime(data::String="data_33_14", MAXI::Int64=15, DAMP::Int64=10, PROB::Float64=0.6, ALPH::Float64=0.4, MEAN::Float64=0.0, VARI::Float64=1e3)
#     xbp, fgraph, init, infe, solu = bpn(data, MAXI, DAMP, PROB, ALPH, MEAN, VARI)
#     bp_time(fgraph, init, infe, solu)
#
#     return xbp
# end

#############
# BP vs WLS #
#############




end # SimplyGBP
