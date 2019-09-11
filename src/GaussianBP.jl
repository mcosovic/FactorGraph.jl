module GaussianBP

export bp

using SparseArrays
using HDF5
using Random
using PrettyTables


##############
#  Includes  #
##############
include("input.jl")
include("factorgraph.jl")
include("auxiliary.jl")
include("initialize.jl")
include("inference.jl")
include("summation.jl")
include("evaluation.jl")
include("simplybp.jl")
include("neumaierbp.jl")


#################
#  Run package  #
#################
function bp(
    DATA::String = "data33_14",
    MAXI::Int64 = 30,
    DAMP::Int64 = 10,
    BUMP::Int64 = MAXI,
    PROB::Float64 = 0.6,
    ALPH::Float64 = 0.4,
    MEAN::Float64 = 0.0,
    VARI::Float64 = 1e5;
    ALGORITHM::String = "sum",
    TIME::String = "off",
    ERROR::String = "off",
    PATH::String = "src/data/")

    H, b, v = model(DATA, PATH)

    if ALGORITHM == "sum"
        xbp = bps(H, b, v, MAXI, DAMP, BUMP, PROB, ALPH, MEAN, VARI, TIME)
    end
    if ALGORITHM == "kahan"
        xbp = bpn(H, b, v, MAXI, DAMP, BUMP, PROB, ALPH, MEAN, VARI, TIME)
    end

    errors(H, b, v, xbp, ERROR, TIME)

    return xbp
end

end # SimplyGBP
