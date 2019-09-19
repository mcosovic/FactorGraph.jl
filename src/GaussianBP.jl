module GaussianBP

export bp

using SparseArrays
using HDF5
using CSV
using Random
using PrettyTables

#########################
#  Variable Name Style  #
#########################
# prefix "N" stands for the amount of the stem
# prefix "M" stands for the mean of the stem
# prefix "V" stands for the variance of the stem
# prefix "W" stands for the inverse of variance (weight) of the stem
# abbreviation "fac" stands for factor node
# abbreviation "var" stands for variable node
# abbreviation "dir" stands for the direct (singly-connected) factor nodes
# abbreviation "ind" stands for the indirect factor nodes
# sufix "Inv" stands for the inverse value of the stem


##############
#  Includes  #
##############
include("input.jl")
include("error.jl")
include("auxiliary.jl")
include("evaluation.jl")

include("factorgraph.jl")
include("initialize.jl")
include("inference.jl")
include("inference_mean.jl")

include("bp_simple_passing.jl")
include("bp_kahan_passing.jl")
include("bp_simple_recursion.jl")

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
    METHOD::String = "passing",
    ALGORITHM::String = "sum",
    TIME::String = "off",
    ERROR::String = "off",
    STATISTIC::String = "off",
    PATH::String = "from_package")


    check_iteration_scheme(MAXI, DAMP, BUMP)

    jacobian, observation, noise = model(DATA, PATH)

    if METHOD == "passing" && ALGORITHM == "sum"
        Xbp = bp_simple_passing(jacobian, observation, noise,
                                MAXI, DAMP, BUMP, PROB, ALPH, MEAN, VARI,
                                TIME, ERROR, STATISTIC)
    end
    if METHOD == "passing" && ALGORITHM == "kahan"
        Xbp = bp_kahan_passing(jacobian, observation, noise,
                               MAXI, DAMP, BUMP, PROB, ALPH, MEAN, VARI,
                               TIME, ERROR, STATISTIC)
    end
    if METHOD == "recursion" && ALGORITHM == "sum"
        Xbp = bp_simple_recursion(jacobian, observation, noise,
                                MAXI, DAMP, BUMP, PROB, ALPH, MEAN, VARI,
                                TIME, ERROR, STATISTIC)
    end

    return Xbp
end

end # SimplyGBP
