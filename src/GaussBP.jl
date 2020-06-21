module GaussBP

export bp

using SparseArrays, LinearAlgebra
using HDF5
using CSVFiles, DataFrames
using Random
using PrettyTables
using Printf


### Variable Name Style
# prefix "N" stands for the amount of the stem
# prefix "M" stands for the mean of the stem
# prefix "V" stands for the variance of the stem
# prefix "W" stands for the inverse of variance (weight) of the stem
# abbreviation "fac" stands for factor node
# abbreviation "var" stands for variable node
# abbreviation "dir" stands for the direct (singly-connected) factor nodes
# abbreviation "ind" stands for the indirect factor nodes
# sufix "Inv" stands for the inverse value of the stem


### Includes
include("input.jl")
include("routine.jl")
include("inference.jl")

### Run package
function bp(
    data::String = "data33_14.h5";
    max::Int64 = 30,
    damp::Int64 = 10,
    bump::Int64 = max,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    mean::Float64 = 0.0,
    variance::Float64 = 1e5,
    method::String = "passing",
    algorithm::String = "sum",
    path::String = "from_package",
    wls::String = "no")

    system = model(data, max, damp, bump, method, algorithm, path)

  prep = @elapsed begin
    graph = factors(system, mean, variance)
    bp = initialize(graph, max, damp, bump, prob, alpha)
  end

  infe = @elapsed begin
    if method == "passing" && algorithm == "sum"
        for iter = 1:bp.IterNative
            factor_to_variable(graph, bp)
            variable_to_factor(graph, bp)
        end
        for iter = 1:bp.IterDamp
            factor_to_variable_damp(graph, bp)
            variable_to_factor(graph, bp)
        end
        for iter = 1:bp.IterBump
            factor_to_variable_mean(graph, bp)
            variable_to_factor_mean(graph, bp)
        end
        for iter = 1:bp.IterDampBump
            factor_to_variable_mean_damp(graph, bp)
            variable_to_factor_mean(graph, bp)
        end
    end

    if method == "passing" && algorithm == "kahan"
        for iter = 1:bp.IterNative
            factor_to_variable_kahan(graph, bp)
            variable_to_factor_kahan(graph, bp)
        end
        for iter = 1:bp.IterDamp
            factor_to_variable_kahan_damp(graph, bp)
            variable_to_factor_kahan(graph, bp)
        end
        for iter = 1:bp.IterBump
            factor_to_variable_mean(graph, bp)
            variable_to_factor_mean(graph, bp)
        end
        for iter = 1:bp.IterDampBump
            factor_to_variable_mean_damp(graph, bp)
            variable_to_factor_mean(graph, bp)
        end
    end

    if method == "recursion" && algorithm == "sum"
        factor_to_variable(graph, bp)
        rec = recursionini(graph)
        for iter = 1:bp.IterNative
            factor_recursion(graph, bp, rec)
        end
        for iter = 1:bp.IterDamp
            factor_recursion_damp(graph, bp, rec)
        end
        for iter = 1:bp.IterBump
            factor_recursion_mean(graph, bp, rec)
        end
        for iter = 1:bp.IterDampBump
            factor_recursion_mean_damp(graph, bp, rec)
        end
    end

    Xbp = marginal(graph, bp)
  end

    results(system, graph, bp, Xbp, wls, prep, infe)

    return Xbp, system
end

end # GaussBP
