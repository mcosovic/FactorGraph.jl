module GaussBP

export gbp

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
include("vanillaGBP.jl")
include("efficientGBP.jl")
include("kahanGBP.jl")

### Run package
function gbp(
    args...;
    max::Int64 = 30,
    damp::Int64 = max,
    bump::Int64 = max,
    prob::Float64 = 0.6,
    alpha::Float64 = 0.4,
    mean::Float64 = 0.0,
    variance::Float64 = 1e5,
    algorithm::String = "gbp",
    out::Union{String, Array{String,1}} = "no")
    
    if checkInputs(args)
        system = outJulia(args)
    else
        system = inJulia(args)
    end
    settings = checkKeywords(max, damp, bump, prob, alpha, mean, variance, algorithm, out)

    preprocessing = @elapsed begin
        graph = factors(system, settings)
        bp = initialize(graph, settings)
    end
    results = resultsdata(system, settings, graph.Nvar)

    inference = @elapsed begin
        if algorithm == "gbp" && !settings.iterate
            for iter = 1:settings.IterNative
                vanilla_factor_to_variable(graph, bp)
                vanilla_variable_to_factor(graph, bp)
            end
            for iter = 1:settings.IterDamp
                vanilla_factor_to_variable_damp(graph, bp)
                vanilla_variable_to_factor(graph, bp)
            end
            for iter = 1:settings.IterBump
                vanilla_factor_to_variable_mean(graph, bp)
                vanilla_variable_to_factor_mean(graph, bp)
            end
            for iter = 1:settings.IterDampBump
                vanilla_factor_to_variable_mean_damp(graph, bp)
                vanilla_variable_to_factor_mean(graph, bp)
            end
            marginal(graph, bp, results)
        end

        if algorithm == "gbp" && settings.iterate
            for iter = 1:settings.IterNative
                vanilla_factor_to_variable(graph, bp)
                vanilla_variable_to_factor(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterDamp
                vanilla_factor_to_variable_damp(graph, bp)
                vanilla_variable_to_factor(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterBump
                vanilla_factor_to_variable_mean(graph, bp)
                vanilla_variable_to_factor_mean(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterDampBump
                vanilla_factor_to_variable_mean_damp(graph, bp)
                vanilla_variable_to_factor_mean(graph, bp)
                marginal(graph, bp, results)
            end
        end

        if algorithm == "efficient" && !settings.iterate
            for iter = 1:settings.IterNative
                efficient_factor_to_variable(graph, bp)
                efficient_variable_to_factor(graph, bp)
            end
            for iter = 1:settings.IterDamp
                efficient_factor_to_variable_damp(graph, bp)
                efficient_variable_to_factor(graph, bp)
            end
            for iter = 1:settings.IterBump
                efficient_factor_to_variable_mean(graph, bp)
                efficient_variable_to_factor_mean(graph, bp)
            end
            for iter = 1:settings.IterDampBump
                efficient_factor_to_variable_mean_damp(graph, bp)
                efficient_variable_to_factor_mean(graph, bp)
            end
            marginal(graph, bp, results)
        end
    
        if algorithm == "efficient" && settings.iterate
            for iter = 1:settings.IterNative
                efficient_factor_to_variable(graph, bp)
                efficient_variable_to_factor(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterDamp
                efficient_factor_to_variable_damp(graph, bp)
                efficient_variable_to_factor(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterBump
                efficient_factor_to_variable_mean(graph, bp)
                efficient_variable_to_factor_mean(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterDampBump
                efficient_factor_to_variable_mean_damp(graph, bp)
                efficient_variable_to_factor_mean(graph, bp)
                marginal(graph, bp, results)
            end
        end

        if algorithm == "kahan" && !settings.iterate
            for iter = 1:settings.IterNative
                kahan_factor_to_variable(graph, bp)
                kahan_variable_to_factor(graph, bp)
            end
            for iter = 1:settings.IterDamp
                kahan_factor_to_variable_damp(graph, bp)
                kahan_variable_to_factor(graph, bp)
            end
            for iter = 1:settings.IterBump
                kahan_factor_to_variable_mean(graph, bp)
                kahan_variable_to_factor_mean(graph, bp)
            end
            for iter = 1:settings.IterDampBump
                kahan_factor_to_variable_mean_damp(graph, bp)
                kahan_variable_to_factor_mean(graph, bp)
            end
            marginal(graph, bp, results)
        end

        if algorithm == "kahan" && settings.iterate
            for iter = 1:settings.IterNative
                kahan_factor_to_variable(graph, bp)
                kahan_variable_to_factor(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterDamp
                kahan_factor_to_variable_damp(graph, bp)
                kahan_variable_to_factor(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterBump
                kahan_factor_to_variable_mean(graph, bp)
                kahan_variable_to_factor_mean(graph, bp)
                marginal(graph, bp, results)
            end
            for iter = 1:settings.IterDampBump
                kahan_factor_to_variable_mean_damp(graph, bp)
                kahan_variable_to_factor_mean(graph, bp)
                marginal(graph, bp, results)
            end
        end
    end

    stats(system, results, settings)
    if settings.displayshow
        displaydata(graph, bp, results, settings, preprocessing, inference)
    end

    return results, system
end

end # GaussBP



