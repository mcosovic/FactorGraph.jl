module GaussBP

using LinearAlgebra: Algorithm
export gbp

using SparseArrays, LinearAlgebra
using HDF5, XLSX
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
include("results.jl")
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
    algorithm::String = "vanilla",
    out::Union{String, Array{String,1}} = "no")
    

    settings = check_keywords(max, damp, bump, prob, alpha, mean, variance, algorithm, out)
    if check_inputs(args)
        system = julia_out(args, settings)
    else
        system = julia_in(args, settings)
    end
 
    graph = factors(system, settings)
    bp = initialize(graph, settings)
    results = initialize_results(settings, system, graph)

    display_stat(graph, bp, settings, system, algorithm)
    
    if algorithm == "vanilla" && !settings.outIterate
        for iter = 1:bp.iterNative
            vanilla_factor_to_variable(graph, bp)
            vanilla_variable_to_factor(graph, bp)
        end
        for iter = 1:bp.iterDamp
            vanilla_factor_to_variable_damp(graph, bp)
            vanilla_variable_to_factor(graph, bp)
        end
        for iter = 1:bp.iterBump
            vanilla_factor_to_variable_mean(graph, bp)
            vanilla_variable_to_factor_mean(graph, bp)
        end
        for iter = 1:bp.iterDampBump 
            vanilla_factor_to_variable_mean_damp(graph, bp)
            vanilla_variable_to_factor_mean(graph, bp)
        end
        marginal(settings, system, graph, bp, results)
    end

    if algorithm == "vanilla" && settings.outIterate || algorithm == "vanillaDynamic" 
        for iter = 1:bp.iterNative
            graph_dynamic(settings, system, graph, bp)
            vanilla_factor_to_variable(graph, bp)
            vanilla_variable_to_factor(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterDamp
            graph_dynamic(settings, system, graph, bp)
            vanilla_factor_to_variable_damp(graph, bp)
            vanilla_variable_to_factor(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterBump
            graph_dynamic(settings, system, graph, bp)
            vanilla_factor_to_variable_mean(graph, bp)
            vanilla_variable_to_factor_mean(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterDampBump
            graph_dynamic(settings, system, graph, bp)
            vanilla_factor_to_variable_mean_damp(graph, bp)
            vanilla_variable_to_factor_mean(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
    end

    if algorithm == "efficient" && !settings.outIterate
        for iter = 1:bp.iterNative
            efficient_factor_to_variable(graph, bp)
            efficient_variable_to_factor(graph, bp)
        end
        for iter = 1:bp.iterDamp
            efficient_factor_to_variable_damp(graph, bp)
            efficient_variable_to_factor(graph, bp)
        end
        for iter = 1:bp.iterBump
            efficient_factor_to_variable_mean(graph, bp)
            efficient_variable_to_factor_mean(graph, bp)
        end
        for iter = 1:bp.iterDampBump
            efficient_factor_to_variable_mean_damp(graph, bp)
            efficient_variable_to_factor_mean(graph, bp)
        end
        marginal(settings, system, graph, bp, results)
    end

    if algorithm == "efficient" && settings.outIterate || algorithm == "efficientDynamic" 
        for iter = 1:bp.iterNative
            graph_dynamic(settings, system, graph, bp)
            efficient_factor_to_variable(graph, bp)
            efficient_variable_to_factor(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterDamp
            graph_dynamic(settings, system, graph, bp)
            efficient_factor_to_variable_damp(graph, bp)
            efficient_variable_to_factor(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterBump
            graph_dynamic(settings, system, graph, bp)
            efficient_factor_to_variable_mean(graph, bp)
            efficient_variable_to_factor_mean(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterDampBump
            graph_dynamic(settings, system, graph, bp)
            efficient_factor_to_variable_mean_damp(graph, bp)
            efficient_variable_to_factor_mean(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
    end

    if algorithm == "kahan" && !settings.outIterate
        for iter = 1:bp.iterNative
            kahan_factor_to_variable(graph, bp)
            kahan_variable_to_factor(graph, bp)
        end
        for iter = 1:bp.iterDamp
            kahan_factor_to_variable_damp(graph, bp)
            kahan_variable_to_factor(graph, bp)
        end
        for iter = 1:bp.iterBump
            kahan_factor_to_variable_mean(graph, bp)
            kahan_variable_to_factor_mean(graph, bp)
        end
        for iter = 1:bp.iterDampBump
            kahan_factor_to_variable_mean_damp(graph, bp)
            kahan_variable_to_factor_mean(graph, bp)
        end
        marginal(settings, system, graph, bp, results)
    end

    if algorithm == "kahan" && settings.outIterate || algorithm == "kahanDynamic" 
        for iter = 1:bp.iterNative
            graph_dynamic(settings, system, graph, bp)
            kahan_factor_to_variable(graph, bp)
            kahan_variable_to_factor(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterDamp
            graph_dynamic(settings, system, graph, bp)
            kahan_factor_to_variable_damp(graph, bp)
            kahan_variable_to_factor(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterBump
            graph_dynamic(settings, system, graph, bp)
            kahan_factor_to_variable_mean(graph, bp)
            kahan_variable_to_factor_mean(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
        for iter = 1:bp.iterDampBump
            graph_dynamic(settings, system, graph, bp)
            kahan_factor_to_variable_mean_damp(graph, bp)
            kahan_variable_to_factor_mean(graph, bp)
            marginal(settings, system, graph, bp, results)
        end
    end
    
    if settings.outWls
        wls_metrics(system, settings, graph, results)
    end

    if settings.outDisplay
        display_data(graph, bp, results, settings)
    end

    return results, system
end

end # GaussBP

