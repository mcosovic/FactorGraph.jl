using Documenter, FactorGraph

makedocs(
    sitename = "FactorGraph",
    modules = [FactorGraph],
    clean = true,
    doctest = false,
    format = Documenter.HTML(
        assets = ["assets/style.css"],
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1,
        repolink = "https://github.com/mcosovic/FactorGraph.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Gaussian Models" => [
            "Factor Graph" => "gaussian_models/gauss_factor_graph.md",
            "Iterative Belief Propagation" => "gaussian_models/gauss_iterative.md",
            "Forward-Backward Belief Propagation" => "gaussian_models/gauss_forward_backward.md",
        ],
        "Discrete Models" => [
            "Factor Graph" => "discrete_models/discrete_factor_graph.md",
            "Iterative Belief Propagation" => "discrete_models/discrete_iterative.md",
            "Forward-Backward Belief Propagation" => "discrete_models/discrete_forward_backward.md",
        ],
        "Examples" => [
            "DC State Estimation" => "examples/dc_state_estimation.md",
            "PMU-Based State Estimation" => "examples/pmu_state_estimation.md",
            "Robot Localization on a Chain" => "examples/robot_localization.md",
            "Protection Alarm Diagnosis" => "examples/protection_alarm_diagnosis.md",
            "Binary Pump Diagnosis" => "examples/binary_pump_diagnosis.md",
            "Discrete Robot Localization" => "examples/discrete_robot_localization.md",
            "LDPC-Style Decoding" => "examples/ldpc_decoding.md",
        ],
        "API" => [
            "Types" => "api/types.md",
            "Factor Graphs" => "api/factor_graphs.md",
            "Inference" => "api/inference.md",
            "Validation" => "api/validation.md",
            "Print Functions" => "api/printing.md",
        ],
        "Release Notes" => "release_notes.md",
    ],
    repo = "https://github.com/mcosovic/FactorGraph.jl/{commit}{path}#L{line}",
    authors = "Mirsad Cosovic",
)

deploydocs(
    repo = "github.com/mcosovic/FactorGraph.jl.git",
    target = "build",
    devbranch = "master",
)
