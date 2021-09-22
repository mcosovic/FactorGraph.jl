using Documenter, FactorGraph

makedocs(
    sitename = "FactorGraph",
    modules = [FactorGraph],
    clean = false,
    doctest = false,
    format = Documenter.HTML(
        assets = ["assets/style.css"],
        prettyurls = prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Continuous Framework" => [
        "Input Data" => "man/continuousInput.md",
        "General Factor Graph" => [
            "Graphical Model" => "man/continuousModel.md",
            "Inference" => "man/continuousInference.md",
            "Output Data" => "man/continuousOutput.md",
        ],
        "Tree Factor Graph" => [
            "Graphical Model" => "man/continuousTreeModel.md",
            "Inference" => "man/continuousTreeInference.md",
            "Output Data" => "man/continuousTreeOutput.md",
        ],
        "Utility Functions" => "man/utility.md"
        ],
        "Discrete Framework" => [
        "Input Data" => "man/discreteInput.md",
        "Tree Factor Graph" => [
            "Graphical Model" => "man/discreteTreeModel.md",
            "Inference" => "man/discreteTreeInference.md",
            "Output Data" => "man/discreteTreeOutput.md",
        ],
        ],
        "Theoretical Background" => "man/theoretical.md",
    ],
)

deploydocs(
    repo = "github.com/mcosovic/FactorGraph.jl.git",
    target = "build",
)
