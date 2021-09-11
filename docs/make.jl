using Documenter, GaussBP

makedocs(
    sitename = "GaussBP",
    modules = [GaussBP],
    clean = false,
    doctest = false,
    format = Documenter.HTML(
        assets = ["assets/style.css"],
        prettyurls = prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1,
    ),
    pages = [
        "Home" => "index.md",
        "Input Data" => "man/input.md",
        "General Factor Graph" => [
            "Graphical Model" => "man/graphicalmodel.md",
            "Inference" => "man/inference.md",
            "Output Data" => "man/output.md",
        ],
        "Tree Factor Graph" => [
            "Graphical Model" => "man/graphicalmodeltree.md",
            "Inference" => "man/inferencetree.md",
            "Output Data" => "man/outputtree.md",
        ],
        "Utility Functions" => "man/utility.md",
        "Theoretical Background" => "man/theoretical.md",
    ],
)

deploydocs(
    repo = "github.com/mcosovic/GaussBP.jl.git",
    target = "build",
)
