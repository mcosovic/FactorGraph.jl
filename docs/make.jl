using Documenter, GaussBP

makedocs(
    sitename = "GaussBP",
    modules = [GaussBP],
    clean = false,
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Input Data" => "man/input.md",
        "Graphical Model" => "man/graphicalmodel.md",
        "Inference" => "man/inference.md",
        "Utility Functions" => "man/utility.md",
        "Output Data" => "man/output.md",
        "Theoretical Background" => "man/theoretical.md",
    ],
)

deploydocs(
    repo = "github.com/mcosovic/GaussBP.jl.git",
    target = "build",
)
