using Documenter, GaussBP

makedocs(
    sitename = "GaussBP",
    modules = [GaussBP],
    clean = false,
    doctest = false,
    format = Documenter.HTML(assets=["assets/style.css"]),
    pages = [
        "Home" => "index.md",
        "Input Data" => "man/input.md",
        "Run Settings" => "man/runsettings.md",
        "Output Data" => "man/output.md",
        "Theoretical Background" => "man/theoretical.md",
    ],
)

deploydocs(
    repo = "github.com/mcosovic/GaussBP.jl.git",
    target = "build",
)
