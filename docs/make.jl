using Documenter, MetidaNCA, Weave, PrettyTables, CSV, DataFrames
#using DocumenterLaTeX

include("validation.jl")

#v_out_path = joinpath(dirname(@__FILE__), "src", "validation_report.md")

makedocs(
        modules = [MetidaNCA],
        sitename = "MetidaNCA.jl",
        authors = "Vladimir Arnautov",
        pages = [
            "Home" => "index.md",
            "Examples" => "examples.md",
            "Details" => "details.md",
            "Parameters" => "parameters.md",
            "API" => "api.md",
            "Validation" => "validation_report.md"
            ],
        )
##=
makedocs(
    format = Documenter.LaTeX(),
    modules = [MetidaNCA],
    sitename = "MetidaNCA.jl",
    authors = "Vladimir Arnautov",
    pages = [
        "Validation" => "validation_report.md"
        ]
    )

##=
deploydocs(repo = "github.com/PharmCat/MetidaNCA.jl.git", devbranch = "main", forcepush = true
)
