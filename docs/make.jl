using Documenter, MetidaNCA
#using DocumenterLaTeX

makedocs(
        modules = [MetidaNCA],
        sitename = "MetidaNCA.jl",
        authors = "Vladimir Arnautov",
        pages = [
            "Home" => "index.md",
            "Examples" => "examples.md",
            "Details" => "details.md",
            "Parameters" => "parameters.md",
            "API" => "api.md"
            ],
        )
deploydocs(repo = "github.com/PharmCat/MetidaNCA.jl.git", devbranch = "main", forcepush = true
)
