using Documenter, MetidaNCA
#using DocumenterLaTeX

makedocs(
        modules = [MetidaNCA],
        sitename = "MetidaNCA.jl",
        authors = "Vladimir Arnautov",
        pages = [
            "Home" => "index.md",

            ],
        )
deploydocs(repo = "github.com/PharmCat/MetidaNCA.jl.git", devbranch = "main"
)
