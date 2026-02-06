using Documenter, CVC

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://dohyunkim116.github.io/CVC.jl",
        assets = String[],
    ),
    sitename = "CVC.jl",
    authors = "Do Hyun Kim",
    modules = [CVC],
    clean = true,
    doctest = false,
    checkdocs = :none,
    pages = [
        "Introduction" => "index.md",
        "Model Fitting" => "model_fitting.md",
    ]
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/dohyunkim116/CVC.jl.git",
        target = "build",
        branch = "gh-pages",
        devbranch = "main",
        push_preview = true,
    )
end
