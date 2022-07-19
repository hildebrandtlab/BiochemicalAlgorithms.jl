using BiochemicalAlgorithms
using Documenter

DocMeta.setdocmeta!(BiochemicalAlgorithms, :DocTestSetup, :(using BiochemicalAlgorithms); recursive=true)

makedocs(;
    modules=[BiochemicalAlgorithms],
    authors="Andreas Hildebrandt <andreas.hildebrandt@uni-mainz.de> and contributors",
    repo="https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/blob/{commit}{path}#{line}",
    sitename="BiochemicalAlgorithms.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hildebrandtlab/BiochemicalAlgorithms.jl",
    devbranch="develop",
)
