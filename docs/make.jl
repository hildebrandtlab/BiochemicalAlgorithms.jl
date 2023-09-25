using BiochemicalAlgorithms
using Documenter

DocMeta.setdocmeta!(BiochemicalAlgorithms, :DocTestSetup, :(using BiochemicalAlgorithms); recursive=true)

makedocs(;
    modules=[BiochemicalAlgorithms],
    authors="Andreas Hildebrandt <andreas.hildebrandt@uni-mainz.de> and contributors",
    repo=Remotes.GitHub("hildebrandtlab", "BiochemicalAlgorithms.jl"),
    sitename="BiochemicalAlgorithms.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl",
        assets=String[],
    ),
    pages = Any[
        "Home" => "index.md",
        "Library" => Any[
            "System representation" => "public/system.md",
            "Internals" => Any[
                "System representation" => "private/system.md"
            ]
        ]
    ],
)

deploydocs(;
    repo="github.com/hildebrandtlab/BiochemicalAlgorithms.jl",
    devbranch="develop",
)
