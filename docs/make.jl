using BALL
using Documenter

DocMeta.setdocmeta!(BALL, :DocTestSetup, :(using BALL); recursive=true)

makedocs(;
    modules=[BALL],
    authors="Andreas Hildebrandt <andreas.hildebrandt@uni-mainz.de> and contributors",
    repo="https://github.com/hildebrandtlab/BALL.jl/blob/{commit}{path}#{line}",
    sitename="BALL.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hildebrandtlab.github.io/BALL.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hildebrandtlab/BALL.jl",
    devbranch="main",
)
