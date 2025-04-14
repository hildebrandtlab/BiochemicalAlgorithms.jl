using BiochemicalAlgorithms
using Documenter

DocMeta.setdocmeta!(BiochemicalAlgorithms, :DocTestSetup, :(using BiochemicalAlgorithms); recursive=true)

const pages = Any[
    "Home" => "index.md",
    "Tutorials" => [
        "How to get started" => "tutorials/01_QuickStart.md",
        "Introduction to core components" => "tutorials/02_CoreComponents.md",
        "How to iterate" => "tutorials/03_Iterate.md",
        "How to apply force fields" => "tutorials/04_ForceFields.md",
        "How to handle molecules" =>"tutorials/05_HandleMolecules.md",
        "How to read and write" => "tutorials/06_ReadAndWrite.md",

    ],
    "Library" => Any[
        "Biomolecular systems" => "public/system.md",
        "File formats" => "public/fileformats.md",
        "Force fields" => "public/forcefields.md",
        "Mappings" => "public/mappings.md",
        "Molecular preprocessing" => "public/preprocessing.md",
        "Internals" => Any[
            "File Formats" => "private/fileformats.md",
            "Biomolecular systems" => "private/system.md",
            "Mappings" => "private/mappings.md"
        ]
    ]
]

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
    pages = pages,
)

deploydocs(;
    repo="github.com/hildebrandtlab/BiochemicalAlgorithms.jl",
    devbranch="develop",
)
