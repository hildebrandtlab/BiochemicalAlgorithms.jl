#!/usr/bin/env julia
#
#

if "--help" ∈ ARGS
    println(
        """
docs/make.jl

Render the `BiochemicalAlgorithms.jl` documenation with optional arguments

Arguments
* `--exclude-tutorials` - exclude the tutorials from the menu of Documenter,
  this can be used if you do not have Quarto installed to still be able to render the docs
  locally on this machine. This option should not be set on CI.
* `--help`              - print this help and exit without rendering the documentation
* `--prettyurls`        – toggle the prettyurls part to true (which is otherwise only true on CI)
* `--quarto`            – run the Quarto notebooks from the `tutorials/` folder before generating the documentation
  this has to be run locally at least once for the `tutorials/*.md` files to exist that are included in
  the documentation (see `--exclude-tutorials`) for the alternative.
  If they are generated ones they are cached accordingly.
  Then you can spare time in the rendering by not passing this argument.
""",
    )
    exit(0)
end
#
# (a) if docs is not the current active environment, switch to it
# (from https://github.com/JuliaIO/HDF5.jl/pull/1020/) 
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(PackageSpec(; path=(@__DIR__) * "/../"))
    Pkg.resolve()
    Pkg.instantiate()
end

# (b) Did someone say render? Then we render!
if "--quarto" ∈ ARGS
    using CondaPkg
    CondaPkg.withenv() do
        @info "Rendering Quarto"
        tutorials_folder = (@__DIR__) * "/../tutorials"
        # instantiate the tutorials environment if necessary
        Pkg.activate(tutorials_folder)
        Pkg.resolve()
        Pkg.instantiate()
        Pkg.build("IJulia") # build IJulia to the right version.
        Pkg.activate(@__DIR__) # but return to the docs one before
        run(`quarto render $(tutorials_folder)`)
    end
end

tutorials_in_menu = true
if "--exclude-tutorials" ∈ ARGS
    @warn """
    You are excluding the tutorials from the Menu,
    which might be done if you can not render them locally.

    Remember that this should never be done on CI for the full documentation.
    """
    tutorials_in_menu = false
end

# (c) load necessary packages for the docs
using BiochemicalAlgorithms
using Documenter

DocMeta.setdocmeta!(BiochemicalAlgorithms, :DocTestSetup, :(using BiochemicalAlgorithms); recursive=true)

## Build tutorials menu
tutorials_menu =
    "How to..." => [
        "get started" => "tutorials/getting_started.md",
        "iterate" => "tutorials/iterate.md",
        "read and write" => "tutorials/read_and_write.md",
        "handle molecules" =>"tutorials/handle_molecules.md"
    ]

const pages = Any[
    "Home" => "index.md",
    (tutorials_in_menu ? [tutorials_menu] : [])...,
    "Library" => Any[
        "System representation" => "public/system.md",
        "Force fields" => "public/forcefields.md",
        "Mappings" => "public/mappings.md",
        "Internals" => Any[
            "System representation" => "private/system.md"
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

#back to main env
Pkg.activate()
