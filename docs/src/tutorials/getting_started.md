# Welcome to BiochemicalAlgorithms.jl


In this tutorial, you will learn about the basic concepts of the `BiochemicalAlgorithms.jl` library, a complete rewrite of the C++ framework [`BALL`](https://github.com/ball-project/ball) in Julia.

To use `BiochemicalAlgorithms.jl` in your code, add it to your project

``` julia
using Pkg
Pkg.add("BiochemicalAlgorithms")
```

and use it in your code:

``` julia
using BiochemicalAlgorithms
```

## Representing molecular systems

In BiochemicalAlgorithms.jl, all molecules are stored in so-called `Systems`. While a `System` can be created from scratch and filled programmatically, it is commonly created by reading a molecular file, such as a `PDB`-file:

``` julia
s = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
```

    System with 2241 atoms (2ptc.pdb)

You can then run methods on this system, e.g.

``` julia
println("The system $(s.name) contains $(natoms(s)) atoms.")
```

    The system 2ptc.pdb contains 2241 atoms.

### Common preparation steps

The data stored in many molecular file formats is incomplete, or needs to be normalized in certain ways. `PDB`-files, for instance, often omit hydrogen atoms, and don’t usually store bonds that can be inferred otherwise. `BiochemicalAlgorithms.jl` offers a number of methods that perform preparation steps that are common to most molecular modelling applications, such as normalizing atom- and fragment names, computing bonds, adding missing atoms from a library of templates (such as amino acids), or saturating a molecule with hydrogen atoms. The template information used by these methods is stored in the so-called `FragmentDB`. A common series of operations to prepare a system for further processing is as follows:

``` julia
fdb = FragmentDB()
normalize_names!(s, fdb)
reconstruct_fragments!(s, fdb)
build_bonds!(s, fdb)
```

    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for  CA:462
    └ @ BiochemicalAlgorithms ~/git/ball.jl/src/preprocessing/reconstruct_fragments.jl:177
    [ Info: reconstruct_fragments!(): added 2364 atoms.
    ┌ Warning: build_bonds!(): could not find reference fragment for  CA.
    └ @ BiochemicalAlgorithms ~/git/ball.jl/src/preprocessing/build_bonds.jl:14
    [ Info: build_bonds!(): built 4498 bonds
