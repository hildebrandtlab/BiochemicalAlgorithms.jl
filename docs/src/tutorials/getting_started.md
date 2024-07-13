# Welcome to BiochemicalAlgorithms.jl

In this tutorial, you will learn about the basic concepts of the
`BiochemicalAlgorithms.jl`-library (sometimes shortened to `BALL.jl`) –
a complete rewrite of the C++-framework `BALL` in Julia.

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

In BiochemicalAlgorithms.jl, all molecules are stored in so-called
`Systems`. While a `System` can be created from scratch and filled
programmatically, it is commonly created by reading a molecular file,
such as a `PDB`-file:

``` julia
s = load_pdb("data/5PTI.pdb")
```

    System with 1087 atoms (5PTI.pdb)

You can then run methods on this system, e.g.

``` julia
println("The system $(s.name) contains $(natoms(s)) atoms.")
```

    The system 5PTI.pdb contains 1087 atoms.

### Common preparation steps

The data stored in many molecular file formats is incomplete, or needs
to be normalized in certain ways. `PDB`-files, for instance, often omit
hydrogen atoms, and don’t usually store bonds that can be inferred
otherwise. `BiochemicalAlgorithms.jl` offers a number of methods that
perform preparation steps that are common to most molecular modelling
applications, such as normalizing atom- and fragment names, computing
bonds, adding missing atoms from a library of templates (such as amino
acids), or saturating a molecule with hydrogen atoms. The template
information used by these methods is stored in the so-called
`FragmentDB`. A common series of operations to prepare a system for
further processing is as follows:

``` julia
fdb = FragmentDB()

normalize_names!(s, fdb)
reconstruct_fragments!(s, fdb)
build_bonds!(s, fdb)
```
    [ Info: reconstruct_fragments!(): added 109 atoms.
    [ Info: build_bonds!(): built 912 bonds
