# Welcome to BiochemicalAlgorithms.jl


# Welcome to BiochemicalAlgorithms.jl

In this tutorial, you will learn about the basic concepts of the
`BiochemicalAlgorithms.jl`-library (sometimes shortened to `BALL.jl`) –
a complete rewrite of the C++-framework `BALL` in Julia.

To use `BiochemicalAlgorithms.jl` in your code, add it to your project

``` julia
import Pkg
#Pkg.add(url="https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl")
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

    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for PO4:70
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:80
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:101
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:102
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:105
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:110
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:111
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:112
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:113
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:116
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:117
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:119
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:121
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:122
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:125
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:126
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:127
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:129
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:133
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:134
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:138
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:140
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:143
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:144
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:145
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:146
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:156
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:157
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:158
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:159
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:160
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:200
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:201
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:202
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:203
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:204
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:205
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:209
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:210
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:211
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:212
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:214
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:216
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:217
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:218
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:219
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:220
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:223
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:225
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:302
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:304
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:310
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:311
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:312
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:313
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:314
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:315
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:316
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:317
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:318
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:319
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:320
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:321
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:322
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for UNX:324
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\reconstruct_fragments.jl:177
    [ Info: reconstruct_fragments!(): added 109 atoms.
    ┌ Warning: build_bonds!(): could not find reference fragment for PO4.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for UNX.
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\preprocessing\build_bonds.jl:14
    [ Info: build_bonds!(): built 912 bonds
