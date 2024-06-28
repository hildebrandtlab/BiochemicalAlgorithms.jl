# BiochemicalAlgorithms

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl/stable)-->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl/dev)
[![Build Status](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/actions/workflows/CI.yml/badge.svg?branch=develop)](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/actions/workflows/CI.yml?query=branch%3Adevelop)
<!--[![Coverage](https://codecov.io/gh/hildebrandtlab/BiochemicalAlgorithms.jl/branch/develop/graph/badge.svg)](https://codecov.io/gh/hildebrandtlab/BiochemicalAlgorithms.jl)-->


BiochemicalAlgorithms.jl is a redesign of the popular [Biochemical Algorithms Library (BALL)](https://github.com/BALL-Project/ball), the largest open source C++-framework of its kind. We focused on three main design goals: efficiency, ease of use and rapid application development (RAD). Our library provides functionality for file I/O, molecular modelling, molecular mechanics methods, and molecular visualization, and hence can serve as a foundation for developing applications within the Julia ecosystem.

# Installation
BiochemicalAlgorithms is not yet registered as a Julia Language package. To install BiochemicalAlgorithms open a Julia REPL, switch to the package mode by pressing `]` and type


```{julia}
pkg> add https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl.git
``` 

# Usage

Here is a simple impression of what you can do with BiochemicalAlgorithms.jl. 
Central to every application is a `System`, which is filled with structures by reading in atom coordinates from pdb or pubchem json files. The system is preprocessed by the `FragementDB` performing three steps: name normalization, reconstruction of missing atoms and the construction of atomic bonds. The energy of the structure is evaluated using Amber forcefield.  With the help of [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl) the structure can be visualized as a ball_and_stick model.

```{julia}
using BiochemicalAlgorithms
using BiochemicalVisualization

# Read PDB file from the BiochemicalAlgorithms.jl repository
pdb_sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

println("Number of atoms: ", natoms(sys))
println("Number of bonds: ", nbonds(sys))

# Prepare Molecule
fdb = FragmentDB()
normalize_names!(sys, fdb)
reconstruct_fragments!(sys, fdb)
build_bonds!(sys, fdb)
println("Number of bonds: ", nbonds(sys))

# Create Amber ForceField and compute the energy of the system
amber = AmberForceField(sys)
compute_energy(amber)
a_ff.energy["Bond Stretches"] 

ball_and_stick(sys)
``` 


# Documentation

If the previous section whetted your appetite, have a look at our [tutorials](/tutorials/) to get started or look at our [documentation](https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl/dev/) for more details.


# Contributing

You have ideas for improvements, critics or ran into problems?  You are looking for a feature that you know of BALL? 
Feedback and contributions are very welcome. Check out our [guidelines](CONTRIBUTING.md) and contact us via mail [Hildebrandt Group](mailto:hildebrandtlab@uni-mainz.de?subject=BiochemicalAlgorithms.jl)


# Cite

Join our talk on [Juliacon 2024](https://pretalx.com/juliacon2024/talk/DFNXCK/). Hope to see you there! 

Or have a look at our poster presented at [ISMB/ ECCB 2023](https://iscb.junolive.co/ismb2023/live/exhibitor/ismbeccb2023_poster_1989) 
