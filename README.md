# BiochemicalAlgorithms

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl/dev)
[![Build Status](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/actions/workflows/CI.yml/badge.svg?branch=develop)](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/actions/workflows/CI.yml?query=branch%3Adevelop)


BiochemicalAlgorithms.jl is a redesign of the popular [Biochemical Algorithms Library (BALL)](https://github.com/BALL-Project/ball), the largest open source C++-framework of its kind. We focused on three main design goals: efficiency, ease of use and rapid application development (RAD). Our library provides functionality for file I/O, molecular modeling, molecular mechanics methods, and molecular visualization, and hence can serve as a foundation for developing applications within the Julia ecosystem.

# Installation
To install BiochemicalAlgorithms, open a Julia REPL, switch to the package mode by pressing `]`, and type


```julia
pkg> add BiochemicalAlgorithms
``` 

# Usage

Here is a simple impression of what you can do with BiochemicalAlgorithms.jl. 
Central to every application is a `System`, which is filled with structures by reading in atom coordinates from PDB or PubChem JSON files. The system is preprocessed by the `FragmentDB` performing three steps: name normalization, reconstruction of missing atoms, and the construction of atomic bonds. The energy of the structure is evaluated using Amber forcefield.  With the help of [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl) the structure can be visualized as a ball-and-stick model.

```julia
using BiochemicalAlgorithms

# Read PDB file from the BiochemicalAlgorithms.jl repository
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

println("Number of atoms: ", natoms(sys))
println("Number of bonds: ", nbonds(sys))

# Prepare molecule
fdb = FragmentDB()
normalize_names!(sys, fdb)
reconstruct_fragments!(sys, fdb)
build_bonds!(sys, fdb)
println("Number of bonds: ", nbonds(sys))

# Create Amber force field and compute the energy of the system
amber = AmberFF(sys)
compute_energy!(amber)
println(amber.energy)
```


# Documentation

If the previous section whetted your appetite, have a look at our [tutorials](https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl/stable/tutorials/getting_started/) to get started.


# Contributing

You have ideas for improvements, criticism, or ran into problems?  You are looking for a feature that you know from BALL? 
Feedback and contributions are very welcome. Check out our [guidelines](CONTRIBUTING.md) and use our [issue tracker](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/issues) or contact us [via e-mail](mailto:hildebrandtlab@uni-mainz.de?subject=BiochemicalAlgorithms.jl).
