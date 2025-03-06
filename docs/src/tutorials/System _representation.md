# Getting started


BiochemicalAlgorithms.jl\`s design is inspired by its predecessor [`BALL`](https://github.com/ball-project/ball). BALL and BiochemicalAlgorithms.jl both have a set of KERNEL classes or core classes to represent molecular entities such as Molecules, Proteins, Atoms, Bonds etc. The following class diagramm shows the different classes and their relationship to each others.

![UML class diagram *Note*: Only the most important functionalities are shown here.](uml.png)

UML class diagram *Note*: Only the most important functionalities are shown here.

In the center of all classes resides the system. So let’s see it in action with a simple peptide:

``` julia
using BiochemicalAlgorithms

s = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
s
```

    System{Float32}: AlaAla.pdb
      23 atoms
       0 bonds
       1 molecules
       1 chains
       0 secondary structures
       2 fragments
