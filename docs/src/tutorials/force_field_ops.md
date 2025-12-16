# Working with force fields


## Preparation

Force fields in BiochemicalAlgorithms.jl are always set up for a given molecular system. Hence, we first need to prepare a `System` object [as usual](getting_started.md), e.g., by loading a structure from a PDB file and applying our fragment database preprocessors:

``` julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

fdb = FragmentDB()
normalize_names!(sys, fdb)
reconstruct_fragments!(sys, fdb)
build_bonds!(sys, fdb)

sys
```

    [ Info: reconstruct_fragments!(): added 0 atoms.
    [ Info: build_bonds!(): built 22 bonds

    System{Float32}: AA
      23 atoms
      22 bonds
       1 molecules
       1 chains
       1 secondary structures
       2 fragments

## Setting up a force field

Setting up a force field boils down to a simple constructor call:

``` julia
ff = AmberFF(sys)
```

    ┌ Warning: 2 warnings occurred during setup that were suppressed:
    │  - Torsion: 2 warnings
    │ Use print_warnings(ff) to display them.
    └ @ BiochemicalAlgorithms ~/research/BiochemicalAlgorithms.jl/src/forcefields/common/forcefield.jl:231

    AmberFF for 23 atoms with 22 bonds.

This will initialize the force field object using the atoms (and their positions in particular) as seen at time of construction. Later alterations of the same system must be made known to the force field by calling the `update!` function:

``` julia
update!(ff)
```

## Force computation

Force vectors for all atoms can be computed using the `compute_forces!` function:

``` julia
compute_forces!(ff)
```

The force vectors computed here are directly stored in the underlying system:

``` julia
atoms(sys).F   # or, equivalently, atoms(ff.system).F
```

    23-element SystemComponentTableCol{StaticArraysCore.SVector{3, Float32}}:
     [-63.913353, -14.859964, 12.627003]
     [45.423416, 8.423723, -5.2788324]
     [82.277824, 115.931015, -197.63892]
     [3812.7246, 2396.3618, -8713.229]
     [0.41454452, 1.2887306, 0.053740263]
     [-0.24682106, -4.6684403, 6.2538323]
     [2.2589614, 0.26905724, 0.4852221]
     [1.5536233, -0.3578552, 0.102772325]
     [-14.104118, 3.8311667, -20.902191]
     [1.2480733, 0.6318529, -2.2023916]
     ⋮
     [-41.699554, -4.197597, 120.749504]
     [23.966557, -25.913675, 20.226442]
     [0.84398764, -0.20457195, -2.9093812]
     [0.71285784, 1.493221, 0.039408203]
     [1.5202898, -0.209776, 0.41397908]
     [-14.851994, 14.544389, -15.180875]
     [-0.56523114, 8.413886, -0.5651454]
     [-3795.7595, -2608.7483, 8814.81]
     [2.5838406, -12.319646, 4.423078]

The force vectors are given in units of kJ/(mol·). In previous versions, forces were computed in Newton.

## Potential energy computation

Similarly, potential energies can be computed via `compute_energy!`:

``` julia
compute_energy!(ff)
```

    1425.5991f0

This will return the total energy of the system (in kJ/mol). Additionally, the force field object keeps track of individual contributions of force field components, which can be queried like this:

``` julia
ff.energy
```

    Dict{String, Float32} with 7 entries:
      "Improper Torsion" => 3.99017f-6
      "Hydrogen Bonds"   => 0.0
      "Proper Torsion"   => 10.7981
      "Electrostatic"    => -85.1466
      "Bond Stretches"   => 1.36306
      "Angle Bends"      => 5.40767
      "Van der Waals"    => 1493.18

## Structure optimization

The force field object can also be used to find a structure minimizing the total energy of the system:

``` julia
optimize_structure!(ff)
compute_energy!(ff)
```

    -133.99623f0

The opmitization function also updates the atom positions correspondingly such that the structure can be visualized in its optimized state, e.g., through [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl).

We provide a variant of the optimization function above to only optimize the positions of hydrogen atoms:

``` julia
optimize_hydrogen_positions!(ff)
```
