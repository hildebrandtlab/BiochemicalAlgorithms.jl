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

    System{Float32}: AlaAla.pdb
      23 atoms
      22 bonds
       1 molecules
       1 chains
       0 secondary structures
       2 fragments

## Setting up a force field

Setting up a force field boils down to a simple constructor call:

``` julia
ff = AmberFF(sys)
```

    ┌ Warning: 2 warnings occurred during setup that were suppressed:
    │  - Torsion: 2 warnings
    │ Use print_warnings(ff) to display them.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/forcefields/common/forcefield.jl:233

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
     [-8.9053415f-10, -3.7217032f-10, 4.1395054f-10]
     [7.8021095f-10, 5.1211528f-11, -4.522637f-11]
     [2.127164f-9, 7.604934f-10, -1.7229925f-9]
     [6.356721f-8, 4.1203478f-8, -1.4751436f-7]
     [-6.499216f-11, -9.310281f-11, 5.3746518f-11]
     [4.2136676f-11, -7.511121f-11, 1.08262455f-10]
     [8.747193f-12, -6.845735f-11, 8.273652f-11]
     [6.989254f-11, 1.4578987f-11, -3.0845635f-11]
     [-1.959529f-10, -1.1358163f-10, -1.10529204f-10]
     [7.265452f-11, 4.82974f-11, -4.450986f-11]
     ⋮
     [2.466124f-10, -1.1778087f-9, 2.604445f-9]
     [8.2778495f-10, -5.4249055f-10, 8.7658f-10]
     [3.3519604f-11, -1.4663634f-10, 3.1680894f-11]
     [4.1687834f-11, 5.73837f-12, 1.0522514f-10]
     [7.105666f-11, -7.694347f-11, 1.360434f-10]
     [-4.5483453f-10, 7.660909f-10, -6.81682f-10]
     [9.818537f-12, -2.6953412f-10, 2.1820956f-10]
     [-6.470714f-8, -4.2964025f-8, 1.4709109f-7]
     [-2.1067495f-10, 1.2029504f-10, -2.337838f-10]

## Potential energy computation

Similarly, potential energies can be computed via `compute_energy!`:

``` julia
compute_energy!(ff)
```

    1425.6064f0

This will return the total energy of the system (in kJ/mol). Additionally, the force field object keeps track of individual contributions of force field components, which can be queried like this:

``` julia
ff.energy
```

    Dict{String, Float32} with 7 entries:
      "Angle Bends"      => 5.40767
      "Hydrogen Bonds"   => 0.0
      "Bond Stretches"   => 1.36306
      "Van der Waals"    => 1493.18
      "Improper Torsion" => 3.99017f-6
      "Electrostatic"    => -85.1465
      "Proper Torsion"   => 10.7981

## Structure optimization

The force field object can also be used to find a structure minimizing the total energy of the system:

``` julia
optimize_structure!(ff)
compute_energy!(ff)
```

    -374.3089f0

The opmitization function also updates the atom positions correspondingly such that the structure can be visualized in its optimized state, e.g., through [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl).

We provide a variant of the optimization function above to only optimize the positions of hydrogen atoms:

``` julia
optimize_hydrogen_positions!(ff)
```
