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

    System{Float32} with 23 atoms (AlaAla.pdb)

## Setting up a force field

Setting up a force field boils down to a simple constructor call:

``` julia
ff = AmberFF(sys)
```

    ┌ Warning: 2 warnings occurred during setup that were suppressed:
    │ Components:
    │ Torsion: 2 warnings
    │ Use print_warnings(ff) to display them.
    └ @ BiochemicalAlgorithms ~/git/ball.jl/src/forcefields/common/forcefield.jl:164

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
     [-8.905341f-10, -3.7217035f-10, 4.1395054f-10]
     [7.802113f-10, 5.1211924f-11, -4.5226885f-11]
     [2.1271647f-9, 7.6049345f-10, -1.7229925f-9]
     [6.35672f-8, 4.1203478f-8, -1.4751434f-7]
     [-6.499216f-11, -9.310282f-11, 5.374651f-11]
     [4.213668f-11, -7.5111216f-11, 1.0826247f-10]
     [8.747198f-12, -6.8457344f-11, 8.2736526f-11]
     [6.989254f-11, 1.4578986f-11, -3.0845635f-11]
     [-1.9595292f-10, -1.1358164f-10, -1.1052925f-10]
     [7.2654154f-11, 4.829702f-11, -4.4509316f-11]
     ⋮
     [2.4661245f-10, -1.1778086f-9, 2.604445f-9]
     [8.2778495f-10, -5.4249055f-10, 8.765799f-10]
     [3.3519593f-11, -1.4663633f-10, 3.168089f-11]
     [4.1687834f-11, 5.7383577f-12, 1.0522515f-10]
     [7.105667f-11, -7.6943465f-11, 1.3604341f-10]
     [-4.5483453f-10, 7.660909f-10, -6.81682f-10]
     [9.818538f-12, -2.6953412f-10, 2.1820955f-10]
     [-6.4707145f-8, -4.2964025f-8, 1.470911f-7]
     [-2.1067495f-10, 1.2029502f-10, -2.3378385f-10]

## Potential energy computation

Similarly, potential energies can be computed via `compute_energy!`:

``` julia
compute_energy!(ff)
```

    1425.6062f0

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
      "Electrostatic"    => -85.1466
      "Proper Torsion"   => 10.7981

## Structure optimization

The force field object can also be used to find a structure minimizing the total energy of the system:

``` julia
optimize_structure!(ff)
compute_energy!(ff)
```

    -374.31445f0

The opmitization function also updates the atom positions correspondingly such that the structure can be visualized in its optimized state, e.g., through [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl).

We provide a variant of the optimization function above to only optimize the positions of hydrogen atoms:

``` julia
optimize_hydrogen_positions!(ff)
```
