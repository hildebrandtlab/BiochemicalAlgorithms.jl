# Working with force fields


## Preparation

Force fields in BiochemicalAlgorithms.jl are always set up for a given molecular system. Hence, we first need to prepare a `System` object [as usual](getting_started.md), e.g., by loading a structure from a PDB file and applying our fragment database preprocessors:

``` julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

infer_topology!(sys)

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
    └ @ BiochemicalAlgorithms ~/git/ball.jl/src/forcefields/common/forcefield.jl:231

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
     [-53.62923, -22.412622, 24.928682]
     [46.985413, 3.084029, -2.7235892]
     [128.1008, 45.79798, -103.76105]
     [3828.1062, 2481.3313, -8883.522]
     [-3.913918, -5.6067834, 3.236691]
     [2.5375295, -4.5233026, 6.519717]
     [0.52676797, -4.122597, 4.9825106]
     [4.2090273, 0.8779677, -1.8575673]
     [-11.800568, -6.8400474, -6.6562233]
     [4.3753576, 2.9085383, -2.680446]
     ⋮
     [14.851358, -70.929306, 156.84334]
     [49.850372, -32.669544, 52.788876]
     [2.0185974, -8.830647, 1.9078674]
     [2.5105, 0.34557325, 6.336806]
     [4.2791333, -4.6336436, 8.192725]
     [-27.39077, 46.13508, -41.05185]
     [0.5912861, -16.231726, 13.140887]
     [-3896.7554, -2587.3542, 8858.033]
     [-12.687143, 7.2443376, -14.0787945]

The force vectors are given in units of kJ/(mol·Å). Before BiochemicalAlgorithms.jl v0.6, forces were computed in Newton.

## Potential energy computation

Similarly, potential energies can be computed via `compute_energy!`:

``` julia
compute_energy!(ff)
```

    1425.6053f0

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

    -374.31473f0

The opmitization function also updates the atom positions correspondingly such that the structure can be visualized in its optimized state, e.g., through [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl).

We provide a variant of the optimization function above to only optimize the positions of hydrogen atoms:

``` julia
optimize_hydrogen_positions!(ff)
```
