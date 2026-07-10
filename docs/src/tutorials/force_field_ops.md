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
     [-53.629234, -22.412619, 24.928682]
     [46.985428, 3.0840526, -2.7236192]
     [128.1008, 45.797974, -103.76107]
     [3828.1062, 2481.3315, -8883.52]
     [-3.9139192, -5.606783, 3.236692]
     [2.53753, -4.523301, 6.5197163]
     [0.52676713, -4.1225967, 4.9825106]
     [4.2090273, 0.87796724, -1.8575672]
     [-11.800568, -6.840048, -6.6562223]
     [4.375335, 2.908515, -2.6804132]
     ⋮
     [14.85136, -70.92931, 156.84335]
     [49.850376, -32.66954, 52.788876]
     [2.0185974, -8.8306465, 1.9078683]
     [2.5104997, 0.34557226, 6.3368053]
     [4.279133, -4.6336436, 8.192726]
     [-27.39077, 46.13508, -41.05185]
     [0.59128666, -16.231728, 13.140886]
     [-3896.7554, -2587.3542, 8858.033]
     [-12.687141, 7.244336, -14.078792]

The force vectors are given in units of kJ/(mol·Å). Before BiochemicalAlgorithms.jl v0.6, forces were computed in Newton.

## Potential energy computation

Similarly, potential energies can be computed via `compute_energy!`:

``` julia
compute_energy!(ff)
```

    1425.599f0

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

    -374.3113f0

The opmitization function also updates the atom positions correspondingly such that the structure can be visualized in its optimized state, e.g., through [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl).

We provide a variant of the optimization function above to only optimize the positions of hydrogen atoms:

``` julia
optimize_hydrogen_positions!(ff)
```
