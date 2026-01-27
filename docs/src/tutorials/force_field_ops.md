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
     [-53.629234, -22.412619, 24.92868]
     [46.985428, 3.0840528, -2.7236197]
     [128.10083, 45.797966, -103.76105]
     [3828.1062, 2481.3315, -8883.5205]
     [-3.9139185, -5.606783, 3.2366903]
     [2.5375302, -4.523301, 6.5197177]
     [0.5267676, -4.122597, 4.9825087]
     [4.2090273, 0.8779674, -1.8575675]
     [-11.800565, -6.840048, -6.6562223]
     [4.3753357, 2.9085145, -2.680413]
     ⋮
     [14.851362, -70.929306, 156.84334]
     [49.850376, -32.66954, 52.78887]
     [2.0185971, -8.830647, 1.9078681]
     [2.5104997, 0.3455725, 6.3368053]
     [4.2791324, -4.633644, 8.192726]
     [-27.39077, 46.135075, -41.05185]
     [0.59128755, -16.231728, 13.140887]
     [-3896.7554, -2587.3545, 8858.033]
     [-12.687142, 7.2443366, -14.078793]

The force vectors are given in units of kJ/(mol·Å). Before BiochemicalAlgorithms.jl v0.6, forces were computed in Newton.

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

    -374.31235f0

The opmitization function also updates the atom positions correspondingly such that the structure can be visualized in its optimized state, e.g., through [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl).

We provide a variant of the optimization function above to only optimize the positions of hydrogen atoms:

``` julia
optimize_hydrogen_positions!(ff)
```
