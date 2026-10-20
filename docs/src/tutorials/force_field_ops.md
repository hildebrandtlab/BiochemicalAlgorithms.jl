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
     [-53.629223, -22.41262, 24.92869]
     [46.976604, 3.0782008, -2.7205548]
     [128.10509, 45.8008, -103.76251]
     [3828.1155, 2481.3374, -8883.526]
     [-3.913918, -5.6067834, 3.2366881]
     [2.5375302, -4.523299, 6.5197215]
     [0.5267689, -4.1225977, 4.9825115]
     [4.2090273, 0.8779677, -1.8575673]
     [-11.800571, -6.840049, -6.6562233]
     [4.375358, 2.9085383, -2.680446]
     ⋮
     [14.851368, -70.92931, 156.84334]
     [49.868458, -32.65759, 52.782658]
     [2.0185935, -8.830651, 1.9078697]
     [2.5105174, 0.3455763, 6.336817]
     [4.279133, -4.6336427, 8.192726]
     [-27.39079, 46.13508, -41.05186]
     [0.59128666, -16.231726, 13.140888]
     [-3896.7551, -2587.354, 8858.032]
     [-12.687141, 7.24434, -14.078795]

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

    -374.31235f0

The opmitization function also updates the atom positions correspondingly such that the structure can be visualized in its optimized state, e.g., through [BiochemicalVisualization.jl](https://github.com/hildebrandtlab/BiochemicalVisualization.jl).

We provide a variant of the optimization function above to only optimize the positions of hydrogen atoms:

``` julia
optimize_hydrogen_positions!(ff)
```

## Compiled force fields

`AmberFF` is the *reference* implementation: it is easy to inspect and always consistent with the underlying `System`. Its components refer to `Atom`s of the system, and every call to `update!` rebuilds the nonbonded interaction list from scratch. This is fine for single evaluations. Workloads that evaluate the same force field many times for slightly different coordinates are better served by a `CompiledForceField`, which

- flattens all interactions into contiguous arrays that refer to atoms by row number,
- evaluates energies and forces on private coordinate and force buffers instead of the atom table,
- rebuilds the nonbonded pair list only when needed, i.e., when atoms have moved far enough (a *Verlet skin*), at a fixed evaluation interval, or on explicit request, and
- optionally uses multiple threads (or a GPU) and a different accumulation precision.

It computes the same energy terms and forces as the reference implementation. It is purely an evaluator: parameter assignment and topology handling are still done by `AmberFF`.

### Compiling a force field

Compiling an existing force field is a single call:

``` julia
cff = compile(ff)
```

    CompiledForceField{Float32, acc=Float32, SerialBackend}: 23 atoms, 22 stretches, 39 bends, 49+1 torsions, 192 vdW / 0 hbond / 192 electrostatic pairs

`compile` calls `update!` on `ff` once and lowers the resulting interactions, so it reflects the atom positions at the time of the call.

!!! warning
The atom topology (and with it the numbering of atoms) must not change during the lifetime of a compiled force field. If atoms or bonds are added or removed, create a new `AmberFF` and compile it again.

Energies and forces are computed just like for `AmberFF`:

``` julia
compute_energy!(cff)
```

    -374.31317f0

``` julia
compute_forces!(cff)
```

In contrast to `AmberFF`, results are **not** written to the system. Forces live in the compiled force field’s own buffer, and the energy contributions of the individual components are kept in a dictionary:

``` julia
cff.energy
```

    Dict{String, Float32} with 7 entries:
      "Angle Bends"      => 29.8607
      "Hydrogen Bonds"   => 0.0
      "Bond Stretches"   => 1.82887
      "Van der Waals"    => 31.7756
      "Improper Torsion" => 4.46057
      "Electrostatic"    => -487.006
      "Proper Torsion"   => 44.7673

``` julia
cff.F
```

    23-element Vector{StaticArraysCore.SVector{3, Float32}}:
     [0.4492079, 0.13000447, -0.5023289]
     [-0.07742739, -0.41107702, 0.028604984]
     [-0.12038422, 0.08013725, -0.3902855]
     [-0.28648376, -0.16231918, 0.33292007]
     [-0.03247261, -0.059452057, 0.032169342]
     [-0.031589985, 0.03486401, 0.121210754]
     [-0.028962612, 0.012744546, -0.00989151]
     [0.06625263, 0.0018118024, -0.025928617]
     [-0.01142025, 0.0013866425, 0.0315094]
     [-0.042267933, -0.02734971, -0.023148417]
     ⋮
     [0.41654098, -0.824533, 0.38330555]
     [-0.035525322, -0.067474365, -0.010866165]
     [-0.021340847, -0.046815872, -0.017537594]
     [0.0464859, 0.084653854, -0.07345843]
     [0.009778619, 0.0300107, 0.033632398]
     [0.20737839, 0.1842537, 0.23509502]
     [-0.032557487, 0.0028123856, -0.037528038]
     [-0.13605976, 0.3531927, -0.41021776]
     [-0.00066280365, 0.47307587, -0.020980358]

Both implementations agree on the total energy:

``` julia
compute_energy!(ff) ≈ compute_energy!(cff)
```

    true

### Synchronizing with the system

The atom table is only updated on request:

``` julia
BiochemicalAlgorithms.sync_to_table!(cff)

atoms(sys).F
```

    23-element SystemComponentTableCol{StaticArraysCore.SVector{3, Float32}}:
     [0.4492079, 0.13000447, -0.5023289]
     [-0.07742739, -0.41107702, 0.028604984]
     [-0.12038422, 0.08013725, -0.3902855]
     [-0.28648376, -0.16231918, 0.33292007]
     [-0.03247261, -0.059452057, 0.032169342]
     [-0.031589985, 0.03486401, 0.121210754]
     [-0.028962612, 0.012744546, -0.00989151]
     [0.06625263, 0.0018118024, -0.025928617]
     [-0.01142025, 0.0013866425, 0.0315094]
     [-0.042267933, -0.02734971, -0.023148417]
     ⋮
     [0.41654098, -0.824533, 0.38330555]
     [-0.035525322, -0.067474365, -0.010866165]
     [-0.021340847, -0.046815872, -0.017537594]
     [0.0464859, 0.084653854, -0.07345843]
     [0.009778619, 0.0300107, 0.033632398]
     [0.20737839, 0.1842537, 0.23509502]
     [-0.032557487, 0.0028123856, -0.037528038]
     [-0.13605976, 0.3531927, -0.41021776]
     [-0.00066280365, 0.47307587, -0.020980358]

Use `sync_to_table!(cff; forces=false)` to copy only the coordinates. Conversely, changes made to atom positions in the system are *not* picked up automatically; positions are read when the force field is compiled, so recompile after modifying the system by hand.

### Evaluation options

`compile` accepts a few keyword arguments:

``` julia
# use all available threads for energy and force evaluation
cff_threads = compile(ff; backend=:threads)

# accumulate energies and forces in double precision
cff_64 = compile(ff; acc=Float64)
```

- `backend`: `:serial` (default), `:threads`, or `:gpu`. The GPU backend requires `using KernelAbstractions` and a device package such as `Metal` or `CUDA`.
- `acc`: the floating point type used for accumulation (`Float32` or `Float64`).
- `max_displacement`: size of the Verlet skin in Å. The nonbonded pair list is rebuilt once an atom has moved further than this from its position at the last build. By default, a safe value is derived from the force field’s cutoffs.
- `update_frequency`: additionally forces a pair list rebuild every given number of evaluations (`0`, the default, disables this).

A rebuild can also be triggered manually at any time, e.g., after changing `cff.r` yourself:

``` julia
rebuild_pairlist!(cff)
```

### Structure optimization

`optimize_structure!` uses a compiled force field internally, and forwards `backend`, `max_displacement`, and `update_frequency` to `compile` (and `accumulation` as `acc`):

``` julia
optimize_structure!(ff; accumulation=Float64)
```

If you need more control, e.g., for your own optimizer or a dynamics loop, drive a compiled force field directly. `prepare_eval!` loads a flat coordinate vector (x₁, y₁, z₁, x₂, …) and refreshes the pair list if needed:

``` julia
r = collect(Float64, Iterators.flatten(atoms(sys).r))

BiochemicalAlgorithms.prepare_eval!(cff, r)
compute_energy!(cff)
```

    -374.31525f0

### When to use which

|  | `AmberFF` | `CompiledForceField` |
|----|----|----|
| Setup | `AmberFF(sys)` | `compile(ff)` on top of an `AmberFF` |
| Results | written to the system | kept in own buffers; `sync_to_table!` on demand |
| Changes to the system | `update!(ff)` | recompile |
| Nonbonded pair list | rebuilt by `update!` | rebuilt when atoms leave the Verlet skin, every `update_frequency` evaluations if set, or manually via `rebuild_pairlist!(cff)` |
| Parallelism, precision | serial, system precision | threads/GPU, selectable accumulation type |
| Best for | single evaluations, inspecting individual interactions | repeated evaluations: optimization, sampling |
