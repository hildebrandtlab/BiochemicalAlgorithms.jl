# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BiochemicalAlgorithms.jl is a Julia package redesigning the C++ Biochemical Algorithms Library (BALL). It provides molecular modeling capabilities including file I/O (PDB, mmCIF, PubChem JSON, SDF), core molecular data structures, molecular mechanics/force fields (AMBER, MMFF94), and preprocessing pipelines.

Documentation: https://hildebrandtlab.github.io/BiochemicalAlgorithms.jl

## Development Commands

```bash
# Activate and install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run all tests
julia --project -e 'using TestItemRunner; @run_package_tests'

# Run tests with verbose output (same as CI)
julia --project -e 'using TestItemRunner; @run_package_tests verbose=true'

# Run a specific test file
julia --project test/core/test_system.jl

# Build documentation
julia --project=docs docs/make.jl
```

Tests use TestItemRunner with the `@testitem` macro. Tests tagged with `:skip_ci` are excluded in CI environments.

## Architecture

### Molecular Hierarchy

```
System (root container)
  └─ Molecule(s)
      └─ Chain(s)
          └─ Fragment(s) [residues/nucleotides]
```

Each level contains Atoms and Bonds. The type parameter `T` (Float32 or Float64) controls numeric precision throughout.

### Key Abstractions

- **System{T}**: Root container with columnar table storage (`_atom_table`, `_bond_table`, etc.)
- **AbstractSystemComponent{T}**: Base trait for Atom, Bond, Molecule, Chain, Fragment
- **AtomContainer{T}**: Trait for types that hold atoms/bonds (System, Molecule, Chain, Fragment, Substructure)
- **Vector3{T}**: `SVector{3,T}` from StaticArrays for coordinates, velocities, forces

### Module Organization

- `src/core/`: Type definitions, system components, table implementations
- `src/fileformats/`: PDB, mmCIF, PubChem JSON, SDF, HIN file I/O
- `src/forcefields/`: ForceField base class and AMBER/MMFF94 implementations
- `src/preprocessing/`: FragmentDB, name normalization, bond building, hydrogen addition
- `src/substructures/`: SMARTS matching, SSSR (smallest set of smallest rings)
- `src/mappings/`: Atom bijections and rigid structure alignment
- `src/optimization/`: Structure minimization via L-BFGS-B

### Force Field Architecture

```
ForceField{T}
  ├─ parameters (from INI files in data/forcefields/)
  ├─ atom_type_templates
  └─ components: [StretchComponent, BendComponent, TorsionComponent, NonBondedComponent]
```

Workflow: `assign_typenames_and_charges!()` → `setup!()` → `compute_energy!()` / `compute_forces!()`

### Standard Preprocessing Pipeline

```julia
fdb = FragmentDB()
normalize_names!(sys, fdb)       # Standardize atom/residue names
reconstruct_fragments!(sys, fdb) # Add missing atoms from templates
build_bonds!(sys, fdb)           # Infer bonds from coordinates
add_hydrogens!(sys, fdb)         # Add hydrogen atoms
```

## Contributing

All commits must include a DCO sign-off:
```bash
git commit -s -m "Your message"
```

This adds: `Signed-off-by: Name <email>` to the commit message.

Default branch for PRs is `develop`.
