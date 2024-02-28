# System representation
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["system.md"]
```

## Abstract types
```@docs
AbstractSystemComponent
AbstractAtomContainer
AbstractMolecule
```

## Common functions
```@docs
has_property
get_property
set_property!
has_flag
set_flag!
unset_flag!
```

## Systems
```@docs
System
default_system
Base.parent(::System)
parent_system
```

## Atoms
```@docs
Atom
AtomTuple
atom_by_idx
atom_by_name
atoms
is_bound_to
is_geminal
is_vicinal
natoms
Base.push!(::System{T}, ::AtomTuple{T}) where T
```

## Bonds
```@docs
Bond
BondTuple
bond_by_idx
bonds
nbonds
Base.push!(::System{T}, ::BondTuple) where T
```

## Molecules
```@docs
Molecule
MoleculeTuple
molecule_by_idx
molecules
nmolecules
parent_molecule
```

## Chains
```@docs
Chain
ChainTuple
chain_by_idx
chains
nchains
parent_chain
Base.push!(::Molecule, ::ChainTuple)
```

## Fragments
```@docs
Fragment
FragmentTuple
fragment_by_idx
fragments
nfragments
parent_fragment
Base.push!(::Chain, ::FragmentTuple)
```

## Nucleotides
```@docs
Nucleotide
NucleotideTuple
nnucleotides
nucleotide_by_idx
nucleotides
parent_nucleotide
```

## Residues
```@docs
Residue
ResidueTuple
nresidues
parent_residue
residue_by_idx
residues
Base.push!(::Chain, ::ResidueTuple)
```
