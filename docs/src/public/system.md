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
atoms_df
eachatom
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
bonds_df
eachbond
nbonds
Base.push!(::System{T}, ::BondTuple) where T
```

## Molecules
```@docs
Molecule
MoleculeTuple
eachmolecule
molecule_by_idx
molecules
molecules_df
nmolecules
parent_molecule
```

## Chains
```@docs
Chain
ChainTuple
chain_by_idx
chains
chains_df
eachchain
nchains
parent_chain
Base.push!(::Molecule, ::ChainTuple)
```

## Fragments
```@docs
Fragment
FragmentTuple
eachfragment
fragment_by_idx
fragments
fragments_df
nfragments
parent_fragment
Base.push!(::Chain, ::FragmentTuple)
```

## Nucleotides
```@docs
Nucleotide
NucleotideTuple
eachnucleotide
nnucleotides
nucleotide_by_idx
nucleotides
nucleotides_df
parent_nucleotide
```

## Residues
```@docs
Residue
ResidueTuple
eachresidue
nresidues
parent_residue
residue_by_idx
residues
residues_df
Base.push!(::Chain, ::ResidueTuple)
```
