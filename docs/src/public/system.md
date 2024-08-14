# System representation
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["system.md"]
```

## Abstract types
```@docs
AbstractColumnTable
AbstractSystemComponentTable
AbstractSystemComponent
AbstractAtomContainer
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
parent_system
Base.parent(::System)
Base.empty!(::System)
```

## Atoms
```@docs
Atom
AtomTable
atom_by_idx
atom_by_name
atoms
is_bound_to
is_geminal
is_vicinal
natoms
Base.delete!(::Atom)
Base.push!(::System{T}, ::Atom{T}) where T
```

## Bonds
```@docs
Bond
BondTable
bond_by_idx
bonds
nbonds
Base.delete!(::Bond)
Base.push!(::System{T}, ::Bond{T}) where T
```

## Molecules
```@docs
Molecule
MoleculeTable
MoleculeVariant
MoleculeVariantType
molecule_by_idx
molecules
nmolecules
parent_molecule
Base.delete!(::Molecule)
Base.push!(::System{T}, ::Molecule{T}) where T
```

### Proteins (molecule variant)
```@docs
Protein
isprotein
nproteins
parent_protein
protein_by_idx
proteins
```

## Chains
```@docs
Chain
ChainTable
chain_by_idx
chains
nchains
parent_chain
Base.delete!(::Chain)
Base.push!(::Molecule{T}, ::Chain{T}) where T
```

## Fragments
```@docs
Fragment
FragmentTable
FragmentVariant
FragmentVariantType
fragment_by_idx
fragments
nfragments
parent_fragment
Base.delete!(::Fragment)
Base.push!(::Chain{T}, ::Fragment{T}) where T
```

### Nucleotides (fragment variant)
```@docs
Nucleotide
isnucleotide
nnucleotides
nucleotide_by_idx
nucleotides
parent_nucleotide
```

### Residues (fragment variant)
```@docs
Residue
isresidue
nresidues
parent_residue
residue_by_idx
residues
```
