# Biomolecular systems
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["system.md"]
```

## Systems
```@docs
System
default_system
parent_system
Base.parent(::System)
Base.empty!(::System)
```

## System components
```@docs
AbstractAtomContainer
AbstractColumnTable
AbstractSystemComponent
AbstractSystemComponentTable
SystemComponentTable
SystemComponentTableCol
get_property
has_flag
has_property
set_flag!
set_property!
unset_flag!
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
MoleculeVariant
MoleculeVariantType
```

### Molecules (all variants)
```@docs
Molecule
MoleculeTable
molecule_by_idx
molecules
nmolecules
parent_molecule
Base.delete!(::Molecule)
Base.push!(::System{T}, ::Molecule{T}) where T
```

### Proteins (`MoleculeVariant.Protein`)
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
FragmentVariant
FragmentVariantType
```

### Fragments (all variants)
```@docs
Fragment
FragmentTable
fragment_by_idx
fragments
nfragments
parent_fragment
Base.delete!(::Fragment)
Base.push!(::Chain{T}, ::Fragment{T}) where T
```

### Nucleotides (`FragmentVariant.Nucleotide`)
```@docs
Nucleotide
isnucleotide
nnucleotides
nucleotide_by_idx
nucleotides
parent_nucleotide
```

### Nucleotides (`FragmentVariant.Residue`)
```@docs
Residue
isresidue
nresidues
parent_residue
residue_by_idx
residues
```
