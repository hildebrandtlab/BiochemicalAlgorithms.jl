# System representation
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["system.md"]
```

## Abstract types
```@docs
AbstractAtomContainer
AbstractMolecule
```

## Systems
```@docs
System
SystemDataFrame
default_system
Base.parent(::System)
parent_system
```

## Atoms
```@docs
Atom
AtomTuple
atoms
atoms_df
eachatom
natoms
Base.push!(::System{T}, ::AtomTuple{T}) where T
```

## Bonds
```@docs
Bond
BondTuple
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
molecules
molecules_df
nmolecules
parent_molecule
```

## Chains
```@docs
Chain
ChainTuple
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
residues
residues_df
Base.push!(::Chain, ::ResidueTuple)
```
