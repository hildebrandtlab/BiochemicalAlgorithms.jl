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
