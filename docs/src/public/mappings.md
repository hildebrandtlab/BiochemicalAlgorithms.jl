# Mappings
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["mappings.md"]
```

## Atom bijections
```@docs
AbstractAtomBijection
TrivialAtomBijection
atoms(::AbstractAtomBijection)
```

## Rigid mapping
```@docs
AbstractRMSDMinimizer
RMSDMinimizerCoutsias
RMSDMinimizerKabsch
RigidTransform
compute_rmsd
compute_rmsd_minimizer
map_rigid!
rigid_transform!
translate!
```
