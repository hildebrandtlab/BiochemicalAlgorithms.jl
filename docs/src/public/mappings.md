# Mappings
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["mappings.md"]
```

## Atom bijection
```@docs
AbstractAtomBijection
TrivialAtomBijection
atoms(::AbstractAtomBijection)
```

## Rigid Mapping
```@docs
AbstractRMSDMinimizer
RMSDMinimizerCoutsias
RMSDMinimizerKabsch
RigidTransform
compute_rmsd
compute_rmsd_minimizer
map_rigid!
match_points
rigid_transform!
translate!
```
