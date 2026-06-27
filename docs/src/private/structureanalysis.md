# Structure analysis
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["structureanalysis.md"]
```

## Molecular surfaces

These types and functions are internal to the molecular-surface pipeline. They
are not exported, but remain documented and callable via the qualified name
(e.g. `BiochemicalAlgorithms.make_watertight`).

### Reduced-surface graph records
```@docs
RSVertex
RSEdge
RSFace
```

### Solvent-accessible-surface graph records
```@docs
SASVertex
SASEdge
SASFace
```

### Solvent-excluded-surface graph records
```@docs
SESVertex
SESEdge
SESFace
```

### Mesh repair helpers
```@docs
weld_close_vertices
extract_manifold
fill_holes
make_watertight
```

### SES cleanup sub-steps
```@docs
split_spheric_faces!
resolve_probe_intersections!
```
