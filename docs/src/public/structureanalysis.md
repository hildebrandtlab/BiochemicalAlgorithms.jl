# Structure analysis
```@meta
CurrentModule = BiochemicalAlgorithms
```

```@index
Pages = ["structureanalysis.md"]
```

## Bounding box computation
```@docs
compute_bounding_box
```

## Molecular surfaces

### Types
```@docs
AbstractMolecularSurface
Circle3
```

### Radii and mesh helpers
```@docs
vdw_radius
assign_radii!
assign_ball_radii!
load_ball_radii_table
icosphere
surface_area
read_msms
```

### Numerical SAS (Eisenhaber/Argos)
```@docs
NumericalSASResult
compute_numerical_sas
sas_area
sas_volume
```

### Reduced surface
```@docs
ReducedSurface
compute_reduced_surface
```

### Analytical solvent-accessible surface
```@docs
SolventAccessibleSurface
compute_sas
```

### Analytical solvent-excluded surface
```@docs
SolventExcludedSurface
compute_ses
ses_area
```

### Triangulation
```@docs
triangulate
triangulate_sas
triangulate_ses
```

### Mesh cleanup and validation
```@docs
clean_ses!
check_rs
check_ses
```
