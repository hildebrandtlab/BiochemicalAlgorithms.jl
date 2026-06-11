# 3D-printing molecules

BiochemicalAlgorithms.jl can export molecular geometry to formats that any
desktop 3D printer slicer understands (PrusaSlicer, OrcaSlicer, Bambu Studio,
Cura). Two workflows are supported:

1. **Surface printing** — export a triangulated molecular surface
   (`triangulate_ses`, `triangulate_sas`, …) as a single watertight STL or
   3MF model.
2. **Construction-kit printing** — turn a small molecule into a set of
   colour-coded atom spheres and bond cylinders that snap together with
   friction-fit pegs (default) or magnets.

## Exporting a surface

```julia
using BiochemicalAlgorithms
sys = load_pdb("bpti.pdb")
assign_radii!(sys)
mesh = triangulate_ses(sys; probe_radius=1.5, density=2.0)

# Single-mesh STL — every slicer reads this.
export_stl(mesh, "bpti_ses.stl")

# Or a 3MF with units + a base colour (slicers display the colour preview).
export_3mf(mesh, "bpti_ses.3mf"; color=(0.7, 0.8, 1.0), name="BPTI SES")
```

The mesh comes out of the surface routines in **Ångströms**. Because 3D
printer slicers default to **millimetres**, you may want to scale before
exporting (e.g. `2.0` gives a model in cm/Å). The construction-kit
orchestrator below does this automatically.

### Coloured SES print

For multi-material printers (Prusa XL, Bambu A1/X1/H2D, Prusa MMU3, Bambu
AMS), each triangle of the SES can carry the CPK colour of its underlying
atom. The 3MF specification stores this as a `<colorgroup>` resource that
slicers read into the per-face filament map:

```julia
sys = load_pdb("bpti.pdb")
assign_radii!(sys)

# Two-call form: lets you reuse the SES for other work.
ses  = compute_ses(sys; probe_radius=1.5)
mesh = triangulate_ses(ses; density=2.0)
fc   = ses_face_colors_by_atom(mesh, ses, sys)
# Build a PrintablePart manually if you want a custom name/scale.

# One-call convenience: triangulate, colour, and write a coloured 3MF.
export_ses_3mf(sys, "bpti_ses_colored.3mf"; probe_radius=1.5, density=2.0)
```

Each triangle is mapped to its nearest atom's CPK colour ([`cpk_color`](@ref))
— a centroid-to-atom-centre lookup, which is accurate for Connolly SES
meshes because every triangle sits on or directly above a single atom.
Slicers that do not understand the `m:colorgroup` extension fall back to
the per-object base colour silently. STL has no place to store the colour
data, so for coloured prints you must use 3MF.

## Construction kit

[`construction_kit`](@ref) generates one printable [`PrintablePart`](@ref) per
atom and per bond. Each atom is a sphere with cylindrical pegs pointing
toward every bonded neighbour; each bond is a cylinder with matching sockets
at each end. Bond orders 2 and 3 produce visually parallel (fused) shafts on
one printed piece, with one peg per atom.

```julia
sys = load_sdfile("ethanol.sdf")[1]
assign_radii!(sys)

parts = construction_kit(sys;
    scale       = 10,    # Å → mm (10 mm/Å gives ~1 cm spheres)
    atom_scale  = 0.4,   # sphere radius as fraction of vdW
    peg_radius  = 1.5,   # mm
    peg_length  = 8.0,   # mm
    bond_radius = 1.5,   # mm
    joint       = :peg,  # or :magnet for 3×1 mm neodymium discs
)

# Write as a single 3MF — slicers auto-arrange all parts on the bed.
export_3mf(parts, "ethanol_kit.3mf")
# Or one STL per part:
export_stl(parts, "ethanol_kit/")
```

### Joint mechanism

The default `joint = :peg` produces friction-fit cylindrical pegs. Print
with PETG or PLA at 0.2 mm layer height — the default 3 mm Ø × 8 mm peg
seats firmly without glue. If parts come out too tight, scale the peg radius
down by 0.05 mm and reprint; if too loose, scale up.

`joint = :magnet` switches to a 3 mm Ø × 1 mm flat-magnet socket geometry,
suitable for 3 mm Ø × 1 mm N52 neodymium discs glued in after print. The
magnetic connector is stronger than friction-fit and survives more
disassembly cycles, but requires sourcing the magnets.

### CPK colours

Each atom part is tagged with a Corey–Pauling–Koltun colour
([`cpk_color`](@ref)). 3MF carries the colour as a base-material attribute;
slicers display it in the preview and (on multi-material printers) drive the
filament change. STL has no colour data, so the colour is lost on STL
export — print each atom from a separate STL file with the matching
filament instead, or use the 3MF.

### Tuning the geometry for your printer

* `scale` — Å → mm. `10` (default) gives carbon spheres of ~7 mm and
  C–C bonds of ~15 mm overall. Increase to `15` or `20` for larger models;
  decrease to `5` if you want a compact kit.
* `atom_scale` — fraction of vdW radius. `0.4` matches commercial molecular
  kits; raise toward `0.6` for chunkier spheres if pegs collide.
* `subdivisions` — icosphere refinement (`2` = 162 vertices/sphere). Raise
  to `3` for smoother spheres (642 verts) at print sizes ≥ 2 cm.
* `segments` — cylinder facet count. `24` looks smooth at typical print
  sizes; lower it to `12` for faster export of huge molecules.
