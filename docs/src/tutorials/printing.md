# 3D-printing molecules

BiochemicalAlgorithms.jl can export molecular geometry to formats that any
desktop 3D-printer slicer reads (PrusaSlicer, OrcaSlicer, Bambu Studio,
Cura). Two workflows are supported:

1. **Surface printing** — export the analytical solvent-excluded surface as
   a single watertight STL or 3MF model, optionally with per-face CPK
   colouring driven by the underlying atoms.
2. **Construction-kit printing** — turn a small molecule into a set of
   colour-coded atom spheres (with friction-fit sockets) and bond cylinders
   (with matching pegs), ready to snap together.

All examples use molecules shipped with the test data, so you can paste
them into a REPL with the package activated:

```julia
using BiochemicalAlgorithms
ball_data_path("../test/data/AlaAla.pdb")   # absolute path string
```

## Surface printing

```julia
using BiochemicalAlgorithms

sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), Float32)
assign_radii!(sys)

mesh = triangulate_ses(sys; probe_radius=1.5, density=2.0)
export_stl(mesh, "bpti_ses.stl")
```

The mesh comes out of [`triangulate_ses`](@ref) in **Ångströms**. Slicers
default to **millimetres**, so by default the raw STL will appear ~10× too
small. Either scale the mesh manually before exporting, or use
[`export_ses_3mf`](@ref) which handles the unit conversion for you.

### Coloured SES with one call

[`export_ses_3mf`](@ref) triangulates the SES, tags every triangle with the
CPK colour of its nearest atom, scales to mm, and writes a 3MF in one step.
Multi-material slicers (Prusa XL, Bambu A1/X1, Prusa MMU3, …) honour the
`<m:colorgroup>` extension and map the colours to filaments; mono-material
slicers fall back to the per-object base colour silently. STL has no place
for the colour data, so coloured prints must use 3MF.

```julia
sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), Float32)
assign_radii!(sys)
export_ses_3mf(sys, "bpti.3mf"; probe_radius=1.5, density=2.0)
```

The two new keyword arguments worth knowing about:

#### `max_size_mm` — fit the build volume

If you set `max_size_mm = N`, the longest dimension of the molecule's
atom-coordinate bounding box maps to **N mm**; the mesh as a whole then
fits inside an `N × N × N` mm cube (give or take the probe inflation,
which adds ~3 mm × scale per axis).

```julia
# Bambu H2C has a 256 × 256 × 256 mm build volume; 200 mm leaves margin
export_ses_3mf(sys, "bpti_fit.3mf"; density=2.0, max_size_mm=200)
```

When `max_size_mm` is set, it overrides the explicit `scale` kwarg.

#### `by_chain` — one object per chain

For proteins, you usually want each chain printed as a separate physical
object. `by_chain = true` iterates `chains(sys)`, computes the SES of each
chain independently, applies the same Å → mm scale to all of them, and
emits one `<object>` per chain into the 3MF. The chains are grid-laid-out
on the build plate so slicers that don't auto-arrange import them already
spread out.

```julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"), Float32)
infer_topology!(sys)          # reconstruct missing atoms + build bonds
assign_radii!(sys)
export_ses_3mf(sys, "2ptc_split.3mf";
    density      = 2.0,
    max_size_mm  = 200,       # whole-molecule fit, applied to every chain
    by_chain     = true,
)
```

`infer_topology!` is the right preprocessing step for PDB inputs that lack
hydrogen records — it fills them in from the fragment database. The
1-arg method requires a `System{Float32}`; for `Float64` pass a
`FragmentDB` explicitly.

### Surface quality

The triangulation density controls how polygonal the surface looks. For an
SES that's going to be displayed in the slicer preview, density 3–4 is
visually smooth; density 6.5 is overkill outside very small molecules
because file size grows linearly with the triangle count. Approximate
triangle counts on 2ptc (4587 atoms after `infer_topology!`):

| density | triangles | 3MF size |
|---|---:|---:|
| 2.0 |  59 k | 7.0 MB |
| 4.0 | 104 k | 12.3 MB |
| 6.5 | 159 k | 18.9 MB |

3MF carries triangle geometry only, no per-vertex normals — that's a
limitation of the spec, not this exporter. Slicers compute face normals
and shade flat; the only way to make a smooth preview is to triangulate
finer.

## Construction kit

[`construction_kit`](@ref) generates one [`PrintablePart`](@ref) per atom
(a sphere with cylindrical sockets drilled toward each bonded neighbour)
and one per bond (a uniform-diameter cylinder with a conic taper to a
slightly-narrower peg at each end). All parts are settled to z = 0 and
laid out on a grid so the slicer doesn't need to rearrange them.

```julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
assign_radii!(sys)

parts = construction_kit(sys; scale = 20)   # 1 Å = 20 mm
export_3mf(parts, "AlaAla_kit.3mf")
```

### Tuning the kit

Defaults are tuned to look like the kits sold by Molymod / HGS / Conatex:

| kwarg | default | meaning |
|---|---|---|
| `scale` | `10` | Å → mm conversion factor |
| `atom_scale` | `0.22` | sphere radius as fraction of vdW radius |
| `peg_radius` | `2.0` mm | socket / peg radius |
| `peg_length` | `8.0` mm | how deep the peg seats into the socket |
| `bond_radius` | `2.8` mm | visible bond shaft radius (gives a 0.8 mm chamfer at each end) |
| `subdivisions` | `4` | icosphere subdivisions (2562 verts/sphere) |
| `segments` | `24` | cylinder facet count |
| `joint` | `:peg` | `:peg` for friction-fit or `:magnet` for 3 mm × 1 mm neodymium discs |

Choose `scale` so the carbon ball comes out 7–10 mm in diameter for a kit
that fits a desk; for a teaching demo, `scale = 20` gives ~15 mm balls.

#### Joint mechanism

The default `:peg` is a friction fit: print in PETG or PLA at 0.2 mm layer
height; the 2 mm × 8 mm peg seats firmly without glue. If parts come out
too tight, reduce `peg_radius` by 0.05 mm; if too loose, raise it.

`:magnet` switches to a 1.55 mm × 1.2 mm socket geometry sized for 3 mm Ø ×
1 mm N52 neodymium discs (glued in after printing). Stronger and more
re-usable than friction-fit but needs sourcing the magnets.

#### CPK colours

Every atom part is tagged with the CPK colour returned by
[`cpk_color`](@ref); bonds are grey. 3MF carries the colour via the
materials extension that PrusaSlicer / OrcaSlicer / Bambu Studio read; you
need a filament for each distinct colour you want printed. STL has no
colour data, so STL kits print mono-coloured.

#### Multi-bond rendering

Bond order ≥ 2 produces the "banana"/"curved-stick" visualisation seen
in commercial molecular kits:

* **Double bond** → two curved cylinders bowing outward in opposite
  perpendicular directions. Each cylinder is a separate printed
  `PrintablePart` (named `bond-N-cyl-1` and `bond-N-cyl-2`).
* **Triple bond** → three curved cylinders at 120° around the bond
  axis.
* **Single bond** stays a single straight cylinder.

Atom spheres on each end of a multi-bond receive multiple sockets — one
per parallel cylinder — at ±17° tilt from the bond axis, so each peg
seats cleanly into its own radial socket. Each curved cylinder uses a
cubic-Bezier sweep whose tangent at both endpoints matches the socket
axis, giving a smooth fit at the atom surface.

## What to check in the slicer

When you load the file:

- **Imported objects** should match what the call produced (one per kit
  part, one per chain when `by_chain=true`, one for a whole-molecule SES).
- **Bounding box** in the slicer's status bar should be at most
  `max_size_mm × max_size_mm × max_size_mm` (probe inflation can push the
  longest axis past `max_size_mm` by 1.5 Å × scale; reduce `max_size_mm`
  if you need a hard cap).
- **Watertightness**: the slicer might report a handful of non-manifold
  edges at high SES densities — they come from the surface
  triangulator's saturation at faces with very high curvature. Slicers
  accept them and print fine.
- **Multi-material colour preview**: turn on multi-material in your
  slicer to see the per-atom CPK colours; otherwise they appear monocrome.
