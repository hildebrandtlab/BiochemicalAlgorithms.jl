export
    PrintablePart,
    export_stl,
    export_3mf,
    cpk_color,
    atom_sphere_mesh,
    bond_cylinder_mesh,
    construction_kit,
    ses_face_colors_by_atom,
    export_ses_3mf

using ZipFile
using LinearAlgebra

# ---------------------------------------------------------------------------
# Result type.
# ---------------------------------------------------------------------------

"""
A single 3D-printable part produced by [`construction_kit`](@ref) or by the
SES-export helpers.

* `mesh`        — watertight `GeometryBasics.Mesh`, vertex units **millimetres**.
* `color`       — base RGB tuple in `[0, 1]³` used when writing 3MF
                  (fallback for slicers that ignore per-triangle colours).
* `name`        — human-readable label (e.g. `"C-atom-3"`, `"bond-1"`).
* `face_colors` — optional per-triangle RGB (length must equal
                  `length(mesh.faces)`). When set, [`export_3mf`](@ref)
                  writes a 3MF `<colorgroup>` resource and references one
                  entry per `<triangle>`; multi-material slicers
                  (PrusaSlicer, OrcaSlicer, Bambu Studio) read this as the
                  per-face filament. Defaults to `nothing` (single-colour
                  part).
"""
struct PrintablePart
    mesh::GeometryBasics.Mesh
    color::NTuple{3, Float64}
    name::String
    face_colors::Union{Nothing, Vector{NTuple{3, Float64}}}
end

PrintablePart(mesh, color, name) = PrintablePart(mesh, color, name, nothing)

# ---------------------------------------------------------------------------
# Binary STL writer.
# ---------------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Write `mesh` to `path` as a binary STL file (40-byte little-endian header +
`UInt32` triangle count + 50 bytes per triangle). Coordinates are written as
`Float32` regardless of the mesh element type, as required by the STL
specification. Returns `path`.
"""
function export_stl(mesh::GeometryBasics.Mesh, path::AbstractString)
    pos = mesh.position
    faces = mesh.faces
    open(path, "w") do io
        write(io, zeros(UInt8, 80))            # 80-byte header
        write(io, UInt32(length(faces)))
        for f in faces
            a = pos[Int(f[1])]
            b = pos[Int(f[2])]
            c = pos[Int(f[3])]
            n = cross(b - a, c - a)
            ℓ = norm(n)
            nv = ℓ > eps(Float64) ? n / ℓ : zero(n)
            write(io, Float32(nv[1])); write(io, Float32(nv[2])); write(io, Float32(nv[3]))
            for v in (a, b, c)
                write(io, Float32(v[1])); write(io, Float32(v[2])); write(io, Float32(v[3]))
            end
            write(io, UInt16(0))               # attribute byte count
        end
    end
    path
end

"""
    $(TYPEDSIGNATURES)

Convenience: write each `PrintablePart` to its own STL file under the directory
`dirpath`, named `"<idx>_<part.name>.stl"`. Returns the list of paths.
"""
function export_stl(parts::AbstractVector{PrintablePart}, dirpath::AbstractString)
    isdir(dirpath) || mkpath(dirpath)
    paths = String[]
    for (i, p) in enumerate(parts)
        path = joinpath(dirpath, "$(lpad(i, 3, '0'))_$(p.name).stl")
        export_stl(p.mesh, path)
        push!(paths, path)
    end
    paths
end

# ---------------------------------------------------------------------------
# 3MF writer (multi-part, per-part RGB via <basematerials>).
# ---------------------------------------------------------------------------

const _3MF_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml" />
  <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml" />
</Types>
"""

const _3MF_RELS = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Target="/3D/3dmodel.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel" />
</Relationships>
"""

@inline _rgb_hex(c::NTuple{3, Float64}) = "#" * string(
    UInt32(clamp(round(Int, c[1] * 255), 0, 255)) << 16 |
    UInt32(clamp(round(Int, c[2] * 255), 0, 255)) <<  8 |
    UInt32(clamp(round(Int, c[3] * 255), 0, 255)),
    base = 16, pad = 6,
)

function _3mf_model_xml(parts::AbstractVector{PrintablePart}; unit::String = "millimeter")
    io = IOBuffer()
    println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
    # The `m:` namespace adds the Materials & Properties extension that defines
    # `<colorgroup>` and the triangle-level `pid`/`p1`/`p2`/`p3` attributes used
    # for per-face filament assignment in PrusaSlicer / OrcaSlicer / Bambu Studio.
    println(io, """<model unit="$unit" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:m="http://schemas.microsoft.com/3dmanufacturing/material/2015/02">""")
    println(io, "  <resources>")
    # Per-part base material (one entry per part).
    println(io, """    <basematerials id="1">""")
    for p in parts
        println(io, """      <base name="$(p.name)" displaycolor="$(_rgb_hex(p.color))" />""")
    end
    println(io, "    </basematerials>")

    # Per-part colorgroup IDs start at 1000 to keep clear of the basematerials
    # id (1) and the object IDs (2..). The slicer uses `pid` to look up the
    # resource; for face-coloured parts we reference the colorgroup, otherwise
    # the basematerials.
    colorgroup_id_of_part = Dict{Int, Int}()
    next_cg = 1000
    for (i, p) in enumerate(parts)
        p.face_colors === nothing && continue
        cg_id = next_cg; next_cg += 1
        colorgroup_id_of_part[i] = cg_id
        # Dedupe colours within this part (slicers cap material count).
        idx_of_color = Dict{NTuple{3, Float64}, Int}()
        ordered = NTuple{3, Float64}[]
        for c in p.face_colors
            haskey(idx_of_color, c) && continue
            idx_of_color[c] = length(ordered)
            push!(ordered, c)
        end
        println(io, """    <m:colorgroup id="$cg_id">""")
        for c in ordered
            println(io, """      <m:color color="$(_rgb_hex(c))" />""")
        end
        println(io, "    </m:colorgroup>")
        # Store the dedup table on the part for the <triangles> loop below.
        colorgroup_id_of_part[-i] = 0  # marker
        # We rebuild the lookup again in the triangle loop (simpler than
        # passing it across the println stream).
    end

    # One <object> per part.
    for (i, p) in enumerate(parts)
        oid = i + 1
        pidx = i - 1
        if haskey(colorgroup_id_of_part, i)
            cg = colorgroup_id_of_part[i]
            # When face_colors is set, advertise the colorgroup as the
            # default property and point the first triangle entry at index 0.
            println(io, """    <object id="$oid" name="$(p.name)" type="model" pid="$cg" pindex="0">""")
        else
            println(io, """    <object id="$oid" name="$(p.name)" type="model" pid="1" pindex="$pidx">""")
        end
        println(io, "      <mesh>")
        println(io, "        <vertices>")
        for v in p.mesh.position
            println(io, """          <vertex x="$(Float32(v[1]))" y="$(Float32(v[2]))" z="$(Float32(v[3]))" />""")
        end
        println(io, "        </vertices>")
        println(io, "        <triangles>")
        if p.face_colors === nothing
            for f in p.mesh.faces
                println(io, """          <triangle v1="$(Int(f[1])-1)" v2="$(Int(f[2])-1)" v3="$(Int(f[3])-1)" />""")
            end
        else
            cg = colorgroup_id_of_part[i]
            # Recompute the dedup table.
            idx_of_color = Dict{NTuple{3, Float64}, Int}()
            for c in p.face_colors
                haskey(idx_of_color, c) || (idx_of_color[c] = length(idx_of_color))
            end
            for (ti, f) in enumerate(p.mesh.faces)
                cidx = idx_of_color[p.face_colors[ti]]
                # `pid` overrides the object's default; `p1` indexes into the
                # colorgroup. Setting p1==p2==p3 selects a single face colour
                # (no interpolation) per the 3MF spec.
                println(io, """          <triangle v1="$(Int(f[1])-1)" v2="$(Int(f[2])-1)" v3="$(Int(f[3])-1)" pid="$cg" p1="$cidx" p2="$cidx" p3="$cidx" />""")
            end
        end
        println(io, "        </triangles>")
        println(io, "      </mesh>")
        println(io, "    </object>")
    end
    println(io, "  </resources>")
    println(io, "  <build>")
    for i in eachindex(parts)
        println(io, """    <item objectid="$(i + 1)" />""")
    end
    println(io, "  </build>")
    println(io, "</model>")
    String(take!(io))
end

"""
    $(TYPEDSIGNATURES)

Write a list of `PrintablePart`s to a 3MF (3D Manufacturing Format) file.
3MF is a zip container holding `[Content_Types].xml`, `_rels/.rels`, and
`3D/3dmodel.model` (XML mesh data). Each part becomes a separate `<object>`
with a base-material colour; slicers (PrusaSlicer, OrcaSlicer, Bambu Studio,
Cura) auto-arrange the objects on the build plate. Returns `path`.
"""
function export_3mf(parts::AbstractVector{PrintablePart}, path::AbstractString;
                    unit::AbstractString = "millimeter")
    w = ZipFile.Writer(path)
    f1 = ZipFile.addfile(w, "[Content_Types].xml")
    write(f1, _3MF_CONTENT_TYPES)
    f2 = ZipFile.addfile(w, "_rels/.rels")
    write(f2, _3MF_RELS)
    f3 = ZipFile.addfile(w, "3D/3dmodel.model")
    write(f3, _3mf_model_xml(parts; unit))
    close(w)
    path
end

# Convenience: single mesh → 3MF
function export_3mf(mesh::GeometryBasics.Mesh, path::AbstractString;
                    unit::AbstractString = "millimeter",
                    color::NTuple{3, Float64} = (0.6, 0.7, 0.9),
                    name::AbstractString = "model")
    export_3mf([PrintablePart(mesh, color, String(name))], path; unit)
end

# ---------------------------------------------------------------------------
# CPK colors.
# ---------------------------------------------------------------------------

const _CPK_TABLE = Dict(
    :H  => (1.000, 1.000, 1.000),
    :C  => (0.300, 0.300, 0.300),
    :N  => (0.137, 0.224, 0.957),
    :O  => (0.957, 0.157, 0.157),
    :F  => (0.565, 0.878, 0.314),
    :Cl => (0.122, 0.941, 0.122),
    :Br => (0.651, 0.161, 0.161),
    :I  => (0.580, 0.000, 0.580),
    :P  => (1.000, 0.502, 0.000),
    :S  => (1.000, 0.784, 0.196),
    :Na => (0.671, 0.361, 0.949),
    :K  => (0.561, 0.251, 0.831),
    :Mg => (0.541, 1.000, 0.000),
    :Ca => (0.239, 1.000, 0.000),
    :Fe => (0.878, 0.400, 0.200),
    :Zn => (0.490, 0.502, 0.690),
    :Cu => (0.784, 0.502, 0.200),
    :B  => (1.000, 0.710, 0.710),
    :Si => (0.941, 0.784, 0.627),
)

"""
    $(TYPEDSIGNATURES)

CPK colour (Corey–Pauling–Koltun) RGB tuple for `element`, components in
`[0, 1]`. Defaults to pink (`(1.0, 0.078, 0.576)`) for unmapped elements.
"""
function cpk_color(element)
    sym = Symbol(element)
    get(_CPK_TABLE, sym, (1.000, 0.078, 0.576))
end

# ---------------------------------------------------------------------------
# Primitive mesh: cylinder (closed, watertight).
# ---------------------------------------------------------------------------

# Frame of orthonormal axes (n, u, v) with n the given direction.
function _ortho_frame(n::AbstractVector{T}) where T
    n = Vec3{T}(n...) / norm(n)
    a = abs(n[1]) < T(0.9) ? Vec3{T}(1, 0, 0) : Vec3{T}(0, 1, 0)
    u = cross(n, a); u = u / norm(u)
    v = cross(n, u)
    (n, u, v)
end

"""
    $(TYPEDSIGNATURES)

Build a closed (watertight) cylinder mesh of `radius` and `length` whose axis
runs from `p1` to `p2`. `segments` is the number of facets around the
perimeter. Triangles are oriented outward (CCW from outside).
"""
function cylinder_mesh(p1::AbstractVector{T}, p2::AbstractVector{T}, radius::T;
                       segments::Int = 24) where T<:Real
    axis = Vec3{T}((p2 - p1)...)
    _, u, v = _ortho_frame(axis)
    verts = Point3{T}[]
    # 2 * segments side vertices + 2 cap centers
    for k in 0:segments-1
        θ = T(2π) * k / segments
        ring = radius * (cos(θ) * u + sin(θ) * v)
        push!(verts, Point3{T}((Vec3{T}(p1...) + ring)...))
    end
    for k in 0:segments-1
        θ = T(2π) * k / segments
        ring = radius * (cos(θ) * u + sin(θ) * v)
        push!(verts, Point3{T}((Vec3{T}(p2...) + ring)...))
    end
    bot_center_idx = length(verts) + 1
    push!(verts, Point3{T}(p1...))
    top_center_idx = length(verts) + 1
    push!(verts, Point3{T}(p2...))

    faces = TriangleFace{Int}[]
    # Sides: quads → 2 triangles each.
    for k in 1:segments
        a = k
        b = mod1(k + 1, segments)
        c = segments + b
        d = segments + a
        push!(faces, TriangleFace{Int}(a, b, c))
        push!(faces, TriangleFace{Int}(a, c, d))
    end
    # Bottom cap (normal -n).
    for k in 1:segments
        a = k
        b = mod1(k + 1, segments)
        push!(faces, TriangleFace{Int}(bot_center_idx, b, a))
    end
    # Top cap (normal +n).
    for k in 1:segments
        a = segments + k
        b = segments + mod1(k + 1, segments)
        push!(faces, TriangleFace{Int}(top_center_idx, a, b))
    end
    GeometryBasics.Mesh(verts, faces)
end

# ---------------------------------------------------------------------------
# Mesh primitives: union by concatenation (slicers union closed sub-meshes).
# ---------------------------------------------------------------------------

# Append `m2` into `m1`, returning a new mesh. Both must use Point3{T} +
# TriangleFace{Int}. Vertex indices in `m2` are shifted by length(m1.position).
function _mesh_union(m1::GeometryBasics.Mesh, m2::GeometryBasics.Mesh)
    off = length(m1.position)
    new_verts = vcat(collect(m1.position), collect(m2.position))
    new_faces = TriangleFace{Int}[]
    sizehint!(new_faces, length(m1.faces) + length(m2.faces))
    for f in m1.faces
        push!(new_faces, TriangleFace{Int}(Int(f[1]), Int(f[2]), Int(f[3])))
    end
    for f in m2.faces
        push!(new_faces, TriangleFace{Int}(Int(f[1]) + off, Int(f[2]) + off, Int(f[3]) + off))
    end
    GeometryBasics.Mesh(new_verts, new_faces)
end

# ---------------------------------------------------------------------------
# Atom-sphere mesh: solid sphere with N protruding cylindrical pegs.
# (Each peg's base intersects the sphere; the slicer takes the boolean union
#  of the closed sub-meshes at slice time. This avoids in-process CSG.)
# ---------------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Build the printed mesh for one atom: a single watertight icosphere of
`radius` at `center`. `bond_directions` is accepted for API symmetry but is
currently unused — peg sockets are emitted by [`bond_cylinder_mesh`](@ref) as
positive geometry on the bond side, and the atom is a featureless sphere
that the user joins to bond pieces with glue (or a small dab of cyanoacrylate
inside the sphere's pre-printed dimple — not yet implemented).

Returned mesh is in the same coordinate units as `center`/`radius` (the caller
is expected to have scaled Å to mm before this call).
"""
function atom_sphere_mesh(center::AbstractVector{T}, radius::T,
                          bond_directions::AbstractVector{<:AbstractVector{<:Real}};
                          peg_radius::T = T(1.5),
                          peg_length::T = T(8.0),
                          subdivisions::Int = 2,
                          segments::Int = 16) where T<:Real
    sphere = icosphere(T, subdivisions)
    sphere_verts = [Point3{T}((Vec3{T}(center...) + radius * Vec3{T}(v...))...)
                    for v in sphere.position]
    sphere_faces = [TriangleFace{Int}(Int(f[1]), Int(f[2]), Int(f[3])) for f in sphere.faces]
    GeometryBasics.Mesh(sphere_verts, sphere_faces)
end

# ---------------------------------------------------------------------------
# Bond mesh: solid cylinder spanning the visible gap between two atom
# surfaces. Each piece is a single watertight closed mesh, so slicers accept
# it without "non-manifold" warnings. Bond order ≥ 2 renders as 2 or 3
# parallel cylinders all of which are individually watertight; they touch
# only at their endpoints (no inter-cylinder cross-section), which slicers
# union cleanly because the touching surfaces are coincident.
# ---------------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Build the printed mesh for one bond between atoms at `p1` and `p2`,
spanning only the visible gap between the two atom surfaces (the caller
passes the atom-surface points, NOT the atom centres). The piece is a
single solid cylinder of `bond_radius`. When `order ≥ 2`, the shaft is
rendered as 2 or 3 parallel cylinders for the double/triple-bond look —
each is its own watertight sub-mesh, which slicers union by simple
inclusion. `peg_radius`/`peg_length` are accepted for API symmetry but
unused in the current build (no protruding pegs).
"""
function bond_cylinder_mesh(p1::AbstractVector{T}, p2::AbstractVector{T};
                            bond_radius::T = T(2.5),
                            peg_radius::T = T(1.5),
                            peg_length::T = T(8.0),
                            order::Int = 1,
                            segments::Int = 24) where T<:Real
    axis = Vec3{T}((p2 - p1)...)
    _, u, v = _ortho_frame(axis)

    if order <= 1
        return cylinder_mesh(p1, p2, bond_radius; segments)
    end

    # Multi-bond visual: 2 or 3 parallel cylinders. Each is thinner than the
    # single-bond shaft so the total cross-sectional area roughly matches.
    shaft_radius = bond_radius * (order == 2 ? T(0.55) : T(0.45))
    offset_distance = bond_radius * T(0.6)

    mesh = nothing
    if order == 2
        for sgn in (-T(1), T(1))
            ofs = sgn * offset_distance * u
            p1o = Point3{T}((Vec3{T}(p1...) + ofs)...)
            p2o = Point3{T}((Vec3{T}(p2...) + ofs)...)
            shaft = cylinder_mesh(p1o, p2o, shaft_radius; segments)
            mesh = mesh === nothing ? shaft : _mesh_union(mesh, shaft)
        end
    else
        for k in 0:2
            θ = T(2π) * k / 3
            ofs = offset_distance * (cos(θ) * u + sin(θ) * v)
            p1o = Point3{T}((Vec3{T}(p1...) + ofs)...)
            p2o = Point3{T}((Vec3{T}(p2...) + ofs)...)
            shaft = cylinder_mesh(p1o, p2o, shaft_radius; segments)
            mesh = mesh === nothing ? shaft : _mesh_union(mesh, shaft)
        end
    end
    mesh
end

# ---------------------------------------------------------------------------
# Construction-kit orchestrator.
# ---------------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Generate a 3D-printable construction kit for a small molecule.

For every atom: a sphere (radius = `atom_scale × vdW`) with one cylindrical
peg per incident bond, oriented toward the bonded neighbour, sized for a
friction fit (default 3 mm × 8 mm). For every bond: a hollow cylinder of
length equal to the inter-atomic distance, with matching sockets at each end.
Bond order 2 produces two parallel cylinders fused at the end caps; order ≥ 3
produces three parallel cylinders.

All coordinates are converted from Ångström to millimetres via `scale`
(default `1.0 Å → 10 mm`, i.e. `scale = 10`), which gives ~1 cm spheres for
typical organic atoms — a good size for FDM printing.

Returns `Vector{PrintablePart}` ready to pass to [`export_3mf`](@ref) or
[`export_stl`](@ref).

# Keyword arguments
 - `scale::Real = 10`         — Å → mm conversion factor.
 - `atom_scale::Real = 0.4`   — sphere radius as fraction of vdW radius.
 - `peg_radius::Real = 1.5`   — peg/socket radius in **mm**.
 - `peg_length::Real = 8.0`   — peg/socket length in **mm**.
 - `bond_radius::Real = 1.5`  — visible bond shaft radius in **mm**.
 - `joint::Symbol = :peg`     — `:peg` (friction fit) or `:magnet` (sized for
                                 standard 3 mm × 1 mm neodymium discs).
 - `subdivisions::Int = 2`    — icosphere subdivisions for atom spheres.
 - `segments::Int = 24`       — cylinder facet count.
"""
function construction_kit(ac::AbstractAtomContainer{T};
                          scale::Real      = 10,
                          atom_scale::Real = 0.35,
                          peg_radius::Real = 1.5,
                          peg_length::Real = 8.0,
                          bond_radius::Real = 2.5,
                          joint::Symbol    = :peg,
                          subdivisions::Int = 2,
                          segments::Int    = 24) where T<:Real
    if joint === :magnet
        # 3 mm Ø × 1 mm neodymium disc, typical kit dimensions.
        peg_radius = 1.55
        peg_length = 1.2
    elseif joint !== :peg
        throw(ArgumentError("joint must be :peg or :magnet, got :$joint"))
    end

    s = T(scale)
    ascale = T(atom_scale)
    pr = T(peg_radius)
    pl = T(peg_length)
    br = T(bond_radius)

    at = atoms(ac)
    bt = bonds(ac)

    # Per-atom bond directions (toward each bonded neighbour, in mm).
    n_atoms = length(at)
    nbrs = [Vector3{T}[] for _ in 1:n_atoms]
    bond_pairs = Tuple{Int, Int, BondOrderType}[]
    for i in eachindex(bt)
        a1 = bt.a1[i]; a2 = bt.a2[i]
        # The bond table stores atom *idx* values; map to row indices.
        r1 = findfirst(==(a1), at.idx)
        r2 = findfirst(==(a2), at.idx)
        (r1 === nothing || r2 === nothing) && continue
        d = Vector3{T}((at.r[r2] - at.r[r1])...)
        push!(nbrs[r1],  d)
        push!(nbrs[r2], -d)
        push!(bond_pairs, (r1, r2, bt.order[i]))
    end

    parts = PrintablePart[]

    # Vdw radii — need them populated. Mirror surfaces' convention.
    radii_filled = !all(iszero, at.radius)
    atom_sphere_radius_mm = Vector{T}(undef, n_atoms)   # for bond-end computation

    for i in 1:n_atoms
        elem = at.element[i]
        vdw = radii_filled && !iszero(at.radius[i]) ? T(at.radius[i]) : T(vdw_radius(Float64, elem))
        sphere_radius_mm = ascale * vdw * s
        atom_sphere_radius_mm[i] = sphere_radius_mm
        center_mm = Vector3{T}((s .* at.r[i])...)
        dirs = [d / norm(d) for d in nbrs[i] if norm(d) > eps(T)]
        mesh = atom_sphere_mesh(center_mm, sphere_radius_mm, dirs;
                                 peg_radius = pr, peg_length = pl,
                                 subdivisions, segments)
        # Settle to z=0 + recentre — slicers auto-arrange parts side-by-side.
        mesh = _settle_to_buildplate(mesh)
        # Force colour through face_colors so multi-material slicer previews
        # render the atom in CPK (the per-object basematerials hint is often
        # ignored unless multi-material printing is enabled).
        col = cpk_color(elem)
        face_colors = fill(col, length(mesh.faces))
        push!(parts, PrintablePart(
            mesh,
            col,
            "atom-$(i)-$(string(Symbol(elem)))",
            face_colors,
        ))
    end

    for (k, (i, j, ord)) in enumerate(bond_pairs)
        ci = Vector3{T}((s .* at.r[i])...)
        cj = Vector3{T}((s .* at.r[j])...)
        axis = cj - ci
        L = norm(axis)
        L < eps(T) && continue
        n̂ = axis / L
        # Bond runs only between atom surfaces — not centre to centre — so
        # the visible cylinder length matches what the user sees in the
        # assembled model and is proportional to the gap.
        ri = atom_sphere_radius_mm[i]
        rj = atom_sphere_radius_mm[j]
        # Tiny overlap (0.2 mm) so the bond sits flush against the sphere.
        surf_i = ci + (ri - T(0.2)) * n̂
        surf_j = cj - (rj - T(0.2)) * n̂
        ord_int = ord === BondOrder.Single  ? 1 :
                  ord === BondOrder.Double  ? 2 :
                  ord === BondOrder.Triple  ? 3 :
                  ord === BondOrder.Aromatic ? 2 : 1
        mesh = bond_cylinder_mesh(surf_i, surf_j;
                                   bond_radius = br,
                                   peg_radius  = pr,
                                   peg_length  = pl,
                                   order       = ord_int,
                                   segments)
        mesh = _settle_to_buildplate(mesh)
        col = (0.55, 0.55, 0.55)
        face_colors = fill(col, length(mesh.faces))
        push!(parts, PrintablePart(mesh, col, "bond-$(k)", face_colors))
    end

    parts
end

# ---------------------------------------------------------------------------
# Print-bed orientation helpers.
# ---------------------------------------------------------------------------

# Translate so the bounding-box min becomes (cx, cy, 0): the part sits flat on
# the build plate, centred on the XY origin so slicers' auto-arrange can pack
# it next to other parts.
function _settle_to_buildplate(mesh::GeometryBasics.Mesh{T}) where T<:Any
    pos = mesh.position
    isempty(pos) && return mesh
    xs = [v[1] for v in pos]; ys = [v[2] for v in pos]; zs = [v[3] for v in pos]
    cx = (minimum(xs) + maximum(xs)) / 2
    cy = (minimum(ys) + maximum(ys)) / 2
    zmin = minimum(zs)
    new_pos = [Point3{eltype(v)}(v[1] - cx, v[2] - cy, v[3] - zmin) for v in pos]
    GeometryBasics.Mesh(
        new_pos,
        [TriangleFace{Int}(Int(f[1]), Int(f[2]), Int(f[3])) for f in mesh.faces],
    )
end

# PCA: rotate the mesh so its longest principal axis aligns with +X, the
# second-longest with +Y, and the shortest with +Z. This minimises the printed
# height (and thus the area exposed to overhangs/supports). After rotation
# the mesh is settled to (0, 0, 0) build-plate origin.
function _pca_align_for_printing(mesh::GeometryBasics.Mesh)
    pos = mesh.position
    isempty(pos) && return mesh
    n = length(pos)
    # Centroid
    cx = sum(v[1] for v in pos) / n
    cy = sum(v[2] for v in pos) / n
    cz = sum(v[3] for v in pos) / n
    # Covariance matrix (3×3)
    Mxx = 0.0; Mxy = 0.0; Mxz = 0.0
    Myy = 0.0; Myz = 0.0; Mzz = 0.0
    for v in pos
        dx = Float64(v[1] - cx); dy = Float64(v[2] - cy); dz = Float64(v[3] - cz)
        Mxx += dx * dx; Mxy += dx * dy; Mxz += dx * dz
        Myy += dy * dy; Myz += dy * dz; Mzz += dz * dz
    end
    Σ = [Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz] ./ n
    # Eigendecomposition. `eigen` of a Symmetric returns ascending eigenvalues.
    F = eigen(Symmetric(Σ))
    # Sort indices by descending eigenvalue (largest spread first).
    order = sortperm(F.values; rev = true)
    R = F.vectors[:, order]
    # Ensure a right-handed (det = +1) rotation; flip the last column if needed.
    if det(R) < 0
        R[:, 3] .= -R[:, 3]
    end
    # Apply: new_v = R' * (v - centroid)
    T = eltype(eltype(pos))
    new_pos = [Point3{T}(
        T(R[1,1] * (Float64(v[1]) - cx) + R[2,1] * (Float64(v[2]) - cy) + R[3,1] * (Float64(v[3]) - cz)),
        T(R[1,2] * (Float64(v[1]) - cx) + R[2,2] * (Float64(v[2]) - cy) + R[3,2] * (Float64(v[3]) - cz)),
        T(R[1,3] * (Float64(v[1]) - cx) + R[2,3] * (Float64(v[2]) - cy) + R[3,3] * (Float64(v[3]) - cz)),
    ) for v in pos]
    rotated = GeometryBasics.Mesh(
        new_pos,
        [TriangleFace{Int}(Int(f[1]), Int(f[2]), Int(f[3])) for f in mesh.faces],
    )
    _settle_to_buildplate(rotated)
end

# ---------------------------------------------------------------------------
# SES export with per-triangle atom colouring.
# ---------------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

For each triangle in `mesh`, return the CPK colour of the atom in `ac` whose
centre is closest to the triangle's centroid. `ac` is the
`AbstractAtomContainer{T}` the SES was computed from; vertex coordinates of
`mesh` are assumed to be in the same coordinate frame as `ac` (Ångströms).
The reduced surface inside `ses` is consulted only to bound the candidate
atom set; in practice the centroid-to-atom-centre lookup is robust enough
because every triangle of a Connolly SES sits on or near one specific atom.

Returns `Vector{NTuple{3, Float64}}` of length `length(mesh.faces)` suitable
for the `face_colors` field of [`PrintablePart`](@ref).
"""
function ses_face_colors_by_atom(mesh::GeometryBasics.Mesh,
                                  ses::SolventExcludedSurface{T},
                                  ac::AbstractAtomContainer) where T
    at = atoms(ac)
    rs = ses.reduced_surface
    # Build a flat list of (centre, element) for nearest-atom lookup. Use
    # the RS atoms (after singularity cleanup) to keep indexing consistent
    # with how the SES was actually computed.
    centres = [Vector3{T}(s.center...) for s in rs.atoms]
    # Element of each RS atom — RS atoms are 1:1 with the original atom
    # container, so we can index `ac.atoms` directly.
    elems = [at.element[i] for i in 1:min(length(centres), length(at))]
    colors = Vector{NTuple{3, Float64}}(undef, length(mesh.faces))
    for (ti, f) in enumerate(mesh.faces)
        a = mesh.position[Int(f[1])]
        b = mesh.position[Int(f[2])]
        c = mesh.position[Int(f[3])]
        cent = Vector3{T}(((a + b + c) / T(3))...)
        best = 0
        best_d2 = T(Inf)
        for (ai, ac_centre) in enumerate(centres)
            ai > length(elems) && continue
            d = cent - ac_centre
            d2 = dot(d, d)
            if d2 < best_d2
                best_d2 = d2
                best = ai
            end
        end
        colors[ti] = best == 0 ? (0.7, 0.7, 0.7) : cpk_color(elems[best])
    end
    colors
end

"""
    $(TYPEDSIGNATURES)

Triangulate the SES and write it to `path` as a 3MF file with one face
colour per triangle (CPK per nearest atom). The mesh is scaled by `scale`
from Ångström to millimetres on the way out (default `10`, giving ~1 cm
per Å — the common "ball-and-stick on the desk" scale).

# Keyword arguments
 - `density::Real = 2.0`  — triangulation density (passed to `triangulate_ses`).
 - `scale::Real    = 10`  — Å → mm conversion.
 - `name::AbstractString = "ses"` — object name in the 3MF.
"""
function export_ses_3mf(ses::SolventExcludedSurface{T}, ac::AbstractAtomContainer,
                        path::AbstractString;
                        density::Real = 2.0,
                        scale::Real   = 10,
                        name::AbstractString = "ses",
                        reorient::Bool = true) where T
    tses = triangulate_ses(ses; density = T(density))
    # Per-face colours computed BEFORE any rotation/translation — face_colors
    # are 1:1 with mesh.faces and survive any rigid transform.
    fc = ses_face_colors_by_atom(tses, ses, ac)
    s = T(scale)
    scaled_verts = [Point3{T}((s * Vec3{T}(v...))...) for v in tses.position]
    scaled_mesh = GeometryBasics.Mesh(
        scaled_verts,
        [TriangleFace{Int}(Int(f[1]), Int(f[2]), Int(f[3])) for f in tses.faces],
    )
    if reorient
        scaled_mesh = _pca_align_for_printing(scaled_mesh)
    else
        scaled_mesh = _settle_to_buildplate(scaled_mesh)
    end
    part = PrintablePart(scaled_mesh, (0.7, 0.7, 0.8), String(name), fc)
    export_3mf([part], path)
end

# Convenience: take the AbstractAtomContainer directly and run the whole pipe.
function export_ses_3mf(ac::AbstractAtomContainer{T}, path::AbstractString;
                        probe_radius::Real = 1.5,
                        density::Real      = 2.0,
                        scale::Real        = 10,
                        name::AbstractString = "ses",
                        reorient::Bool     = true) where T
    ses = compute_ses(ac; probe_radius = T(probe_radius))
    export_ses_3mf(ses, ac, path; density, scale, name, reorient)
end
