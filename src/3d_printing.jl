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

Build the printed mesh for one atom: an icosphere of `radius` at `center`
with one cylindrical socket of `peg_radius × peg_length` drilled radially
outward along every entry in `bond_directions`. The result is a **single
watertight mesh** — no overlapping closed sub-meshes — so slicers accept it
without "non-manifold" warnings.

The cut is procedural (no CSG library): for each socket direction `axis_k`,
let `θ_k = asin(peg_radius / radius)` be the half-angle of the spherical
cap that fits the peg. Vertices with `dot(v_unit, axis_k) > cos(θ_k)` are
classified as "in-cap"; in-cap triangles are dropped, and triangles
straddling the cap boundary are clipped against the cone by bisecting their
edges. The resulting hole rim is connected to an inner ring (translated
along `-axis_k` by `peg_length`) by a quad strip, then capped at the
bottom with a fan. Each socket adds one outer rim + one inner ring + one
bottom centre vertex.

Pass an empty `bond_directions` to get a plain solid sphere.

Use enough `subdivisions` that several icosphere triangles fall inside each
cap — for the default `peg_radius=1.5 mm` and `radius ≈ 6–11 mm` the cap
half-angle is ~8–14° and `subdivisions = 3` (642 verts) gives 12–20-sided
rims, which prints smoothly.

Returned mesh is in the same coordinate units as `center`/`radius` (the
caller is expected to have scaled Å to mm before this call).
"""
function atom_sphere_mesh(center::AbstractVector{T}, radius::T,
                          bond_directions::AbstractVector{<:AbstractVector{<:Real}};
                          peg_radius::T = T(1.5),
                          peg_length::T = T(8.0),
                          subdivisions::Int = 3,
                          segments::Int = 16) where T<:Real
    # Empty / oversized-peg fall-throughs: just a plain sphere.
    if isempty(bond_directions) || peg_radius >= radius
        sphere = icosphere(T, subdivisions)
        verts = [Point3{T}((Vec3{T}(center...) + radius * Vec3{T}(v...))...) for v in sphere.position]
        faces = [TriangleFace{Int}(Int(f[1]), Int(f[2]), Int(f[3])) for f in sphere.faces]
        return GeometryBasics.Mesh(verts, faces)
    end

    # Normalise socket directions.
    K = length(bond_directions)
    axes_n = Vector{Vec3{T}}(undef, K)
    for k in 1:K
        a = Vec3{T}(bond_directions[k][1], bond_directions[k][2], bond_directions[k][3])
        axes_n[k] = a / norm(a)
    end
    # sin(θ_k) = peg_radius / radius  ⇒  cos(θ_k) for cap classification.
    cos_θ = T(sqrt(1 - (Float64(peg_radius) / Float64(radius))^2))

    # Step 1: icosphere on the unit sphere.
    ico = icosphere(T, subdivisions)
    pos = [Vec3{T}(p...) for p in ico.position]  # will grow
    n_orig = length(pos)

    # Step 2: classify each original vertex by socket (0 = outside all).
    in_socket = zeros(Int, n_orig)
    for i in 1:n_orig
        v = pos[i]
        for k in 1:K
            if dot(v, axes_n[k]) > cos_θ
                in_socket[i] = k
                break
            end
        end
    end

    # Bisect edge (a,b) against cone-k. Returns new vertex on the unit sphere.
    # (One endpoint must be inside socket k, the other outside it.)
    function bisect_boundary(a::Int, b::Int, k::Int)
        va = pos[a]; vb = pos[b]
        ax = axes_n[k]
        lo, hi = 0.0, 1.0
        # f(t) = dot(normalize((1-t)*va + t*vb), ax) - cos_θ
        f_lo = Float64(dot(va, ax)) - Float64(cos_θ)
        for _ in 1:40
            mid = (lo + hi) / 2
            vm = (1 - mid) * va + mid * vb
            vm = vm / norm(vm)
            f_mid = Float64(dot(vm, ax)) - Float64(cos_θ)
            if (f_lo > 0) == (f_mid > 0)
                lo = mid; f_lo = f_mid
            else
                hi = mid
            end
        end
        t = (lo + hi) / 2
        vm = (1 - T(t)) * va + T(t) * vb
        vm = vm / norm(vm)
        push!(pos, vm)
        length(pos)
    end

    boundary_vert_of = Dict{Tuple{Int, Int, Int}, Int}()
    boundary_verts_by_socket = [Int[] for _ in 1:K]

    function get_boundary_vert(a::Int, b::Int, k::Int)
        key = a < b ? (a, b, k) : (b, a, k)
        if haskey(boundary_vert_of, key)
            return boundary_vert_of[key]
        end
        idx = bisect_boundary(a, b, k)
        boundary_vert_of[key] = idx
        push!(boundary_verts_by_socket[k], idx)
        idx
    end

    # Step 3: process each triangle. We assume each triangle straddles AT MOST
    # one socket (true when sockets don't overlap on the sphere — i.e. bonds
    # don't share a direction, which is always the case for real molecules).
    new_faces = TriangleFace{Int}[]
    for f in ico.faces
        a, b, c = Int(f[1]), Int(f[2]), Int(f[3])
        sa, sb, sc = in_socket[a], in_socket[b], in_socket[c]

        if sa == 0 && sb == 0 && sc == 0
            push!(new_faces, TriangleFace{Int}(a, b, c))
            continue
        end
        if sa == sb == sc != 0
            continue  # whole triangle in a socket → drop
        end

        # Pick the socket this triangle is partly inside.
        k = max(sa, sb, sc)
        ia = sa == k; ib = sb == k; ic = sc == k
        cnt_in = (ia ? 1 : 0) + (ib ? 1 : 0) + (ic ? 1 : 0)

        # Walk the triangle in its original cyclic order (a → b → c) and find
        # the two boundary crossings. Producing the OUTSIDE part with the same
        # cyclic order preserves outward winding.
        if cnt_in == 1
            # One vertex in, two out. Outside region is a quad.
            verts = (a, b, c)
            ins = (ia, ib, ic)
            # Find the in-vertex index in (1,2,3)
            iv = ins[1] ? 1 : (ins[2] ? 2 : 3)
            o1 = mod1(iv + 1, 3); o2 = mod1(iv + 2, 3)
            vi = verts[iv]; vo1 = verts[o1]; vo2 = verts[o2]
            nb_in_to_o1 = get_boundary_vert(vi, vo1, k)
            nb_o2_to_in = get_boundary_vert(vo2, vi, k)
            # Original winding: ... → vi → vo1 → vo2 → vi → ...
            # Outside quad in same cyclic order: nb_in_to_o1 → vo1 → vo2 → nb_o2_to_in
            push!(new_faces, TriangleFace{Int}(nb_in_to_o1, vo1, vo2))
            push!(new_faces, TriangleFace{Int}(nb_in_to_o1, vo2, nb_o2_to_in))
        else  # cnt_in == 2
            verts = (a, b, c)
            ins = (ia, ib, ic)
            ov = ins[1] ? (ins[2] ? 3 : 2) : 1   # the one out-vertex index
            i1 = mod1(ov + 1, 3); i2 = mod1(ov + 2, 3)
            vo = verts[ov]; vi1 = verts[i1]; vi2 = verts[i2]
            nb_o_to_i1 = get_boundary_vert(vo, vi1, k)
            nb_i2_to_o = get_boundary_vert(vi2, vo, k)
            # Outside is a triangle in the same cyclic order:
            # vo → nb_o_to_i1 → nb_i2_to_o (preserves outward normal).
            push!(new_faces, TriangleFace{Int}(vo, nb_o_to_i1, nb_i2_to_o))
        end
    end

    # Step 4: for each socket, build cylindrical wall + bottom cap.
    for k in 1:K
        bverts = unique!(sort(boundary_verts_by_socket[k]))
        length(bverts) < 3 && continue   # cap too small / no usable rim

        ax = axes_n[k]
        _, u, v = _ortho_frame(ax)
        # Sort by azimuth in the (u, v) plane perpendicular to the axis.
        function azi(idx)
            p = pos[idx]
            atan(Float64(dot(p, v)), Float64(dot(p, u)))
        end
        sorted = sort(bverts; by = azi)
        n = length(sorted)

        # Inner ring lies at axial position `cos(θ) - peg_length/radius` (in
        # unit-sphere coords). Radial position is `peg_radius/radius` along
        # the same azimuth as the corresponding outer vertex.
        axial_inner = T(Float64(cos_θ) - Float64(peg_length) / Float64(radius))
        radial = T(Float64(peg_radius) / Float64(radius))
        inner_ring = Int[]
        for bv in sorted
            p = pos[bv]
            perp = p - dot(p, ax) * ax
            perp_n = norm(perp)
            unit_perp = perp_n > eps(T) ? perp / perp_n : u
            push!(pos, radial * unit_perp + axial_inner * ax)
            push!(inner_ring, length(pos))
        end

        # Wall: quad strip, normals point TOWARD the axis (into the cylinder
        # interior — "outside" of the solid is the cylinder bore).
        # Order verified: (out_i, inner_j, out_j) gives the correct outward
        # winding away from the bore (i.e. into the bore, since the bore is
        # part of the "negative space"). Adjust if normals come out wrong.
        for i in 1:n
            j = mod1(i + 1, n)
            oi = sorted[i]; oj = sorted[j]
            ii = inner_ring[i]; ij = inner_ring[j]
            push!(new_faces, TriangleFace{Int}(oi, ij, ii))
            push!(new_faces, TriangleFace{Int}(oi, oj, ij))
        end

        # Bottom cap centre vertex. The cap is at axial position
        # `axial_inner * ax` (deeper than the sphere surface along the bond
        # axis); its outward normal must point INTO the bore — i.e. toward
        # the sphere surface, along +ax. Visiting (in_i → in_j → centre) is
        # CCW when viewed from +ax (since inner_ring is sorted CCW around
        # +ax), so the resulting normal points along +ax, which is what we
        # want. The previous winding (in_i, centre, in_j) traces CW and
        # gave a flipped normal — the cap rendered as inside-out and
        # slicers (Bambu Studio in particular) saw through to the opposite
        # side of the sphere.
        push!(pos, axial_inner * ax)
        center_idx = length(pos)
        for i in 1:n
            j = mod1(i + 1, n)
            push!(new_faces, TriangleFace{Int}(inner_ring[i], inner_ring[j], center_idx))
        end
    end

    # Step 5: scale to radius and translate to centre.
    cv = Vec3{T}(center...)
    final_verts = [Point3{T}((cv + radius * p)...) for p in pos]
    mesh = GeometryBasics.Mesh(final_verts, new_faces)
    # Step 6: compact away vertices that aren't referenced by any face. The
    # 1-in-2-out clip path replaces the in-cap icosphere vertex with two new
    # boundary verts; the original is then orphaned in `pos` even though it
    # never appears in any triangle. Slicers tolerate unreferenced verts but
    # the resulting Euler characteristic is wrong — compact for cleanliness.
    _compact_vertices(mesh)
end

# Lay out a list of (already-settled) parts in a roughly square grid on the
# build plate so the slicer doesn't have them stacked at (0, 0, 0). Each part
# is shifted in XY only; z stays untouched. `padding` is added between cells.
function _layout_grid(parts::AbstractVector{PrintablePart}; padding::Real = 5.0)
    n = length(parts)
    n == 0 && return parts
    # Footprint of each part: (xspan, yspan).
    spans = Vector{Tuple{Float64, Float64}}(undef, n)
    for (i, p) in enumerate(parts)
        xs = [Float64(v[1]) for v in p.mesh.position]
        ys = [Float64(v[2]) for v in p.mesh.position]
        spans[i] = (maximum(xs) - minimum(xs), maximum(ys) - minimum(ys))
    end
    cell_w = maximum(s[1] for s in spans) + padding
    cell_h = maximum(s[2] for s in spans) + padding
    cols = max(1, ceil(Int, sqrt(n)))
    out = PrintablePart[]
    sizehint!(out, n)
    for (i, p) in enumerate(parts)
        col = (i - 1) % cols
        row = (i - 1) ÷ cols
        dx = col * cell_w
        dy = row * cell_h
        T = eltype(eltype(p.mesh.position))
        new_pos = [Point3{T}(v[1] + T(dx), v[2] + T(dy), v[3]) for v in p.mesh.position]
        new_mesh = GeometryBasics.Mesh(
            new_pos,
            [TriangleFace{Int}(Int(f[1]), Int(f[2]), Int(f[3])) for f in p.mesh.faces],
        )
        push!(out, PrintablePart(new_mesh, p.color, p.name, p.face_colors))
    end
    out
end

# Drop vertices not referenced by any face; renumber faces accordingly.
function _compact_vertices(mesh::GeometryBasics.Mesh)
    used = falses(length(mesh.position))
    for f in mesh.faces
        used[Int(f[1])] = true
        used[Int(f[2])] = true
        used[Int(f[3])] = true
    end
    remap = zeros(Int, length(used))
    new_pos = eltype(mesh.position)[]
    for i in eachindex(used)
        if used[i]
            push!(new_pos, mesh.position[i])
            remap[i] = length(new_pos)
        end
    end
    new_faces = [TriangleFace{Int}(remap[Int(f[1])], remap[Int(f[2])], remap[Int(f[3])])
                 for f in mesh.faces]
    GeometryBasics.Mesh(new_pos, new_faces)
end

# ---------------------------------------------------------------------------
# Bond mesh: stepped cylinder with two pegs (one at each end) that
# friction-fit into the atom's sockets. Single watertight closed mesh.
# Bond order is encoded in the central-shaft radius (single bond uses the
# nominal `bond_radius`, double scales it up, triple scales further), so
# each printed bond piece always has exactly one peg per atom socket.
# ---------------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Build the printed mesh for one bond connecting atom surfaces at `p1` and
`p2`. The geometry is a **stepped cylinder**:

```
peg_r → shaft_r → peg_r       (radius profile)
└──╴   ────────   ╶──┘
peg ↑   ↑shaft↑   ↑ peg
```

i.e. a central cylindrical shaft of `bond_radius` (multiplied by an
order-dependent factor) running between `p1` and `p2`, plus a thinner
cylindrical peg of `peg_radius × peg_length` protruding from each shoulder
along the bond axis. The peg ends fit into atom sockets produced by
[`atom_sphere_mesh`](@ref). The whole part is a single watertight closed
mesh (six rings of vertices + two annular washers where the radius changes
+ two flat end disks).

`order = 1` uses `bond_radius` as-is. `order = 2` scales the shaft radius by
1.15 and `order = 3` by 1.3; pegs stay at `peg_radius` so every bond fits
the same atom socket regardless of order. (Parallel multi-cylinder visuals
would require CSG-merging at the shoulders.)
"""
function bond_cylinder_mesh(p1::AbstractVector{T}, p2::AbstractVector{T};
                            bond_radius::T = T(2.5),
                            peg_radius::T = T(1.5),
                            peg_length::T = T(8.0),
                            order::Int = 1,
                            segments::Int = 24) where T<:Real
    p1v = Vec3{T}(p1...)
    p2v = Vec3{T}(p2...)
    axis = p2v - p1v
    L = norm(axis)
    n, u, v = _ortho_frame(axis)

    # Order-dependent shaft radius (single peg always).
    shaft_r = bond_radius * (order == 1 ? T(1.0) : order == 2 ? T(1.15) : T(1.3))

    # Six ring profiles. Each ring has `segments` vertices.
    # Profile, from left peg tip to right peg tip:
    #   peg_tip → peg_base_at_surface → taper → shaft → taper → peg_base_at_surface → peg_tip
    # The conic taper sits in the visible bond gap (NOT inside the atom
    # socket, whose bore is straight). Matches the look of commercial
    # molecular-kit bond sticks (Molymod / HGS).
    # Cap chamfer at L/3 so short bonds still have a visible shaft.
    chamfer = min(shaft_r - peg_radius, L / T(3))
    chamfer = max(chamfer, T(0))
    rings = [
        (-peg_length,         peg_radius),   # 1: left peg tip (at socket bottom)
        (zero(T),             peg_radius),   # 2: peg ends at atom surface
        (chamfer,             shaft_r),      # 3: taper ends, shaft begins
        (L - chamfer,         shaft_r),      # 4: shaft ends, taper begins
        (L,                   peg_radius),   # 5: peg starts at far atom surface
        (L + peg_length,      peg_radius),   # 6: right peg tip
    ]

    verts = Point3{T}[]
    ring_start = Int[]   # first vertex index of each ring (1-based)
    for (z, r) in rings
        push!(ring_start, length(verts) + 1)
        center_z = p1v + z * n
        for k in 0:segments-1
            θ = T(2π) * k / segments
            push!(verts, Point3{T}((center_z + r * (cos(θ) * u + sin(θ) * v))...))
        end
    end

    faces = TriangleFace{Int}[]

    # Connect ring i to ring i+1 with a quad strip. The rings are oriented
    # CCW around +n (because the per-vertex angle θ_k = 2π k/segments goes
    # CCW in the right-handed (u, v) frame). For an outward-facing normal
    # (away from the axis), the panel's triangle must be wound so that
    # going (a, bnext, b) is CCW when viewed from OUTSIDE the cylinder.
    # `a` is at ring i, angle θ_k; `b` at ring i+1, angle θ_k; etc.
    for ri in 1:length(rings)-1
        a0 = ring_start[ri]
        b0 = ring_start[ri + 1]
        for k in 0:segments-1
            knext = mod(k + 1, segments)
            a     = a0 + k
            anext = a0 + knext
            b     = b0 + k
            bnext = b0 + knext
            push!(faces, TriangleFace{Int}(a, bnext, b))
            push!(faces, TriangleFace{Int}(a, anext, bnext))
        end
    end

    # End caps at ring 1 (left peg tip, normal -n) and ring 6 (right peg tip, +n).
    bot_centre_idx = length(verts) + 1
    push!(verts, Point3{T}((p1v + rings[1][1] * n)...))
    a0 = ring_start[1]
    for k in 0:segments-1
        knext = mod(k + 1, segments)
        push!(faces, TriangleFace{Int}(bot_centre_idx, a0 + knext, a0 + k))
    end

    top_centre_idx = length(verts) + 1
    push!(verts, Point3{T}((p1v + rings[end][1] * n)...))
    a0 = ring_start[end]
    for k in 0:segments-1
        knext = mod(k + 1, segments)
        push!(faces, TriangleFace{Int}(top_centre_idx, a0 + k, a0 + knext))
    end

    GeometryBasics.Mesh(verts, faces)
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
                          subdivisions::Int = 3,
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
        # Bond shoulders sit exactly at atom surfaces — the peg geometry
        # then extends INTO the atom's pre-drilled socket. Each printed
        # part is a separate solid; no volumetric overlap in the slicer.
        surf_i = ci + ri * n̂
        surf_j = cj - rj * n̂
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

    # Grid-layout the parts on the build plate. Each part has already been
    # settled to z = 0 with its xy bbox centred at the origin; the grid shifts
    # them to non-overlapping cells. Slicers that DO auto-arrange will simply
    # rearrange; slicers that don't (e.g. OrcaSlicer's 3MF importer) honour
    # this layout and the parts come in already spread out on the plate.
    _layout_grid(parts; padding = T(5.0))
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
