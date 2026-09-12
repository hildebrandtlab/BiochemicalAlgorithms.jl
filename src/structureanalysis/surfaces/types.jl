export
    Sphere,
    Circle3,
    AbstractMolecularSurface,
    icosphere,
    surface_area,
    nvertices,
    ntriangles,
    read_msms

"""
    Circle3{T<:Real}

Oriented circle in 3D: center `p`, unit normal `n` (the circle's axis), radius `r`.

`GeometryBasics.Circle` is 2D only; `Circle3` is used throughout the analytical
surface algorithms (Connolly) where probe-atom contact circles must carry an
orientation.
"""
struct Circle3{T<:Real}
    p::Vector3{T}
    n::Vector3{T}
    r::T
end

"""
    AbstractMolecularSurface{T<:Real}

Common supertype for analytical molecular surfaces (reduced surface, SAS, SES).
"""
abstract type AbstractMolecularSurface{T<:Real} end

@inline nvertices(m::GeometryBasics.Mesh) = length(m.position)
@inline ntriangles(m::GeometryBasics.Mesh) = length(m.faces)

"""
    $(TYPEDSIGNATURES)

Total surface area of a triangulated mesh, computed as `½ Σ ‖(b-a) × (c-a)‖`
over all triangles.
"""
function surface_area(m::GeometryBasics.Mesh)
    pos = m.position
    s = zero(eltype(eltype(pos)))
    @inbounds for f in m.faces
        a, b, c = pos[f[1]], pos[f[2]], pos[f[3]]
        s += norm(cross(b - a, c - a))
    end
    s / 2
end

"""
    extract_manifold(m::GeometryBasics.Mesh)

Post-pass that returns a maximal 2-manifold sub-mesh of `m`. For each
non-manifold edge (3+ incident triangles), keeps a pair with
consistent traversal (one direction a→b, other b→a) chosen by
largest dot product of triangle normals (= smoothest dihedral). Drops
the rest. Useful for simulation: the triangulator's advancing-front
occasionally over-fills edges with sliver triangles, and simulation
codes assume 2-manifold input.

Returns a fresh mesh with the kept triangle subset and unchanged
vertices/normals.
"""
function extract_manifold(m::GeometryBasics.Mesh)
    pos = m.position
    nfaces = length(m.faces)
    # Edge → list of (tri_idx, direction-bit) entries. Direction-bit: 1
    # if the triangle traverses the edge a→b (a<b sorted), 0 if b→a.
    EdgeKey = Tuple{Int, Int}
    edge_tris = Dict{EdgeKey, Vector{Tuple{Int, Bool}}}()
    for (ti, f) in enumerate(m.faces)
        for (x, y) in ((f[1], f[2]), (f[2], f[3]), (f[3], f[1]))
            a, b = min(x, y), max(x, y)
            forward = x == a
            push!(get!(edge_tris, (a, b), Tuple{Int, Bool}[]), (ti, forward))
        end
    end
    # Per-triangle normals (for dihedral scoring).
    tri_normal = Vector{eltype(pos)}(undef, nfaces)
    for (ti, f) in enumerate(m.faces)
        a = pos[f[1]]; b = pos[f[2]]; c = pos[f[3]]
        n = cross(b - a, c - a)
        ℓ = norm(n)
        tri_normal[ti] = ℓ > 0 ? n / ℓ : eltype(pos)(0, 0, 1)
    end
    # For each NM edge, pick the best pair of triangles to keep.
    keep = trues(nfaces)
    for (k, lst) in edge_tris
        length(lst) <= 2 && continue
        # Split into forward and backward sets.
        fwd = [t for t in lst if t[2]]
        bwd = [t for t in lst if !t[2]]
        # Score every (forward, backward) pair by normal dot product.
        # Higher = smoother dihedral.
        best_pair = (0, 0)
        best_score = -Inf
        for (tf, _) in fwd, (tb, _) in bwd
            score = dot(tri_normal[tf], tri_normal[tb])
            if score > best_score
                best_score = score
                best_pair = (tf, tb)
            end
        end
        # If no forward+backward pair exists, fall back to any two.
        if best_pair == (0, 0)
            best_pair = (lst[1][1], lst[2][1])
        end
        keepset = Set{Int}((best_pair[1], best_pair[2]))
        for (t, _) in lst
            t in keepset || (keep[t] = false)
        end
    end
    new_faces = [m.faces[i] for i in 1:nfaces if keep[i]]
    GeometryBasics.Mesh(m.position, new_faces; normal = m.normal)
end

"""
    fill_holes(m::GeometryBasics.Mesh; max_hole_size = 100)

Fill every boundary loop (closed cycle of count-1 edges) in `m` with
a triangle fan from a newly inserted centroid vertex. Returns a fresh
mesh where every edge is shared by ≥2 triangles (watertight on a
2-manifold input).

The triangle orientation of each fan is aligned with the local
outward normal estimated from the loop's existing boundary triangles.
`max_hole_size` caps the number of edges in any loop processed (very
large loops are usually a sign of corrupt input topology).

Use after `extract_manifold` for simulation-quality watertight
output. For boundary element / FEM codes that require closed
manifolds, this is the final mesh post-pass.
"""
function fill_holes(m::GeometryBasics.Mesh; max_hole_size::Int = 100)
    pos = m.position
    nrm = m.normal
    nfaces = length(m.faces)
    T = eltype(eltype(pos))
    # Build edge → triangle list (directed). For each triangle's edge
    # (x → y), record direction. A boundary edge has only one direction.
    EdgeKey = Tuple{Int, Int}
    edge_dirs = Dict{EdgeKey, Vector{Tuple{Int, Bool}}}()
    for (ti, f) in enumerate(m.faces)
        for (x, y) in ((f[1], f[2]), (f[2], f[3]), (f[3], f[1]))
            a, b = min(x, y), max(x, y)
            forward = x == a
            push!(get!(edge_dirs, (a, b), Tuple{Int, Bool}[]), (ti, forward))
        end
    end
    # Collect boundary edges (count == 1). Store with their ORIGINAL
    # traversal direction (= the direction the single incident triangle
    # walks the edge). The "hole" loop walks the edge in the OPPOSITE
    # direction (manifold fill convention).
    boundary_edges = Dict{Int, Vector{Tuple{Int, Int}}}()  # vertex → list of (next_vertex, tri_idx)
    for (k, lst) in edge_dirs
        length(lst) == 1 || continue
        ti, fwd = lst[1]
        a, b = k
        # Fill direction is opposite of triangle direction.
        from, to = fwd ? (b, a) : (a, b)
        push!(get!(boundary_edges, from, Tuple{Int, Int}[]), (to, ti))
    end
    new_pos = collect(pos)
    new_nrm = collect(nrm)
    new_faces = collect(m.faces)
    visited = Set{Tuple{Int, Int}}()
    n_loops = 0
    for (from_start, _) in boundary_edges
        haskey(boundary_edges, from_start) || continue
        # Find an unvisited edge starting at from_start.
        for (to_start, ti_start) in boundary_edges[from_start]
            (from_start, to_start) in visited && continue
            # Walk the loop.
            loop = Int[from_start, to_start]
            push!(visited, (from_start, to_start))
            cur = to_start
            guard = 0
            closed = false
            while guard < max_hole_size && haskey(boundary_edges, cur)
                guard += 1
                # Pick an unvisited outgoing edge from cur.
                next_v = 0
                for (nv, _) in boundary_edges[cur]
                    (cur, nv) in visited && continue
                    next_v = nv; break
                end
                next_v == 0 && break
                push!(visited, (cur, next_v))
                if next_v == from_start
                    closed = true; break
                end
                push!(loop, next_v)
                cur = next_v
            end
            (closed && length(loop) >= 3) || continue
            # Compute centroid and outward normal hint.
            centroid_pt = zero(eltype(pos))
            navg = zero(eltype(nrm))
            for v in loop
                centroid_pt = centroid_pt + pos[v]
                navg = navg + nrm[v]
            end
            centroid_pt = centroid_pt / T(length(loop))
            ℓ = norm(navg)
            navg = ℓ > 0 ? navg / ℓ : eltype(nrm)(0, 0, 1)
            push!(new_pos, centroid_pt)
            push!(new_nrm, navg)
            c_idx = length(new_pos)
            # Fan triangles. Wind so the cross product agrees with navg.
            for k in 1:length(loop)
                a = loop[k]
                b = loop[mod1(k + 1, length(loop))]
                pa = pos[a]; pb = pos[b]
                n_tri = cross(pb - pa, centroid_pt - pa)
                if dot(n_tri, navg) >= 0
                    push!(new_faces, GeometryBasics.TriangleFace{Int}(a, b, c_idx))
                else
                    push!(new_faces, GeometryBasics.TriangleFace{Int}(a, c_idx, b))
                end
            end
            n_loops += 1
        end
    end
    GeometryBasics.Mesh(new_pos, new_faces; normal = new_nrm)
end

"""
    make_watertight(m::GeometryBasics.Mesh; tol = 1e-3, max_hole_size = 100)

Convenience pipeline for simulation-quality watertight meshes:
1. `weld_close_vertices` to merge coincident duplicates.
2. `extract_manifold` to eliminate non-manifold edges.
3. `fill_holes` to close all remaining boundary loops.

Returns a closed 2-manifold mesh suitable for boundary element /
finite element / volumetric simulation codes that assume watertight
input.
"""
function make_watertight(m::GeometryBasics.Mesh;
                         tol::Real = 1e-3, max_hole_size::Int = 100)
    m1 = weld_close_vertices(m; tol = tol)
    m2 = extract_manifold(m1)
    fill_holes(m2; max_hole_size = max_hole_size)
end

"""
    weld_close_vertices(m::GeometryBasics.Mesh; tol = 1e-3)

Post-pass that merges vertex pairs whose Euclidean distance is below
`tol` (default 0.001 Å). Triangles whose 3 vertices collapse to fewer
than 3 distinct points after welding are dropped (degenerate slivers
— zero area). Returns a fresh `GeometryBasics.Mesh` with welded
vertex pool.
"""
function weld_close_vertices(m::GeometryBasics.Mesh; tol::Real = 1e-3)
    pos = m.position
    nrm = m.normal
    T = eltype(eltype(pos))
    tol2 = T(tol)^2
    n = length(pos)
    # Spatial grid hash: cell size = tol. Each vertex maps to a (gx, gy, gz)
    # cell; neighbours within tol live in this or adjacent cells.
    cell_size = T(tol)
    grid = Dict{NTuple{3, Int}, Vector{Int}}()
    canonical = collect(1:n)  # canonical[i] = representative vertex for vertex i
    for i in 1:n
        p = pos[i]
        gx = floor(Int, p[1] / cell_size)
        gy = floor(Int, p[2] / cell_size)
        gz = floor(Int, p[3] / cell_size)
        found = 0
        # Check this cell and 26 neighbours.
        for dz in -1:1, dy in -1:1, dx in -1:1
            cell = (gx + dx, gy + dy, gz + dz)
            haskey(grid, cell) || continue
            for j in grid[cell]
                jp = pos[j]
                d2 = (p[1] - jp[1])^2 + (p[2] - jp[2])^2 + (p[3] - jp[3])^2
                if d2 < tol2
                    found = j
                    break
                end
            end
            found != 0 && break
        end
        if found != 0
            canonical[i] = canonical[found]
        else
            push!(get!(grid, (gx, gy, gz), Int[]), i)
        end
    end
    # Build the new vertex/normal list (unique representatives only).
    remap = zeros(Int, n)
    new_pos = eltype(pos)[]
    new_nrm = eltype(nrm)[]
    for i in 1:n
        if canonical[i] == i
            push!(new_pos, pos[i])
            push!(new_nrm, nrm[i])
            remap[i] = length(new_pos)
        end
    end
    for i in 1:n
        if canonical[i] != i
            remap[i] = remap[canonical[i]]
        end
    end
    # Rebuild face list, dropping degenerate triangles.
    new_faces = eltype(m.faces)[]
    for f in m.faces
        a = remap[f[1]]; b = remap[f[2]]; c = remap[f[3]]
        (a == b || b == c || a == c) && continue
        push!(new_faces, GeometryBasics.TriangleFace{Int}(a, b, c))
    end
    GeometryBasics.Mesh(new_pos, new_faces; normal = new_nrm)
end

"""
    icosphere(::Type{T<:Real}, subdivisions::Int = 0) -> GeometryBasics.Mesh

Generate a geodesic unit sphere by `subdivisions` rounds of icosahedral
midpoint refinement.

Starts from a regular icosahedron (12 vertices, 20 faces), subdivides each
triangle into four by splitting at edge midpoints, and projects every new
vertex back onto the unit sphere. The resulting mesh has
`10·4^subdivisions + 2` vertices.

Returns a `GeometryBasics.Mesh` with `Point3{T}` positions, per-vertex unit
`Vec3{T}` outward normals (equal to the position on the unit sphere), and
`TriangleFace{Int}` faces.
"""
function icosphere(::Type{T}, subdivisions::Int=0) where {T<:Real}
    # BALL `TriangulatedSphere::icosaeder` (triangulatedSurface.C:1395-1424)
    # uses a POLAR-AXIS icosahedron: vertex 0 = north pole (0,0,1),
    # vertex 11 = south pole (0,0,-1), 5 vertices each at z = ±0.4472.
    # My port previously used the (rotated) golden-ratio coordinates;
    # after clipping against the contact-face SES-edge planes the
    # surviving icosphere subsets differed from BALL's, which fed extra
    # candidates into the advancing-front and inflated the triangle count
    # by ~2.85× on BPTI.
    # Vertex positions are taken byte-for-byte from BALL.
    raw = Point3{T}[
        Point3{T}( 0.0     ,  0.0     ,  1.0      ),
        Point3{T}( 0.894427,  0.0     ,  0.4472135),
        Point3{T}( 0.276393,  0.850651,  0.4472135),
        Point3{T}(-0.723607,  0.525731,  0.4472135),
        Point3{T}(-0.723607, -0.525731,  0.4472135),
        Point3{T}( 0.276393, -0.850651,  0.4472135),
        Point3{T}( 0.723607,  0.525731, -0.4472135),
        Point3{T}(-0.276393,  0.850651, -0.4472135),
        Point3{T}(-0.894427,  0.0     , -0.4472135),
        Point3{T}(-0.276393, -0.850651, -0.4472135),
        Point3{T}( 0.723607, -0.525731, -0.4472135),
        Point3{T}( 0.0     ,  0.0     , -1.0      ),
    ]
    vertices = [v / norm(v) for v in raw]
    # Triangles derived from BALL's `triangles_tmp` definitions. Each row is
    # a CCW-from-outside triple of vertex indices (1-based) so the cross
    # product `(b-a) × (c-a)` of the triangle (a, b, c) gives an outward
    # normal. The BALL list had 9 of the 20 base triangles wound CW (mostly
    # the lower-middle band and the bottom cap); we corrected those by
    # swapping the last two indices so signed-volume of the icosphere
    # comes out positive (~ 4π/3 for the unit sphere).
    faces = TriangleFace{Int}[
        ( 1,  2,  3), ( 1,  3,  4), ( 1,  4,  5), ( 1,  5,  6), ( 1,  6,  2),
        ( 2,  7,  3), ( 3,  7,  8), ( 3,  8,  4), ( 4,  8,  9), ( 4,  9,  5),
        ( 5,  9, 10), ( 5, 10,  6), ( 6, 10, 11), ( 2, 11,  7), ( 2,  6, 11),
        ( 7, 12,  8), ( 8, 12,  9), ( 9, 12, 10), (10, 12, 11), ( 7, 11, 12),
    ]

    for _ in 1:subdivisions
        midpoint_cache = Dict{NTuple{2,Int},Int}()
        new_faces = TriangleFace{Int}[]
        sizehint!(new_faces, 4 * length(faces))

        get_midpoint = (a::Int, b::Int) -> begin
            key = a < b ? (a, b) : (b, a)
            haskey(midpoint_cache, key) && return midpoint_cache[key]
            m = (vertices[a] + vertices[b]) / 2
            m = m / norm(m)
            push!(vertices, m)
            idx = length(vertices)
            midpoint_cache[key] = idx
            idx
        end

        for f in faces
            v1, v2, v3 = Int(f[1]), Int(f[2]), Int(f[3])
            a = get_midpoint(v1, v2)
            b = get_midpoint(v2, v3)
            c = get_midpoint(v3, v1)
            push!(new_faces, TriangleFace{Int}(v1, a, c))
            push!(new_faces, TriangleFace{Int}(v2, b, a))
            push!(new_faces, TriangleFace{Int}(v3, c, b))
            push!(new_faces, TriangleFace{Int}(a, b, c))
        end
        faces = new_faces
    end

    normals = [Vec3{T}(v...) for v in vertices]
    GeometryBasics.Mesh(vertices, faces; normal=normals)
end

"""
    read_msms(
        ::Type{T<:Real},
        vert_path::AbstractString,
        face_path::AbstractString
    ) -> GeometryBasics.Mesh

Read a triangulated surface from a pair of MSMS files (Sanner et al.).

`vert_path` provides per-vertex coordinates and outward normals; `face_path`
provides 1-based triangle vertex indices. Both files use a three-line header
followed by whitespace-separated data lines. Returns a `GeometryBasics.Mesh`.
"""
function read_msms(::Type{T}, vert_path::AbstractString, face_path::AbstractString) where {T<:Real}
    vertices = Point3{T}[]
    normals  = Vec3{T}[]

    open(vert_path, "r") do io
        for _ in 1:3
            eof(io) && break
            readline(io)
        end
        for line in eachline(io)
            fields = split(line)
            length(fields) >= 6 || continue
            push!(vertices, Point3{T}(parse(T, fields[1]), parse(T, fields[2]), parse(T, fields[3])))
            push!(normals,  Vec3{T}(parse(T, fields[4]),  parse(T, fields[5]),  parse(T, fields[6])))
        end
    end

    faces = TriangleFace{Int}[]
    open(face_path, "r") do io
        for _ in 1:3
            eof(io) && break
            readline(io)
        end
        for line in eachline(io)
            fields = split(line)
            length(fields) >= 3 || continue
            push!(faces, TriangleFace{Int}(parse(Int, fields[1]), parse(Int, fields[2]), parse(Int, fields[3])))
        end
    end

    GeometryBasics.Mesh(vertices, faces; normal=normals)
end
