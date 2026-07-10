export
    ReducedSurface,
    compute_reduced_surface

# ---------------------------------------------------------------------------
# Records.
#
# The reduced surface is built by direct enumeration: every (atom_a, atom_b,
# atom_c) triple in the neighbor list yields up to two candidate probe
# positions, each of which is a candidate RS face. A candidate is kept iff
# the probe does not overlap any other atom (Connolly's validity criterion).
# Edges and vertices are derived from the resulting face set.
# ---------------------------------------------------------------------------

"""
    RSVertex

Reduced-surface vertex sitting on a single atom. Carries the indices of
the incident edges and faces (graph back-references).
"""
mutable struct RSVertex
    atom::Int               # 1-based atom index in source container
    edges::Set{Int}
    faces::Set{Int}
end
RSVertex(atom::Int) = RSVertex(atom, Set{Int}(), Set{Int}())

"""
    RSEdge{T}

Reduced-surface edge connecting two RS vertices (on two distinct atoms). The
torus geometry describes the path swept by the probe center while rolling
over the two atoms between the two adjacent faces.
"""
mutable struct RSEdge{T<:Real}
    v1::Int
    v2::Int
    f1::Int                 # face indices; 0 = none
    f2::Int
    center_of_torus::Vector3{T}
    radius_of_torus::T
    angle::T                # rolling angle between f1 and f2 (0 if f2 = 0)
    contact_circle1::Circle3{T}   # on atom of v1
    contact_circle2::Circle3{T}   # on atom of v2
    singular::Bool
    # Singular-edge geometry: when the probe sphere intersects the (a,b)
    # line, two extra "singular vertices" appear on the SES toric face.
    # These are the line/probe-sphere intersection points, ordered with
    # `intersection1` closer to atom of v1. Populated only when `singular`.
    intersection1::Vector3{T}
    intersection2::Vector3{T}
    # BALL-faithful index assignment. BALL's `rs_->insert(edge)`
    # (reducedSurface.C:250) assigns `edge->index_ = number_of_edges_++`
    # at the moment `treatEdge` succeeds (= the second face attaches).
    # My port creates edges eagerly during rolling, so we track the
    # BALL-order index here and PERMUTE `rs.edges` at the end of RS
    # construction to match. `ball_idx == 0` means "not yet assigned"
    # (= edge was created but not treated, or is a deleted-sentinel).
    ball_idx::Int
end
function RSEdge{T}(v1::Int, v2::Int) where {T<:Real}
    z = zero(Vector3{T})
    c = Circle3{T}(z, z, zero(T))
    RSEdge{T}(v1, v2, 0, 0, z, zero(T), zero(T), c, c, false, z, z, 0)
end

"""
    RSFace{T}

Reduced-surface face: a triangle of three vertices (on three distinct atoms)
plus the probe sphere center and outward normal.
"""
mutable struct RSFace{T<:Real}
    v1::Int; v2::Int; v3::Int
    # Incident RS edges. BALL stores these in `std::list<Edge*>` so a face
    # can have more than 3 incident edges when multiple probes share the
    # same atom-pair axis (the "stacked probes on the same torus" case —
    # roughly 260 instances on BPTI). My Julia version uses `Vector{Int}`
    # for the same reason: a fixed e1/e2/e3 trio loses one edge per slot
    # in that case and breaks downstream SES topology.
    edges::Vector{Int}
    center::Vector3{T}      # probe sphere center
    normal::Vector3{T}      # outward normal
    singular::Bool
end

# ---------------------------------------------------------------------------
# ReducedSurface container
# ---------------------------------------------------------------------------

"""
    ReducedSurface{T}

Connolly's reduced surface for a set of atom spheres and a probe radius.

The RS is a graph of vertices (per-atom), edges (per atom pair where the
probe can roll between two adjacent faces), and faces (per atom triple where
the probe rests on three atoms simultaneously).

Build it from an atom container via [`compute_reduced_surface`](@ref); atoms
must have `radius > 0` (use [`assign_radii!`](@ref) first).

# Fields
 - `atoms::Vector{Sphere{T}}`  — input atom spheres (radius does *not* include the probe)
 - `probe_radius::T`
 - `vertices::Vector{RSVertex}`, `edges::Vector{RSEdge{T}}`, `faces::Vector{RSFace{T}}`
 - `bounding_box::BoundingBox{T}` — over atom centers
 - `r_max::T` — maximum atom radius
"""
mutable struct ReducedSurface{T<:Real} <: AbstractMolecularSurface{T}
    atoms::Vector{Sphere{T}}
    probe_radius::T
    vertices::Vector{RSVertex}
    edges::Vector{RSEdge{T}}
    faces::Vector{RSFace{T}}
    bounding_box::BoundingBox{T}
    r_max::T
end

# Counts only ALIVE elements — `atom == 0` (vertices), `v1 == 0`
# (edges/faces) are deleted sentinels left in place to keep stable
# indices for the SES/triangulator.
@inline nvertices(rs::ReducedSurface) = count(v -> v.atom != 0, rs.vertices)
@inline nedges(rs::ReducedSurface)    = count(e -> e.v1 != 0, rs.edges)
@inline nfaces(rs::ReducedSurface)    = count(f -> f.v1 != 0, rs.faces)

# ---------------------------------------------------------------------------
# Type conversion (used by the F32→F64-internal promotion path; see
# `compute_reduced_surface(::AbstractAtomContainer{T})`).
# ---------------------------------------------------------------------------

function _convert_rs(::Type{T_out}, rs::ReducedSurface{T_in}) where {T_in<:Real, T_out<:Real}
    T_in === T_out && return rs
    atoms = Sphere{T_out}[
        Sphere{T_out}(Point3{T_out}(s.center...), T_out(s.r)) for s in rs.atoms
    ]
    pr = T_out(rs.probe_radius)
    edges = RSEdge{T_out}[
        RSEdge{T_out}(e.v1, e.v2, e.f1, e.f2,
                       Vector3{T_out}(T_out.(e.center_of_torus)...),
                       T_out(e.radius_of_torus),
                       T_out(e.angle),
                       Circle3{T_out}(Vector3{T_out}(T_out.(e.contact_circle1.p)...),
                                      Vector3{T_out}(T_out.(e.contact_circle1.n)...),
                                      T_out(e.contact_circle1.r)),
                       Circle3{T_out}(Vector3{T_out}(T_out.(e.contact_circle2.p)...),
                                      Vector3{T_out}(T_out.(e.contact_circle2.n)...),
                                      T_out(e.contact_circle2.r)),
                       e.singular,
                       Vector3{T_out}(T_out.(e.intersection1)...),
                       Vector3{T_out}(T_out.(e.intersection2)...),
                       e.ball_idx)
        for e in rs.edges
    ]
    faces = RSFace{T_out}[
        RSFace{T_out}(f.v1, f.v2, f.v3, copy(f.edges),
                       Vector3{T_out}(T_out.(f.center)...),
                       Vector3{T_out}(T_out.(f.normal)...),
                       f.singular)
        for f in rs.faces
    ]
    bb = BoundingBox{T_out}(Vector3{T_out}(T_out.(rs.bounding_box.min)...),
                              Vector3{T_out}(T_out.(rs.bounding_box.max)...))
    ReducedSurface{T_out}(atoms, pr, rs.vertices, edges, faces, bb, T_out(rs.r_max))
end

@inline function _sort3(a::Int, b::Int, c::Int)
    if a > b; a, b = b, a; end
    if a > c; a, c = c, a; end
    if b > c; b, c = c, b; end
    (a, b, c)
end

@inline _sort2(a::Int, b::Int) = a < b ? (a, b) : (b, a)

# ---------------------------------------------------------------------------
# Neighbor-list pre-processing.
# ---------------------------------------------------------------------------

# Returns (neighbors_per_atom, bounding_box, r_max). `neighbors[i]` lists atom
# indices j (sorted ascending) such that the probe could simultaneously touch
# both i and j: |c_i - c_j| ≤ r_i + r_j + 2·probe_radius.
function _build_neighbors(atoms::Vector{Sphere{T}}, probe_radius::T) where T
    n = length(atoms)
    if n == 0
        return (Vector{Vector{Int}}(), BoundingBox{T}(zero(Vector3{T}), zero(Vector3{T})), zero(T))
    end

    lo = Vector3{T}(atoms[1].center...)
    hi = Vector3{T}(atoms[1].center...)
    rmax = atoms[1].r
    for i in 2:n
        c = atoms[i].center
        lo = Vector3{T}(min(lo[1], c[1]), min(lo[2], c[2]), min(lo[3], c[3]))
        hi = Vector3{T}(max(hi[1], c[1]), max(hi[2], c[2]), max(hi[3], c[3]))
        rmax = max(rmax, atoms[i].r)
    end
    bb = BoundingBox{T}(lo, hi)
    cutoff = 2 * (rmax + probe_radius)
    nb = [Int[] for _ in 1:n]
    if n >= 2
        centers = [Vector3{T}(atoms[i].center...) for i in 1:n]
        pairs = neighborlist(; positions=centers, cutoff=cutoff)
        for (i, j, _) in pairs
            d = centers[i] - centers[j]
            ri = atoms[i].r + probe_radius
            rj = atoms[j].r + probe_radius
            if dot(d, d) <= (ri + rj)^2
                push!(nb[i], j)
                push!(nb[j], i)
            end
        end
        for ns in nb; sort!(ns); end
    end
    (nb, bb, rmax)
end

# Mutual neighbors of two atoms (sorted ascending).
function _common_neighbors(nb::Vector{Vector{Int}}, a::Int, b::Int)
    na = nb[a]; nbb = nb[b]
    out = Int[]
    i = 1; j = 1
    while i <= length(na) && j <= length(nbb)
        x = na[i]; y = nbb[j]
        if x == y
            x != a && x != b && push!(out, x)
            i += 1; j += 1
        elseif x < y
            i += 1
        else
            j += 1
        end
    end
    out
end

# ---------------------------------------------------------------------------
# Per-atom lobe decomposition. Returns `lobes::Vector{Dict{Int,Int}}` where
# `lobes[a][face_idx] = lobe_id` for every face that touches atom `a`. Two
# faces fall into the same lobe iff they are connected by a rolling-probe
# edge sharing atom `a`.
#
# The edge graph mirrors BALL's RS construction: for each face f with
# atom-pair side (a, X), f connects to the face whose probe sits at the
# smallest POSITIVE rotation angle in f's rolling-forward direction around
# (a, X). The forward direction is signed per face by `f.normal × axis`
# tested against the third atom (reducedSurface.C:1099-1106).
#
# When the n faces around a given (a, X) axis split into multiple rolling
# loops that close on themselves, faces at atom `a` decompose into separate
# components -> multiple RSVertex lobes for atom `a`.
# ---------------------------------------------------------------------------
function _build_atom_lobes(face_keys::Vector{NTuple{4,Int}},
                           face_probes::Vector{Vector3{T}},
                           face_normals::Vector{Vector3{T}},
                           atoms::Vector{Sphere{T}}) where T
    n_atoms = length(atoms)
    faces_at_atom = [Int[] for _ in 1:n_atoms]
    for (fi, fk) in enumerate(face_keys)
        a, b, c, _ = fk
        push!(faces_at_atom[a], fi)
        push!(faces_at_atom[b], fi)
        push!(faces_at_atom[c], fi)
    end

    # Group faces by atom pair.
    pair_to_faces = Dict{NTuple{2,Int}, Vector{Int}}()
    for (fi, fk) in enumerate(face_keys)
        a, b, c, _ = fk
        for (x, y) in ((a, b), (a, c), (b, c))
            key = _sort2(x, y)
            push!(get!(() -> Int[], pair_to_faces, key), fi)
        end
    end

    # For each (a, b), compute next_face[k] for k in 1..n using each
    # face's own rolling-forward direction. Returns the connection set:
    # `pair_links[(a, b)]` is a Vector{Tuple{Int,Int}} of (fi, fj) UNDIRECTED
    # links between faces.
    pair_links = Dict{NTuple{2,Int}, Vector{Tuple{Int,Int}}}()
    for (key, faces) in pair_to_faces
        n = length(faces)
        n < 2 && continue
        a, b = key
        pa = Vector3{T}(atoms[a].center...)
        pb = Vector3{T}(atoms[b].center...)
        next_face = zeros(Int, n)
        for k in 1:n
            fk = face_keys[faces[k]]
            ak, bk, ck, _ = fk
            third_k = (ak != a && ak != b) ? ak :
                      (bk != a && bk != b) ? bk : ck
            ax = pa - pb
            test_vec = cross(face_normals[faces[k]], ax)
            third_p = Vector3{T}(atoms[third_k].center...)
            if dot(test_vec, third_p) > dot(test_vec, pa)
                ax = -ax
            end
            ax_len2 = dot(ax, ax)
            rel_k = face_probes[faces[k]] - pa
            proj_k = rel_k - (dot(rel_k, ax) / ax_len2) * ax
            best_angle = T(3 * π)
            best_j = 0
            for j in 1:n
                j == k && continue
                rel_j = face_probes[faces[j]] - pa
                proj_j = rel_j - (dot(rel_j, ax) / ax_len2) * ax
                ang = oriented_angle(proj_k, proj_j, ax)
                if ang > T(1e-6) && ang < best_angle
                    best_angle = ang
                    best_j = j
                end
            end
            next_face[k] = best_j
        end
        links = Tuple{Int,Int}[]
        for k in 1:n
            nk = next_face[k]
            nk == 0 && continue
            fi1 = faces[k]; fi2 = faces[nk]
            push!(links, fi1 < fi2 ? (fi1, fi2) : (fi2, fi1))
        end
        unique!(links)
        pair_links[key] = links
    end

    # Per-atom lobe ids: faces are unioned only by (atom, X)-rolling links.
    atom_lobes = Vector{Dict{Int,Int}}(undef, n_atoms)
    for a in 1:n_atoms
        face_list = faces_at_atom[a]
        if isempty(face_list)
            atom_lobes[a] = Dict{Int,Int}()
            continue
        end
        # Local union-find over faces at atom a.
        parent = Dict{Int,Int}(fi => fi for fi in face_list)
        function lfind(x)
            while parent[x] != x
                parent[x] = parent[parent[x]]
                x = parent[x]
            end
            x
        end
        function lunite!(x, y)
            rx = lfind(x); ry = lfind(y)
            rx != ry && (parent[rx] = ry)
        end
        # Apply links through pairs that include atom a.
        for (pair, links) in pair_links
            (pair[1] == a || pair[2] == a) || continue
            for (fi, fj) in links
                lunite!(fi, fj)
            end
        end
        lobe_map = Dict{Int,Int}()
        lobe_id_for_root = Dict{Int,Int}()
        for fi in face_list
            root = lfind(fi)
            if !haskey(lobe_id_for_root, root)
                lobe_id_for_root[root] = length(lobe_id_for_root) + 1
            end
            lobe_map[fi] = lobe_id_for_root[root]
        end
        atom_lobes[a] = lobe_map
    end
    atom_lobes
end

# ---------------------------------------------------------------------------
# Probe validity: a probe center is valid iff it does not overlap any of the
# three atoms beyond the tangent surface (already guaranteed by construction)
# AND no fourth atom occludes it.
# ---------------------------------------------------------------------------

# Returns true if `probe_center` is *not* inside any neighbour of {a, b, c}.
function _probe_unoccluded(probe_center::Vector3{T}, atoms::Vector{Sphere{T}},
                           probe_radius::T, neighbors::Vector{Vector{Int}},
                           a::Int, b::Int, c::Int, eps::T) where T
    # Candidates: atoms in nbrs(a) ∩ nbrs(b) ∩ nbrs(c), excluding {a,b,c}.
    nab = _common_neighbors(neighbors, a, b)
    isempty(nab) && return true
    nac = _common_neighbors(neighbors, a, c)
    isempty(nac) && return true
    # Intersect the two sorted lists, excluding b/c/a respectively.
    i = 1; j = 1
    while i <= length(nab) && j <= length(nac)
        x = nab[i]; y = nac[j]
        if x == y
            if x != a && x != b && x != c
                rj = atoms[x].r + probe_radius
                d = probe_center - Vector3{T}(atoms[x].center...)
                if dot(d, d) < rj * rj - eps
                    return false
                end
            end
            i += 1; j += 1
        elseif x < y
            i += 1
        else
            j += 1
        end
    end
    true
end

# ---------------------------------------------------------------------------
# Geometric helpers tailored to the RS algorithm.
# ---------------------------------------------------------------------------

# Outward face normal for the triangle (a, b, c) given the probe center.
function _face_normal(rs::ReducedSurface{T}, a::Int, b::Int, c::Int,
                      probe::Vector3{T}) where T
    pa = rs.atoms[a].center; pb = rs.atoms[b].center; pc = rs.atoms[c].center
    n = plane_normal(pa, pb, pc)
    n === nothing && return Vector3{T}(0, 0, 0)
    if dot(n, probe - Vector3{T}(pa[1], pa[2], pa[3])) > 0
        return Vector3{T}(-n[1], -n[2], -n[3])
    end
    n
end

# Contact circles for the probe rolling over two atoms.
function _contact_circles(rs::ReducedSurface{T}, a::Int, b::Int) where T
    pr = rs.probe_radius
    sa = Sphere{T}(rs.atoms[a].center, rs.atoms[a].r + pr)
    sb = Sphere{T}(rs.atoms[b].center, rs.atoms[b].r + pr)
    c1 = intersect_spheres(sa, sb)
    c1 === nothing && return nothing
    ratio_a = rs.atoms[a].r / sa.r
    ratio_b = rs.atoms[b].r / sb.r
    ca = Vector3{T}(sa.center...)
    cb = Vector3{T}(sb.center...)
    p2 = ca + (c1.p - ca) * ratio_a
    p3 = cb + (c1.p - cb) * ratio_b
    c2 = Circle3{T}(p2, c1.n, c1.r * ratio_a)
    c3 = Circle3{T}(p3, c1.n, c1.r * ratio_b)
    (c1, c2, c3)
end

# ---------------------------------------------------------------------------
# Top-level driver: enumerate faces, then derive edges and vertices.
# ---------------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Compute Connolly's reduced surface for an atom container.

The reduced surface (RS) is the topological skeleton produced by a spherical
probe of radius `probe_radius` rolling over the inflated atom spheres: its
vertices sit on individual atoms, its edges on pairs of atoms (where the
probe rolls between two adjacent faces), and its faces on triples of atoms
(where the probe rests tangent to all three).

Atoms with `radius == 0` are skipped; call [`assign_radii!`](@ref) first to
populate van-der-Waals radii.

# Precision

The rolling-probe algorithm has many ULP-sensitive geometric tests
(3-sphere intersection, probe-collision checks, oriented angles).
BALL hardcodes `double` precision throughout its surface stack —
all its `TSphere3`, `TVector3`, `TCircle3` are `<double>` regardless
of atom storage. We follow that convention: when called with a
`Float32` atom container the RS is built internally in `Float64`,
then **converted to `ReducedSurface{Float32}` for return**. This
keeps the user-facing type stable while avoiding the ~50 hole
defects that Float32-internal arithmetic produces on BPTI.

# Keyword arguments
 - `probe_radius::T = 1.5`  — probe sphere radius (Å)
"""
function compute_reduced_surface(ac::AbstractAtomContainer{T};
                                 probe_radius::T = T(1.5)) where T
    if T !== Float64
        # Promote to Float64 for the rolling-probe construction (BALL's
        # convention), then convert back. _atom_spheres_for_rs reads from
        # the AtomContainer; we widen each sphere here.
        spheres64 = Sphere{Float64}[
            Sphere{Float64}(Point3{Float64}(s.center...), Float64(s.r))
            for s in _atom_spheres_for_rs(ac)
        ]
        rs64 = compute_reduced_surface(spheres64;
                                         probe_radius=Float64(probe_radius))
        return _convert_rs(T, rs64)
    end
    compute_reduced_surface(_atom_spheres_for_rs(ac); probe_radius)
end

function compute_reduced_surface(spheres::Vector{Sphere{T}};
                                 probe_radius::Real = 1.5) where T
    if T !== Float64
        spheres64 = Sphere{Float64}[
            Sphere{Float64}(Point3{Float64}(s.center...), Float64(s.r))
            for s in spheres
        ]
        rs64 = compute_reduced_surface(spheres64; probe_radius=Float64(probe_radius))
        return _convert_rs(T, rs64)
    end
    pr = T(probe_radius)
    # BALL-faithful path: rolling-probe construction mirroring
    # `RSComputer::run` (reducedSurface.C:578-700). Degeneracies are
    # corrected in-place via `_rs_correct!`; no outer restart loop.
    return _rs_compute_with_rolling(spheres, pr)
end

# Detect BALL's degeneracy conditions in thirdAtom (reducedSurface.C:1081):
# (1) the "PROBE TOUCHES FOUR ATOMS" case where some 4th atom is at the
# probe's surface (= rotation angle 0 or 2π), OR
# (2) the "MULTIPLE ATOMS TIE" case where 2+ atoms have the same minimum
# rotation angle around the rolling axis (line 1157-1168).
# In both cases, BALL calls correct(atom) on the offending atom. We detect
# these by enumerating atom triples and pairs, looking for the same
# probe-surface or co-circular conditions.
function _find_four_atom_degeneracy(atoms::Vector{Sphere{T}}, pr::T) where T
    isempty(atoms) && return 0
    # Build neighbor list (atoms within 2·(r_max + pr)) using my existing
    # logic. We only need atoms that COULD touch each other through a probe.
    neighbors, _, _ = _build_neighbors(atoms, pr)
    # BALL's `Maths::isZero` epsilon is 1e-6 for the rotation-angle
    # comparison. Use the same on the squared-distance check.
    surface_tol = T(1e-6)
    for a in 1:length(atoms)
        for b in neighbors[a]
            b <= a && continue
            common = _common_neighbors(neighbors, a, b)
            for c in common
                c <= b && continue
                s1 = Sphere{T}(atoms[a].center, atoms[a].r + pr)
                s2 = Sphere{T}(atoms[b].center, atoms[b].r + pr)
                s3 = Sphere{T}(atoms[c].center, atoms[c].r + pr)
                pts = intersect_three_spheres(s1, s2, s3)
                pts === nothing && continue
                for probe in pts
                    # Look for 4th atom d such that
                    # |probe - atom_d.center|² ≈ (atom_d.r + pr)²
                    # (i.e., atom d is essentially at the probe surface).
                    nab = _common_neighbors(neighbors, a, b)
                    nac = _common_neighbors(neighbors, a, c)
                    i = 1; j = 1
                    while i <= length(nab) && j <= length(nac)
                        x = nab[i]; y = nac[j]
                        if x == y
                            if x != a && x != b && x != c
                                rj = atoms[x].r + pr
                                d = probe - Vector3{T}(atoms[x].center...)
                                d_sq = dot(d, d)
                                if abs(d_sq - rj * rj) < surface_tol
                                    return x   # atom to shrink
                                end
                            end
                            i += 1; j += 1
                        elseif x < y
                            i += 1
                        else
                            j += 1
                        end
                    end
                end
            end
        end
    end
    return 0
end

# Apply accumulated radius corrections.
function _apply_radius_corrections(spheres::Vector{Sphere{T}},
                                     shrunk::Dict{Int, T}) where T
    isempty(shrunk) && return spheres
    out = Vector{Sphere{T}}(undef, length(spheres))
    for (i, s) in enumerate(spheres)
        delta = get(shrunk, i, zero(T))
        out[i] = delta == zero(T) ? s :
                 Sphere{T}(s.center, max(s.r - delta, T(0.001)))
    end
    return out
end

# Detect a stacked-probe configuration: 2 RS edges sharing the same atom
# pair on the same RS face. Returns the index (into spheres) of the
# "third atom" of one of these edges that should be shrunk, or 0 if none.
# Mirrors BALL's correct(atom) target: the atom whose presence causes
# the degenerate rolling configuration.
function _detect_stacked_probe_atom(rs::ReducedSurface{T}) where T
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        # For each atom-pair side of the face, count incident RS edges.
        verts = (f.v1, f.v2, f.v3)
        atoms_v = ntuple(k -> rs.vertices[verts[k]].atom, 3)
        for k in 1:3
            a = atoms_v[k]; b = atoms_v[mod1(k + 1, 3)]
            matches = Int[]
            for ei in f.edges
                re = rs.edges[ei]
                a1 = rs.vertices[re.v1].atom
                a2 = rs.vertices[re.v2].atom
                if (a1 == a && a2 == b) || (a1 == b && a2 == a)
                    push!(matches, ei)
                end
            end
            if length(matches) >= 2
                # Pick the THIRD atom of one of the secondary edges. The
                # "primary" edge survives (it owns the slot); the
                # secondary edge's other RS face has a third atom we
                # shrink to eliminate the degenerate stacking.
                e_secondary = matches[2]
                re = rs.edges[e_secondary]
                # Find the OTHER RS face on this edge (not fi).
                other_face_idx = re.f1 == fi ? re.f2 : re.f1
                other_face_idx == 0 && continue
                other_f = rs.faces[other_face_idx]
                other_f.v1 == 0 && continue
                # Third atom of other_f (not a, not b).
                for v in (other_f.v1, other_f.v2, other_f.v3)
                    at = rs.vertices[v].atom
                    if at != a && at != b
                        return at
                    end
                end
            end
        end
    end
    return 0
end

# =========================================================================
# Port of BALL's rolling-probe RS construction (RSComputer in
# reducedSurface.C). The enumeration path below remains as the active
# implementation; this rolling-probe path is being developed in parallel
# (entry point: `_compute_reduced_surface_rolling`).
#
# Key correspondence with BALL:
#   * `preProcessing` — build neighbor relationships
#   * `findFirstFace` / `getStartPosition` — initial valid 3-atom face
#   * `extendComponent` — queue-driven exploration via rolling
#   * `treatFace` + `treatEdge` — per-edge extension
#   * `thirdAtom` — rotate probe around an edge axis to find next face
#   * `correct(atom)` — shrink and reset when "PROBE TOUCHES FOUR ATOMS"
#     detected via zero rotation angle or multiple atoms tied at minimum
# =========================================================================

# Exception type for BALL's correct(atom) restart signal.
struct CorrectAtomException <: Exception
    atom::Int
end

# Computes the two probe-sphere centers tangent to atoms (a, b, c).
# Returns nothing if no such probe exists.
function _rs_probe_positions(atoms::Vector{Sphere{T}}, pr::T,
                              a::Int, b::Int, c::Int) where T
    s1 = Sphere{T}(atoms[a].center, atoms[a].r + pr)
    s2 = Sphere{T}(atoms[b].center, atoms[b].r + pr)
    s3 = Sphere{T}(atoms[c].center, atoms[c].r + pr)
    intersect_three_spheres(s1, s2, s3)
end

# Check if a probe at `probe_p` is unoccluded — i.e., no other atom
# protrudes into the probe sphere. Returns the index of the FIRST
# protruding atom (a "fourth atom" in BALL terminology), or 0 if none.
# The strict-inside check uses tolerance `eps`.
function _rs_first_protruding_atom(atoms::Vector{Sphere{T}}, pr::T,
                                    probe_p::Vector3{T},
                                    a::Int, b::Int, c::Int,
                                    neighbors::Vector{Vector{Int}},
                                    eps::T) where T
    nab = _common_neighbors(neighbors, a, b)
    isempty(nab) && return 0
    nac = _common_neighbors(neighbors, a, c)
    isempty(nac) && return 0
    i = 1; j = 1
    while i <= length(nab) && j <= length(nac)
        x = nab[i]; y = nac[j]
        if x == y
            if x != a && x != b && x != c
                rj = atoms[x].r + pr
                d = probe_p - Vector3{T}(atoms[x].center...)
                if dot(d, d) < rj * rj - eps
                    return x
                end
            end
            i += 1; j += 1
        elseif x < y
            i += 1
        else
            j += 1
        end
    end
    return 0
end

# Top-level rolling-probe wrapper: handles correct(atom) restart loop.
# WIP: not yet active. The enumeration-based path remains the entry point.
function _compute_reduced_surface_rolling(spheres::Vector{Sphere{T}},
                                            pr::T) where T
    # BALL's correct(atom) (reducedSurface.C:982): SHRINKS the atom by
    # `10 * Constants::EPSILON` (= 1e-3 Å during RS construction) and
    # re-runs. This is a small tie-breaking nudge.
    #
    # Additionally, BALL's STATUS_INSIDE marking (thirdAtom, lines
    # 1144-1158) buries atoms that lose the rotation-angle race during
    # rolling. My enumeration-based RS doesn't have rolling, so we
    # approximate STATUS_INSIDE via the SHRINKING repeated until any
    # tight 4-atom-degenerate cluster gets perturbed enough to break.
    # Repeated shrinks accumulate up to ~0.1 Å, effectively burying the
    # most-degenerate atoms — matching BALL's behavior on BPTI.
    shrunk = Dict{Int, T}()
    max_restarts = 128
    for _ in 1:max_restarts
        atoms = _apply_radius_corrections(spheres, shrunk)
        try
            rs = _compute_rs_rolling_inner(atoms, pr)
            # Post-construction: if there's a tight 4-atom-degenerate
            # cluster (probes within ~0.1 Å of each other on a 3-5-face
            # atom — BALL's correct-via-STATUS_INSIDE catches these via
            # rolling), shrink the involved atom by 1e-3 and retry.
            # Repeated shrinks accumulate until the degeneracy breaks.
            clustered_atom = _find_first_clustered_atom(rs)
            if clustered_atom != 0
                shrunk[clustered_atom] = get(shrunk, clustered_atom, zero(T)) + T(1e-3)
                # Cap per-atom shrinkage at 0.1 Å — beyond that, the
                # geometry is significantly altered and we should stop.
                shrunk[clustered_atom] > T(0.1) && return rs
                continue
            end
            return rs
        catch e
            e isa CorrectAtomException || rethrow()
            shrunk[e.atom] = get(shrunk, e.atom, zero(T)) + T(1e-3)
        end
    end
    atoms = _apply_radius_corrections(spheres, shrunk)
    return _compute_reduced_surface_one_pass(atoms, pr)
end

# Find the FIRST atom whose RSFaces are all tightly clustered (= 4-atom
# degeneracy signature). Returns 0 if none. Called repeatedly to drive
# shrink-and-retry.
function _find_first_clustered_atom(rs::ReducedSurface{T}) where T
    atom_faces = Dict{Int, Vector{Int}}()
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        for vi in (f.v1, f.v2, f.v3)
            a = rs.vertices[vi].atom
            push!(get!(atom_faces, a, Int[]), fi)
        end
    end
    cluster_tol2 = T(0.005)   # ~0.07 Å radius
    for (a, faces) in atom_faces
        length(faces) <= 5 || continue
        length(faces) >= 2 || continue
        c = Vector3{T}(0, 0, 0)
        for fi in faces; c = c + rs.faces[fi].center; end
        c = c / T(length(faces))
        cluster = true
        for fi in faces
            d = rs.faces[fi].center - c
            if dot(d, d) > cluster_tol2
                cluster = false
                break
            end
        end
        cluster && return a
    end
    return 0
end

# Detect atoms whose RSFaces ALL have tightly-clustered probe centers (the
# 4-atom-degeneracy signature). BALL's `correct(atom)` buries these
# explicitly via in-rolling angle detection; my enumeration path misses
# them when the angular separation is above the `_rs_third_atom`
# tolerance but the probes still coincide in 3-space within ~0.1 Å.
function _find_clustered_atoms(rs::ReducedSurface{T}) where T
    # Group RSFaces by each touching atom.
    atom_faces = Dict{Int, Vector{Int}}()
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        for vi in (f.v1, f.v2, f.v3)
            a = rs.vertices[vi].atom
            push!(get!(atom_faces, a, Int[]), fi)
        end
    end
    result = Set{Int}()
    cluster_tol2 = T(0.005)    # ~0.07 Å radius from centroid
    for (a, faces) in atom_faces
        length(faces) <= 5 || continue   # tight isolated 4-atom-degenerate cluster only
        length(faces) >= 2 || continue
        c = Vector3{T}(0, 0, 0)
        for fi in faces; c = c + rs.faces[fi].center; end
        c = c / T(length(faces))
        cluster = true
        for fi in faces
            d = rs.faces[fi].center - c
            if dot(d, d) > cluster_tol2
                cluster = false
                break
            end
        end
        cluster && push!(result, a)
    end
    result
end

# Apply burial: zero the radius of atoms in `buried` so they get filtered
# out of subsequent RS construction.
function _bury_atoms(spheres::Vector{Sphere{T}}, buried::Set{Int}) where T
    isempty(buried) && return spheres
    out = Vector{Sphere{T}}(undef, length(spheres))
    for (i, s) in enumerate(spheres)
        out[i] = i in buried ? Sphere{T}(s.center, zero(T)) : s
    end
    out
end

# Find third atom by rolling probe around axis (atom_a, atom_b) starting
# from `start_probe`. Port of BALL's `thirdAtom` (reducedSurface.C:1081).
# Returns (atom_c, probe_p, rotation_angle) or throws CorrectAtomException
# on degeneracy.
function _rs_third_atom(atoms::Vector{Sphere{T}}, pr::T,
                          neighbors::Vector{Vector{Int}},
                          atom_a::Int, atom_b::Int,
                          start_probe::Vector3{T},
                          third_face_atom::Int,
                          face_normal::Vector3{T}) where T
    # BALL's logic: enumerate atoms touched by probe at the (a, b) rolling
    # axis. For each candidate, compute rotation angle from start_probe.
    # Pick the candidate with the smallest positive angle. Detect:
    #   * Rotation angle == 0 (touches 4 atoms simultaneously) → correct.
    #   * Multiple atoms tie at same minimum angle → correct each but first.
    # Set up rotation axis and reference circle.
    pa = Vector3{T}(atoms[atom_a].center...)
    pb = Vector3{T}(atoms[atom_b].center...)
    axis = pb - pa
    third_p = Vector3{T}(atoms[third_face_atom].center...)
    test_vec = cross(face_normal, axis)
    if dot(test_vec, third_p) < dot(test_vec, pa)
        axis = -axis
    end
    # Probe rolling: probe traces a circle perpendicular to (a, b) axis.
    # Intersection of (atom_a + pr) and (atom_b + pr) spheres = a circle.
    # Reference vector from this circle's center to start_probe.
    s1 = Sphere{T}(atoms[atom_a].center, atoms[atom_a].r + pr)
    s2 = Sphere{T}(atoms[atom_b].center, atoms[atom_b].r + pr)
    # Find intersection circle of s1 and s2 (a circle perpendicular to
    # axis passing through valid probe positions).
    pts = intersect_three_spheres(
        Sphere{T}(atoms[atom_a].center, atoms[atom_a].r + pr),
        Sphere{T}(atoms[atom_b].center, atoms[atom_b].r + pr),
        Sphere{T}(start_probe, T(0)),  # degenerate — start_probe IS on both spheres
    )
    # Simpler: circle center = projection of probe onto axis line.
    # axis is from atom_a to atom_b. Project start_probe onto axis line.
    ax_len = sqrt(dot(axis, axis))
    ax_hat = ax_len > eps(T) ? axis / ax_len : Vector3{T}(0, 0, 1)
    t = dot(start_probe - pa, ax_hat)
    circle_p = pa + t * ax_hat
    v1 = start_probe - circle_p
    # Enumerate candidate third atoms.
    common = _common_neighbors(neighbors, atom_a, atom_b)
    min_angle = T(3 * π)
    candidates_at_min = Tuple{Int, Vector3{T}, T}[]
    for c in common
        c == atom_a && continue
        c == atom_b && continue
        c == third_face_atom && continue
        # Probe positions touching (a, b, c)
        probe_pts = _rs_probe_positions(atoms, pr, atom_a, atom_b, c)
        probe_pts === nothing && continue
        for probe in probe_pts
            # Rotation angle from start_probe to this probe around axis.
            v2 = probe - circle_p
            ang = oriented_angle(v1, v2, axis)
            # BALL: ignore start position (angle 0 from same probe).
            if (probe[1] ≈ start_probe[1] && probe[2] ≈ start_probe[2] &&
                probe[3] ≈ start_probe[3])
                continue
            end
            # BALL's degeneracy: angle ≈ 0 (probe touches 4 atoms). BALL
            # uses `Maths::isZero` with Constants::EPSILON = 1e-4 during
            # RS construction (reducedSurface.C:677).
            if ang < T(1e-4) || ang > 2 * T(π) - T(1e-4)
                throw(CorrectAtomException(c))
            end
            if ang < min_angle - T(1e-4)
                empty!(candidates_at_min)
                min_angle = ang
                push!(candidates_at_min, (c, probe, ang))
            elseif abs(ang - min_angle) < T(1e-4)
                push!(candidates_at_min, (c, probe, ang))
            end
        end
    end
    isempty(candidates_at_min) && return nothing
    # BALL: if multiple candidates tie at min angle, fire correct on all
    # but the first.
    if length(candidates_at_min) > 1
        for k in 2:length(candidates_at_min)
            throw(CorrectAtomException(candidates_at_min[k][1]))
        end
    end
    c, probe, ang = candidates_at_min[1]
    return (c, probe, ang)
end

# Inner rolling-probe construction. Throws CorrectAtomException to
# signal a needed atom shrink + restart (BALL's exception pattern).
function _compute_rs_rolling_inner(atoms::Vector{Sphere{T}}, pr::T) where T
    rs = _compute_reduced_surface_one_pass(atoms, pr)
    isempty(rs.faces) && return rs
    neighbors, _, _ = _build_neighbors(atoms, pr)
    for f in rs.faces
        f.v1 == 0 && continue
        atoms_v = (rs.vertices[f.v1].atom,
                   rs.vertices[f.v2].atom,
                   rs.vertices[f.v3].atom)
        for k in 1:3
            a = atoms_v[k]; b = atoms_v[mod1(k + 1, 3)]
            third_atom = atoms_v[mod1(k + 2, 3)]
            _rs_third_atom(atoms, pr, neighbors, a, b, f.center,
                            third_atom, f.normal)
        end
    end
    return rs
end

# Simulate BALL's STATUS_INSIDE marking (thirdAtom lines 1144-1158) via
# BFS over the enumerated RS faces. For each rolling step (along an
# atom-pair edge), the candidate with smallest positive rotation angle
# wins; all losers get their unique "third atom" marked STATUS_INSIDE
# and never considered as a winner in subsequent rollings. After BFS,
# faces involving any STATUS_INSIDE atom are removed.
#
# This catches the buried-atom configurations BALL's rolling naturally
# eliminates but my enumeration includes (e.g. atoms 13, 304 on BPTI).
function _ball_simulate_status_inside!(rs::ReducedSurface{T}) where T
    isempty(rs.faces) && return rs
    n_atoms = length(rs.atoms)

    # Build pair → faces map for fast candidate lookup.
    pair_faces = Dict{NTuple{2,Int}, Vector{Int}}()
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        a = rs.vertices[f.v1].atom
        b = rs.vertices[f.v2].atom
        c = rs.vertices[f.v3].atom
        for (x, y) in ((a, b), (a, c), (b, c))
            k = x < y ? (x, y) : (y, x)
            push!(get!(() -> Int[], pair_faces, k), fi)
        end
    end

    # 0=UNKNOWN, 1=ON_SURFACE, 2=INSIDE. Indexed by sphere index.
    status = zeros(Int, n_atoms)

    # Seed: BALL picks an extreme atom (most negative x); we pick the
    # face whose centroid has the smallest x-coordinate (a proxy for
    # "extreme"). This deterministically picks a true surface face.
    seed = 0
    best_x = typemax(T)
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        xc = f.center[1]
        if xc < best_x
            best_x = xc
            seed = fi
        end
    end
    seed == 0 && return rs
    sf = rs.faces[seed]
    for vi in (sf.v1, sf.v2, sf.v3)
        status[rs.vertices[vi].atom] = 1
    end

    visited = Set{Int}(seed)
    # Process each (a, b) atom-pair rolling exactly ONCE (mirrors BALL's
    # treatEdge which is called per RSEdge, not per face-side incidence).
    # Without this, a single pair gets thirdAtom-processed multiple times
    # (once per face referencing it), spuriously marking same-angle
    # candidates INSIDE in later passes.
    processed_pairs = Set{NTuple{2,Int}}()
    queue = Int[seed]
    while !isempty(queue)
        fi = popfirst!(queue)
        f = rs.faces[fi]
        f.v1 == 0 && continue
        ar = rs.vertices[f.v1].atom
        br = rs.vertices[f.v2].atom
        cr = rs.vertices[f.v3].atom
        for (xa, ya, ta) in ((ar, br, cr), (ar, cr, br), (br, cr, ar))
            key = xa < ya ? (xa, ya) : (ya, xa)
            key in processed_pairs && continue
            push!(processed_pairs, key)
            cand_faces = get(pair_faces, key, Int[])
            length(cand_faces) <= 1 && continue
            # Per-face rolling axis (matches `_build_atom_lobes`).
            pa = Vector3{T}(rs.atoms[xa].center...)
            pb = Vector3{T}(rs.atoms[ya].center...)
            ax = pa - pb
            test_vec = cross(f.normal, ax)
            third_p = Vector3{T}(rs.atoms[ta].center...)
            if dot(test_vec, third_p) > dot(test_vec, pa)
                ax = -ax
            end
            ax_len2 = dot(ax, ax)
            rel_f = f.center - pa
            proj_f = rel_f - (dot(rel_f, ax) / ax_len2) * ax
            # Find smallest-angle candidate (not currently INSIDE).
            best_ang = T(3 * π)
            best_face = 0
            for cf in cand_faces
                cf == fi && continue
                g = rs.faces[cf]
                g.v1 == 0 && continue
                # g's third atom.
                gt = 0
                for vi in (g.v1, g.v2, g.v3)
                    a_ = rs.vertices[vi].atom
                    a_ != xa && a_ != ya && (gt = a_)
                end
                gt == 0 && continue
                status[gt] == 2 && continue  # already INSIDE
                rel_g = g.center - pa
                proj_g = rel_g - (dot(rel_g, ax) / ax_len2) * ax
                ang = oriented_angle(proj_f, proj_g, ax)
                if ang > T(1e-6) && ang < best_ang
                    best_ang = ang
                    best_face = cf
                end
            end
            best_face == 0 && continue
            # Mark winner's third ON_SURFACE; queue it.
            wg = rs.faces[best_face]
            for vi in (wg.v1, wg.v2, wg.v3)
                a_ = rs.vertices[vi].atom
                a_ != xa && a_ != ya && status[a_] == 0 && (status[a_] = 1)
            end
            if best_face ∉ visited
                push!(visited, best_face)
                push!(queue, best_face)
            end
            # Mark all OTHER candidates' thirds as INSIDE (= "lost the
            # rolling race" in this step). BALL's thirdAtom lines 1144,
            # 1155: only mark if status was UNKNOWN.
            for cf in cand_faces
                (cf == fi || cf == best_face) && continue
                g = rs.faces[cf]
                g.v1 == 0 && continue
                for vi in (g.v1, g.v2, g.v3)
                    a_ = rs.vertices[vi].atom
                    a_ != xa && a_ != ya && status[a_] == 0 && (status[a_] = 2)
                end
            end
        end
    end

    # Remove faces involving any INSIDE atom.
    for f in rs.faces
        f.v1 == 0 && continue
        for vi in (f.v1, f.v2, f.v3)
            a_ = rs.vertices[vi].atom
            if status[a_] == 2
                f.v1 = 0
                break
            end
        end
    end
    rs
end

# Core RS construction (formerly the body of compute_reduced_surface).
# Keeps the input sphere indexing intact: zero-radius atoms remain in
# rs.atoms but are skipped when enumerating triples (so external code can
# refer to them by their original sphere index, e.g. for burial decisions
# made across restarts).
function _compute_reduced_surface_one_pass(spheres::Vector{Sphere{T}}, pr::T) where T
    atoms = spheres  # keep all indices stable; treat r <= 0.001 as buried
    if all(s -> s.r <= T(0.001), atoms)
        return ReducedSurface{T}(
            atoms, pr,
            RSVertex[], RSEdge{T}[], RSFace{T}[],
            BoundingBox{T}(zero(Vector3{T}), zero(Vector3{T})),
            zero(T),
        )
    end

    neighbors, bb, rmax = _build_neighbors(atoms, pr)
    rs = ReducedSurface{T}(
        atoms, pr,
        RSVertex[], RSEdge{T}[], RSFace{T}[],
        bb, rmax,
    )

    # Enumerate all valid faces. A face is a tuple (a,b,c, slot, probe).
    # `slot ∈ (1,2)` indexes the two sphere–sphere–sphere solutions.
    eps = T(1e-8)
    face_keys = NTuple{4,Int}[]            # (a, b, c, slot), sorted a<b<c
    face_probes = Vector3{T}[]
    seen_triples = Set{NTuple{4,Int}}()

    for a in 1:length(atoms)
        for b in neighbors[a]
            b <= a && continue
            common = _common_neighbors(neighbors, a, b)
            for c in common
                c <= b && continue
                key3 = (a, b, c)
                # three-sphere intersection (radii inflated by probe)
                s1 = Sphere{T}(atoms[a].center, atoms[a].r + pr)
                s2 = Sphere{T}(atoms[b].center, atoms[b].r + pr)
                s3 = Sphere{T}(atoms[c].center, atoms[c].r + pr)
                pts = intersect_three_spheres(s1, s2, s3)
                pts === nothing && continue
                for slot in 1:2
                    fk = (a, b, c, slot)
                    fk in seen_triples && continue
                    push!(seen_triples, fk)
                    probe = slot == 1 ? pts[1] : pts[2]
                    _probe_unoccluded(probe, atoms, pr, neighbors, a, b, c, eps) || continue
                    push!(face_keys, fk)
                    push!(face_probes, probe)
                end
            end
        end
    end

    # Pre-compute face normals; needed by lobe assignment (rolling-forward
    # direction depends on per-face normal).
    face_normals = Vector{Vector3{T}}(undef, length(face_keys))
    for (i, fk) in enumerate(face_keys)
        a, b, c, _ = fk
        face_normals[i] = _face_normal(rs, a, b, c, face_probes[i])
    end

    # Build RSVertices: one per *lobe* on each atom. Two faces at atom `a`
    # belong to the same lobe iff connected through rolling-probe edges that
    # share atom `a`. Mirrors BALL's `vertices_[atom]` where one atom may
    # carry multiple disconnected RSVertices.
    atom_lobes = _build_atom_lobes(face_keys, face_probes, face_normals, atoms)
    atom_lobe_to_vertex = Dict{NTuple{2,Int}, Int}()
    for a in 1:length(atoms)
        for (_, lobe_id) in atom_lobes[a]
            key = (a, lobe_id)
            if !haskey(atom_lobe_to_vertex, key)
                push!(rs.vertices, RSVertex(a))
                atom_lobe_to_vertex[key] = length(rs.vertices)
            end
        end
    end

    for (i, fk) in enumerate(face_keys)
        a, b, c, _ = fk
        v1 = atom_lobe_to_vertex[(a, atom_lobes[a][i])]
        v2 = atom_lobe_to_vertex[(b, atom_lobes[b][i])]
        v3 = atom_lobe_to_vertex[(c, atom_lobes[c][i])]
        f = RSFace{T}(v1, v2, v3, Int[], face_probes[i],
                      face_normals[i], false)
        # BALL `RSComputer::correctProbePosition` (reducedSurface.C:832):
        # `new_face->singular_ = Maths::isLess(GetDistance(probe.p, plane),
        # rs_->probe_radius_)`. `GetDistance` returns the SIGNED PERPENDICULAR
        # DISTANCE (= dot product with the UNIT normal). My port's
        # plane_distance is just dot(p - a, n) — correct only if `n` is unit-
        # normalized. `f.normal` is NOT normalized at this stage, so the
        # comparison was always against the scaled value. Normalize the
        # normal before the distance check.
        nrm = sqrt(dot(f.normal, f.normal))
        f.singular = nrm > eps(T) && abs(plane_distance(
            face_probes[i],
            Vector3{T}(atoms[a].center...),
            f.normal,
        )) / nrm < pr - T(1e-6)
        push!(rs.faces, f)
        push!(rs.vertices[v1].faces, length(rs.faces))
        push!(rs.vertices[v2].faces, length(rs.faces))
        push!(rs.vertices[v3].faces, length(rs.faces))
    end

    _build_edges!(rs, neighbors)
    _add_free_edges_and_isolated_vertices!(rs, neighbors)
    _detect_singular_edges!(rs)
    # Note: `_delete_similar_faces!` is implemented (see below) but
    # is NOT called here — BALL only invokes deleteSimilarFaces for a
    # specific SES degeneracy (9-edge SESFaces in
    # solventExcludedSurface.C:1687), not as a routine RS post-pass.
    # For normal inputs, same-atom-set RSFace pairs are a LEGITIMATE
    # part of the topology (probe rolling above and below an atom
    # triple plane), and merging them would destroy the closed surface.
    rs
end

# Port of BALL's `ReducedSurface::deleteSimilarFaces`
# (reducedSurface.C:272). Two RSFaces are "similar" when they have the
# same 3-atom set (in any permutation). The cleanup merges them: similar
# vertices are joined (edge/face incidence lists union), similar edges
# are corrected (either deleted when both faces share them or rerouted),
# and the two faces themselves are removed.
#
# Operates in-place on `rs`. Deleted faces/edges/vertices are marked
# with sentinel field values: face.v1 = 0, edge.v1 = 0, vertex.atom = 0.
# Downstream code (compute_ses, triangulator) must skip these entries.
function _delete_similar_faces!(rs::ReducedSurface{T}) where T
    n_faces = length(rs.faces)
    # Build a map: sorted atom triple → list of face indices with that set.
    by_atoms = Dict{NTuple{3, Int}, Vector{Int}}()
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        atoms = (rs.vertices[f.v1].atom,
                 rs.vertices[f.v2].atom,
                 rs.vertices[f.v3].atom)
        key = NTuple{3, Int}(sort(collect(atoms)))
        push!(get!(by_atoms, key, Int[]), fi)
    end
    n_merged = 0
    for (_, face_list) in by_atoms
        length(face_list) >= 2 || continue
        # Process pairs (we expect exactly 2 in non-degenerate similar cases).
        i = 1
        while i + 1 <= length(face_list)
            f1_idx = face_list[i]
            f2_idx = face_list[i + 1]
            if rs.faces[f1_idx].v1 != 0 && rs.faces[f2_idx].v1 != 0
                _merge_similar_rs_faces!(rs, f1_idx, f2_idx)
                n_merged += 1
            end
            i += 2
        end
    end
    n_merged
end

# Helper: merge two similar RSFaces. Mirrors BALL's joinVertices +
# correctEdges sequence (reducedSurface.C:351-418).
function _merge_similar_rs_faces!(rs::ReducedSurface{T},
                                   f1_idx::Int, f2_idx::Int) where T
    f1 = rs.faces[f1_idx]
    f2 = rs.faces[f2_idx]
    # findSimilarVertices: for each of face1's vertices, find face2's
    # vertex with the same atom (= same vertex if both faces share the
    # same RSVertex objects; otherwise a different RSVertex for that
    # atom).
    v1_arr = (f1.v1, f1.v2, f1.v3)
    v2_arr = (f2.v1, f2.v2, f2.v3)
    pair_v = ntuple(3) do k
        atom_target = rs.vertices[v1_arr[k]].atom
        for j in 1:3
            if rs.vertices[v2_arr[j]].atom == atom_target
                return (v1_arr[k], v2_arr[j])
            end
        end
        return (v1_arr[k], 0)
    end
    # findSimilarEdges: for each of face1's edges, find face2's edge
    # with the same atom pair. Both faces may have more than 3 edges
    # under the stacked-probes-on-one-axis case; pair every entry of f1
    # to its same-atom-pair counterpart in f2 (or 0 if absent).
    function edge_atoms(eid)
        eid == 0 && return (0, 0)
        re = rs.edges[eid]
        a = rs.vertices[re.v1].atom
        b = rs.vertices[re.v2].atom
        a < b ? (a, b) : (b, a)
    end
    pair_e = Tuple{Int,Int}[]
    for e1_id in f1.edges
        eatoms_1 = edge_atoms(e1_id)
        match = 0
        for e2_id in f2.edges
            if edge_atoms(e2_id) == eatoms_1
                match = e2_id; break
            end
        end
        push!(pair_e, (e1_id, match))
    end
    # joinVertices: merge edge/face lists of similar vertex pairs.
    for (v1_id, v2_id) in pair_v
        v2_id == 0 && continue
        if v1_id != v2_id
            v1 = rs.vertices[v1_id]
            v2 = rs.vertices[v2_id]
            union!(v1.edges, v2.edges)
            union!(v1.faces, v2.faces)
            # Substitute v2 → v1 in all edges referencing it.
            for eidx in v2.edges
                re = rs.edges[eidx]
                if re.v1 == v2_id
                    re.v1 = v1_id
                elseif re.v2 == v2_id
                    re.v2 = v1_id
                end
            end
            # Substitute v2 → v1 in all faces referencing it.
            for fidx in v2.faces
                rf = rs.faces[fidx]
                if rf.v1 == v2_id
                    rf.v1 = v1_id
                elseif rf.v2 == v2_id
                    rf.v2 = v1_id
                elseif rf.v3 == v2_id
                    rf.v3 = v1_id
                end
            end
            rs.vertices[v2_id].atom = 0  # mark deleted
            empty!(v2.edges)
            empty!(v2.faces)
        end
        # Erase face1 and face2 from v1's face list.
        delete!(rs.vertices[v1_id].faces, f1_idx)
        delete!(rs.vertices[v1_id].faces, f2_idx)
    end
    # correctEdges
    for (e1_id, e2_id) in pair_e
        e2_id == 0 && continue
        if e1_id == e2_id
            # Same edge shared by both faces.
            re = rs.edges[e1_id]
            if re.singular
                # Delete the singular edge.
                delete!(rs.vertices[re.v1].edges, e1_id)
                delete!(rs.vertices[re.v2].edges, e1_id)
                re.v1 = 0
                re.v2 = 0
            else
                # Free the edge: it'll become a free-rolling toric.
                re.f1 = 0
                re.f2 = 0
                re.angle = T(2 * π)
            end
        else
            # Different edges with the same atom pair: re-route e1
            # through e2's neighbour, delete e2.
            re1 = rs.edges[e1_id]
            re2 = rs.edges[e2_id]
            neighbour2 = re2.f1 == f2_idx ? re2.f2 : re2.f1
            if re1.f1 == f1_idx
                re1.f1 = neighbour2
            elseif re1.f2 == f1_idx
                re1.f2 = neighbour2
            end
            if neighbour2 != 0
                rf_n = rs.faces[neighbour2]
                idx = findfirst(==(e2_id), rf_n.edges)
                idx !== nothing && (rf_n.edges[idx] = e1_id)
            end
            delete!(rs.vertices[re2.v1].edges, e2_id)
            delete!(rs.vertices[re2.v2].edges, e2_id)
            re2.v1 = 0
            re2.v2 = 0
        end
    end
    # Mark the faces as deleted.
    f1.v1 = 0
    f2.v1 = 0
end

# Mark each RS edge as singular when the probe sphere on its torus
# self-intersects (i.e. the torus is a spindle: probe_radius > radius_of_torus
# and the probe at any rolling position crosses the (a,b) line). Compute the
# two singular vertices — the line/probe-sphere intersection points — and
# store them in `intersection1`, `intersection2` (ordered with the one closer
# to atom of v1 first).
function _detect_singular_edges!(rs::ReducedSurface{T}) where T
    pr = rs.probe_radius
    for e in rs.edges
        e.singular = false
        e.radius_of_torus <= 0 && continue
        e.radius_of_torus < pr - T(1e-9) || continue
        # Pick a probe center on the rolling circle; if `f1` is unset use the
        # torus center perturbed along an arbitrary radial direction.
        probe_center = e.f1 != 0 ? rs.faces[e.f1].center :
                                   e.center_of_torus  # degenerate placeholder
        atom_a = rs.atoms[rs.vertices[e.v1].atom].center
        atom_b = rs.atoms[rs.vertices[e.v2].atom].center
        line_dir = Vector3{T}((atom_b - atom_a)...)
        pts = intersect_sphere_line(
            Sphere{T}(probe_center, pr),
            Vector3{T}(atom_a...),
            line_dir,
        )
        pts === nothing && continue
        ip1, ip2 = pts
        # Order so ip1 is closer to atom_a.
        d1 = ip1 - Vector3{T}(atom_a...)
        d2 = ip2 - Vector3{T}(atom_a...)
        if dot(d2, d2) < dot(d1, d1)
            ip1, ip2 = ip2, ip1
        end
        e.intersection1 = ip1
        e.intersection2 = ip2
        e.singular = true
    end
end

# After face/edge construction, any atom whose probe surface is *not* fully
# obstructed by face-generating neighbours can still contribute to the RS as
# an isolated vertex (1 atom) or via a free edge (2-atom torus). This matches
# BALL's `findFirstVertex`/`createFreeEdge` fall-throughs.
function _add_free_edges_and_isolated_vertices!(rs::ReducedSurface{T},
                                                neighbors::Vector{Vector{Int}}) where T
    pr = rs.probe_radius
    n = length(rs.atoms)
    has_vertex = falses(n)
    for v in rs.vertices
        has_vertex[v.atom] = true
    end

    # Free edges: every neighbour pair (a, b) that has no faces around it
    # contributes one free edge with full 2π rolling — provided no third
    # atom can host a probe on this pair (else it would already be a face).
    seen_pairs = Set{NTuple{2,Int}}()
    for a in 1:n, b in neighbors[a]
        b <= a && continue
        (a, b) in seen_pairs && continue
        push!(seen_pairs, (a, b))
        # Does either incident vertex already share a face involving (a, b)?
        # If so, _build_edges! has handled this pair already.
        if has_vertex[a] && has_vertex[b]
            va = 0; vb = 0
            for (vi, v) in enumerate(rs.vertices)
                v.atom == a && (va = vi)
                v.atom == b && (vb = vi)
            end
            already_connected = false
            if va != 0 && vb != 0
                for ei in rs.vertices[va].edges
                    e = rs.edges[ei]
                    if (e.v1 == va && e.v2 == vb) || (e.v1 == vb && e.v2 == va)
                        already_connected = true
                        break
                    end
                end
            end
            already_connected && continue
        end
        # Does some third atom support a probe on (a, b)?  If yes, that's a
        # face-producing pair and we leave it to _build_edges!.
        common = _common_neighbors(neighbors, a, b)
        face_capable = false
        for c in common
            s1 = Sphere{T}(rs.atoms[a].center, rs.atoms[a].r + pr)
            s2 = Sphere{T}(rs.atoms[b].center, rs.atoms[b].r + pr)
            s3 = Sphere{T}(rs.atoms[c].center, rs.atoms[c].r + pr)
            if intersect_three_spheres(s1, s2, s3) !== nothing
                face_capable = true
                break
            end
        end
        face_capable && continue

        # Free edge with full 2π rolling. Requires contact circles to exist.
        ccs = _contact_circles(rs, a, b)
        ccs === nothing && continue
        c1, c2, c3 = ccs
        c1.r <= pr && continue   # torus too small

        # Materialise vertices on both atoms if they don't already exist.
        function ensure_vertex(atom_idx)
            for (vi, v) in enumerate(rs.vertices)
                v.atom == atom_idx && return vi
            end
            push!(rs.vertices, RSVertex(atom_idx))
            length(rs.vertices)
        end
        va = ensure_vertex(a)
        vb = ensure_vertex(b)
        e = RSEdge{T}(va, vb)
        e.center_of_torus = c1.p
        e.radius_of_torus = c1.r
        e.angle           = 2 * T(π)
        e.contact_circle1 = c2
        e.contact_circle2 = c3
        push!(rs.edges, e)
        eidx = length(rs.edges)
        push!(rs.vertices[va].edges, eidx)
        push!(rs.vertices[vb].edges, eidx)
        has_vertex[a] = true
        has_vertex[b] = true
    end

    # Isolated vertices: atoms with positive radius and no neighbours and no
    # existing RS vertex on them.
    for a in 1:n
        has_vertex[a] && continue
        isempty(neighbors[a]) || continue
        push!(rs.vertices, RSVertex(a))
    end
end

# Group faces by sorted atom pair → list of (face_idx, the third atom).
function _build_edges!(rs::ReducedSurface{T},
                       neighbors::Vector{Vector{Int}}) where T
    # Group RSFaces by atom pair (= rolling axis between 2 atoms). For
    # each atom pair, the RSFaces that share it are connected in a cycle
    # (or chain) around the rolling axis.
    pair_to_faces = Dict{NTuple{2,Int}, Vector{Int}}()
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        a = rs.vertices[f.v1].atom
        b = rs.vertices[f.v2].atom
        c = rs.vertices[f.v3].atom
        for (x, y) in ((a, b), (a, c), (b, c))
            key = _sort2(x, y)
            push!(get!(() -> Int[], pair_to_faces, key), fi)
        end
    end

    pr = rs.probe_radius

    for (key, faces) in pair_to_faces
        a, b = key
        n = length(faces)
        if n == 1
            # Free edge: probe can roll freely; no second face.
            f = rs.faces[faces[1]]
            v1, v2 = _face_vertex_indices_for(rs, f, a, b)
            e = RSEdge{T}(v1, v2)
            e.f1 = faces[1]
            ccs = _contact_circles(rs, a, b)
            if ccs !== nothing
                c1, c2, c3 = ccs
                e.center_of_torus = c1.p
                e.radius_of_torus = c1.r
                e.contact_circle1 = c2
                e.contact_circle2 = c3
                e.angle = 2 * T(π)
            end
            _attach_edge!(rs, e, faces[1])
            continue
        end

        # n >= 2: connect each face to its rolling-forward neighbour.
        # Mirrors BALL's thirdAtom (reducedSurface.C:1099-1106): the
        # rolling axis is signed per face so that "forward" means rolling
        # AWAY from the face's third atom. The next face along (a, b) is
        # the one with the smallest positive rotation angle in this
        # direction. When the n faces split into multiple angular clusters
        # whose forward-direction rollings close on themselves, the cycle
        # decomposes into separate sub-cycles -> multiple RSVertex lobes.
        pa = Vector3{T}(rs.atoms[a].center...)
        pb = Vector3{T}(rs.atoms[b].center...)
        axis_ab = pa - pb     # BALL: atom[atom1].p - atom[atom2].p
        ccs = _contact_circles(rs, a, b)
        torus_p = ccs === nothing ? Vector3{T}(0,0,0) : ccs[1].p
        torus_r = ccs === nothing ? zero(T) : ccs[1].r
        cc1 = ccs === nothing ? Circle3{T}(Vector3{T}(0,0,0), Vector3{T}(0,0,0), zero(T)) : ccs[2]
        cc2 = ccs === nothing ? Circle3{T}(Vector3{T}(0,0,0), Vector3{T}(0,0,0), zero(T)) : ccs[3]

        # For each face f_k, compute (signed) axis and the next face index
        # (within the `faces` list) along rolling-forward.
        next_face = zeros(Int, n)
        delta_angle = zeros(T, n)
        # Pre-fetch third atoms and face references.
        thirds = Vector{Int}(undef, n)
        for k in 1:n
            f = rs.faces[faces[k]]
            t = 0
            for vi in (f.v1, f.v2, f.v3)
                x = rs.vertices[vi].atom
                x != a && x != b && (t = x)
            end
            thirds[k] = t
        end
        for k in 1:n
            f_k = rs.faces[faces[k]]
            ax = axis_ab
            test_vec = cross(f_k.normal, ax)
            third_p = Vector3{T}(rs.atoms[thirds[k]].center...)
            if dot(test_vec, third_p) > dot(test_vec, pa)
                ax = -ax
            end
            # Rotation angle from f_k.center to each other face's center
            # around ax (signed, in [0, 2π)).
            rel_k = f_k.center - torus_p
            proj_k = rel_k - (dot(rel_k, ax) / dot(ax, ax)) * ax
            best_angle = T(3 * π)
            best_j = 0
            for j in 1:n
                j == k && continue
                f_j = rs.faces[faces[j]]
                rel_j = f_j.center - torus_p
                proj_j = rel_j - (dot(rel_j, ax) / dot(ax, ax)) * ax
                ang = oriented_angle(proj_k, proj_j, ax)
                if ang > T(1e-6) && ang < best_angle
                    best_angle = ang
                    best_j = j
                end
            end
            next_face[k] = best_j
            delta_angle[k] = best_angle
        end

        # Materialise edges. F→next(F) gives one directed edge per face;
        # dedup so that F→F' and F'→F produce a single undirected RSEdge.
        added = Set{Tuple{Int,Int}}()
        for k in 1:n
            nk = next_face[k]
            nk == 0 && continue
            fi1 = faces[k]
            fi2 = faces[nk]
            key_e = fi1 < fi2 ? (fi1, fi2) : (fi2, fi1)
            key_e in added && continue
            push!(added, key_e)
            f1 = rs.faces[fi1]
            v1, v2 = _face_vertex_indices_for(rs, f1, a, b)
            e = RSEdge{T}(v1, v2)
            e.f1 = fi1
            e.f2 = fi2
            e.center_of_torus = torus_p
            e.radius_of_torus = torus_r
            e.contact_circle1 = cc1
            e.contact_circle2 = cc2
            e.angle = delta_angle[k]
            _attach_edge!(rs, e, fi1, fi2)
        end
    end
end


# Find the two vertex indices of `f` that sit on atoms `a` and `b`.
function _face_vertex_indices_for(rs::ReducedSurface{T}, f::RSFace{T},
                                  a::Int, b::Int) where T
    v1 = v2 = 0
    for vi in (f.v1, f.v2, f.v3)
        atom = rs.vertices[vi].atom
        if atom == a; v1 = vi; end
        if atom == b; v2 = vi; end
    end
    (v1, v2)
end

function _attach_edge!(rs::ReducedSurface{T}, e::RSEdge{T},
                       faces_to_update::Vararg{Int,N}) where {T, N}
    push!(rs.edges, e)
    eidx = length(rs.edges)
    push!(rs.vertices[e.v1].edges, eidx)
    push!(rs.vertices[e.v2].edges, eidx)
    # Append to the face's edge list (BALL stores incident edges as a
    # list; we use Vector{Int}). No slot-overwrite is possible — each
    # call appends one new edge, so stacked-probe configurations
    # accumulate the corresponding multiple incidences naturally.
    for fi in faces_to_update
        f = rs.faces[fi]
        eidx in f.edges || push!(f.edges, eidx)
    end
end

# Build a Vector{Sphere{T}} from an atom container.
function _atom_spheres_for_rs(ac::AbstractAtomContainer{T}) where T
    at = atoms(ac)
    out = Sphere{T}[]
    for i in eachindex(at.r)
        r = at.radius[i]
        r > T(0.001) || continue
        push!(out, Sphere{T}(Point3{T}(at.r[i]...), r))
    end
    out
end

# ===========================================================================
# BALL-faithful rolling-probe RS construction.
#
# Port of BALL's `RSComputer::run()` (reducedSurface.C:670). Replaces the
# enumeration-based path with the actual rolling algorithm: seed face, then
# treatFace/treatEdge propagation with STATUS_INSIDE marking of losers in
# each thirdAtom step. This is what produces BALL's exact RS topology.
# ===========================================================================

# BALL `reducedSurface.h:502-507`:
#   enum AtomStatus { STATUS_ON_SURFACE=0, STATUS_INSIDE=1, STATUS_UNKNOWN=2 }
# Mirrored exactly so numeric comparisons against atom_status[] match BALL.
const _RS_STATUS_ON_SURFACE = 0
const _RS_STATUS_INSIDE     = 1
const _RS_STATUS_UNKNOWN    = 2

# Per-atom triple probe-position cache. Stores up to 2 probe positions
# (sphere-sphere-sphere intersection) per sorted atom triple. status[i]=0
# (untested), 1 (OK), 2 (NOT_OK = occluded).
mutable struct _ProbePos{T<:Real}
    valid::Bool   # true iff intersect_three_spheres returned points
    p0::Vector3{T}
    p1::Vector3{T}
    status0::Int
    status1::Int
end

mutable struct _RSComputer{T<:Real}
    rs::ReducedSurface{T}
    atoms::Vector{Sphere{T}}
    pr::T
    neighbors::Vector{Vector{Int}}
    atom_status::Vector{Int}
    # Per-atom list of RSVertex indices. atom_idx → vector of vertex indices.
    vertices_per_atom::Vector{Vector{Int}}
    # Cached probe positions keyed by sorted atom triple.
    probe_cache::Dict{NTuple{3,Int}, _ProbePos{T}}
    # Cached common neighbors for atom pairs.
    pair_neighbors::Dict{NTuple{2,Int}, Vector{Int}}
    # Face processing queue (treatFace).
    new_faces::Vector{Int}
    # Vertex queue for extendComponent.
    new_vertices::Vector{Int}
    # Removed vertex indices (for dangling-ptr detection).
    rm_vertices::Set{Int}
    # BALL `rs_->insert(edge)` counter — assigns `edge->index_` at the moment
    # `treatEdge` succeeds. We track this per-edge via RSEdge.ball_idx and
    # permute rs.edges at end of construction to match BALL's edge ordering.
    ball_edge_counter::Base.RefValue{Int}
end

function _RSComputer{T}(atoms::Vector{Sphere{T}}, pr::T) where T
    n = length(atoms)
    nb, bb, rmax = _build_neighbors(atoms, pr)
    rs = ReducedSurface{T}(atoms, pr, RSVertex[], RSEdge{T}[], RSFace{T}[], bb, rmax)
    _RSComputer{T}(rs, atoms, pr, nb, fill(_RS_STATUS_UNKNOWN, n),
                   [Int[] for _ in 1:n],
                   Dict{NTuple{3,Int}, _ProbePos{T}}(),
                   Dict{NTuple{2,Int}, Vector{Int}}(),
                   Int[], Int[], Set{Int}(), Ref(0))
end

# Cached common neighbors of (a, b). Returns ATOMS within reach of both.
function _cm_common_neighbors(cm::_RSComputer{T}, a::Int, b::Int) where T
    key = a < b ? (a, b) : (b, a)
    get!(cm.pair_neighbors, key) do
        _common_neighbors(cm.neighbors, a, b)
    end
end

# Cached probe positions for atom triple (a, b, c). Returns (valid, p0, p1).
function _cm_probe_positions(cm::_RSComputer{T}, a::Int, b::Int, c::Int) where T
    key = _sort3(a, b, c)
    pp = get(cm.probe_cache, key, nothing)
    if pp === nothing
        sa = cm.atoms[a]; sb = cm.atoms[b]; sc = cm.atoms[c]
        s1 = Sphere{T}(sa.center, sa.r + cm.pr)
        s2 = Sphere{T}(sb.center, sb.r + cm.pr)
        s3 = Sphere{T}(sc.center, sc.r + cm.pr)
        pts = intersect_three_spheres(s1, s2, s3)
        if pts === nothing
            pp = _ProbePos{T}(false, zero(Vector3{T}), zero(Vector3{T}), 0, 0)
        else
            pp = _ProbePos{T}(true, pts[1], pts[2], 0, 0)
        end
        cm.probe_cache[key] = pp
    end
    pp
end

# BALL's `checkProbe` (reducedSurface.C:1751): probe is OK iff no OTHER
# common neighbor of (a, b, c) is inside the inflated probe sphere.
function _cm_check_probe(cm::_RSComputer{T}, probe::Vector3{T},
                          a::Int, b::Int, c::Int) where T
    nb_a = cm.neighbors[a]; nb_b = cm.neighbors[b]; nb_c = cm.neighbors[c]
    # Iterate atoms in all three lists' intersection.
    for x in nb_a
        x == a || x == b || x == c || continue   # narrow down before inner
    end
    # Brute force: each common neighbor of all three.
    for x in cm.neighbors[a]
        (x == a || x == b || x == c) && continue
        x in cm.neighbors[b] || continue
        x in cm.neighbors[c] || continue
        rj = cm.atoms[x].r + cm.pr
        d = probe - Vector3{T}(cm.atoms[x].center...)
        if dot(d, d) < rj * rj - T(1e-4)
            return false
        end
    end
    true
end

# Add a new RSVertex on `atom`, return its index. Mirrors BALL
# `RSComputer::insert(RSVertex*)` (reducedSurface.C:1954-1960):
# pushes the new vertex to `new_vertices_`, registers it in
# `vertices_[atom]`, and marks the atom STATUS_ON_SURFACE.
function _cm_insert_vertex!(cm::_RSComputer{T}, atom::Int) where T
    push!(cm.rs.vertices, RSVertex(atom))
    vi = length(cm.rs.vertices)
    push!(cm.vertices_per_atom[atom], vi)
    cm.atom_status[atom] = _RS_STATUS_ON_SURFACE
    push!(cm.new_vertices, vi)
    haskey(ENV, "RS_TRACE") &&
        println(stderr, "MARK ON_SURFACE atom=$(atom-1) src=insert-vertex")
    vi
end

# Mark atom as INSIDE if currently UNKNOWN (BALL's status guard).
@inline function _cm_mark_inside!(cm::_RSComputer, atom::Int, src::AbstractString = "?")
    if cm.atom_status[atom] == _RS_STATUS_UNKNOWN
        cm.atom_status[atom] = _RS_STATUS_INSIDE
        haskey(ENV, "RS_TRACE") &&
            println(stderr, "MARK INSIDE atom=$(atom-1) src=$src")
    end
    nothing
end

# BALL `preProcessing` too-close filter (reducedSurface.C:1888-1924).
# Two atoms whose centre-to-centre distance is ≤ 5% of the larger
# radius destabilise the rolling algorithm; BALL removes one. The
# kept atom inherits the larger radius/center if the LATER atom was
# bigger. We mirror that here without changing the atom array length:
# the LATER atom is marked STATUS_INSIDE (skipped by all rolling
# code paths) and, if it was larger, its centre and radius are copied
# onto the EARLIER slot.
function _rs_preprocess_remove_too_close!(cm::_RSComputer{T}) where T
    for i in 1:length(cm.atoms)
        cm.atom_status[i] == _RS_STATUS_INSIDE && continue
        si = cm.atoms[i]
        for j in (i+1):length(cm.atoms)
            cm.atom_status[j] == _RS_STATUS_INSIDE && continue
            sj = cm.atoms[j]
            d2 = let dv = sj.center - si.center
                dot(dv, dv)
            end
            thresh = T(0.05) * max(si.r, sj.r)
            if d2 <= thresh * thresh
                # BALL: replace earlier with later if later is bigger,
                # then delete later. Match that:
                if sj.r > si.r
                    cm.atoms[i] = Sphere{T}(sj.center, sj.r)
                    cm.rs.atoms[i] = cm.atoms[i]
                    si = cm.atoms[i]
                end
                cm.atom_status[j] = _RS_STATUS_INSIDE
            end
        end
    end
    nothing
end

# `_rs_compute_with_rolling(spheres, pr)` — BALL-faithful rolling-probe RS
# construction. Mirrors `RSComputer::run()`:
# 1. Pre-process neighbor lists, init atom statuses.
# 2. Find seed (extreme-atom face / edge / vertex via getStartPosition).
# 3. getRSComponent: treatFace each face, treatEdge each unprocessed edge.
# 4. extendComponent: explore disconnected components via free edges.
#
# `correct(atom)` runs in-place inside thirdAtom whenever a degeneracy is
# detected; the resulting `CorrectAtomException` is caught at treatEdge and
# treatFace, signalling getRSComponent to reset its face iterator. No outer
# restart loop is needed.
function _rs_compute_with_rolling(spheres::Vector{Sphere{T}}, pr::T) where T
    _rs_run_once!(copy(spheres), pr)
end

# One full rolling pass with the current atoms. correct(atom) runs
# in-place; exceptions are caught at treatEdge boundaries so this
# function does not propagate them.
function _rs_run_once!(atoms::Vector{Sphere{T}}, pr::T) where T
    cm = _RSComputer{T}(atoms, pr)
    # Mark zero-radius atoms (= buried) so they're skipped throughout.
    for (i, s) in enumerate(atoms)
        s.r <= T(0.001) && (cm.atom_status[i] = _RS_STATUS_INSIDE)
    end
    # BALL `preProcessing` (reducedSurface.C:1888-1924) removes atoms
    # that are essentially coincident with another atom (`distance ≤
    # 0.05 · max(r_i, r_j)`). With stable atom indices we cannot
    # delete; instead we mark the LATER atom STATUS_INSIDE (= buried)
    # so it is skipped by all subsequent rolling/seed lookups, and if
    # the later atom is larger we transfer its centre+radius onto the
    # earlier slot (mirroring BALL line 1908-1912). The result is the
    # same in-RS topology BALL produces.
    _rs_preprocess_remove_too_close!(cm)
    # Outer loop: find seed, process component, repeat for any disconnected
    # subsurface (BALL's `while (start != 0)` in run()).
    debug = haskey(ENV, "RS_ROLL_TRACE")
    iter = 0
    while true
        iter += 1
        kind = _rs_get_start_position!(cm)
        debug && println("OUTER iter $iter: getStartPosition=$kind, faces=$(length(cm.rs.faces)) verts=$(length(cm.rs.vertices)) on_surface=$(count(==(_RS_STATUS_ON_SURFACE), cm.atom_status)) inside=$(count(==(_RS_STATUS_INSIDE), cm.atom_status))")
        # BALL `RSComputer::run` switch (reducedSurface.C:686-694):
        #   3 (face)   → getRSComponent (which itself ends with extendComponent)
        #   2 (edge)   → extendComponent
        #   1 (vertex) → no-op (the vertex was already inserted)
        #   0          → terminate the outer loop
        if kind == :face
            _rs_get_rs_component!(cm)
        elseif kind == :edge
            _rs_extend_component!(cm)
        elseif kind == :vertex
            # nothing
        else
            break
        end
        iter > 1000 && break
    end
    # BALL-faithful edge reordering. My port assigns `rs.edges` positions
    # at CREATION (e.g., when treatEdge calls `_rs_create_edge!` from a
    # newly-created face's loop), while BALL assigns indices at first
    # successful `treatEdge` (= the second face attaches). To match
    # BALL's ordering, each edge's `ball_idx` was recorded at the moment
    # it would have been inserted in BALL. Now permute rs.edges so that
    # the FIRST entry has ball_idx=1, etc. — and remap every reference.
    _rs_permute_edges_to_ball_order!(cm)
    # Final cleanup is deferred to downstream `compute_ses`; the RS is
    # ready to use.
    cm.rs
end

# Permute `rs.edges` so its order matches BALL's edge-index assignment
# order (recorded in `RSEdge.ball_idx`). Updates `vertex.edges` (Set of
# edge indices) and `face.edges` (Vector of edge indices) to use the new
# indices. Edges with `ball_idx == 0` (deleted-sentinel or never-treated)
# stay at the END in their relative creation order.
function _rs_permute_edges_to_ball_order!(cm::_RSComputer{T}) where T
    rs = cm.rs
    n = length(rs.edges)
    if haskey(ENV, "RS_PERMUTE_DEBUG")
        with_bidx = count(e -> e.ball_idx > 0, rs.edges)
        without = n - with_bidx
        alive = count(e -> e.v1 != 0, rs.edges)
        println(stderr, "RS_PERMUTE n=$n alive=$alive with_ball_idx=$with_bidx without=$without max_bidx=$(cm.ball_edge_counter[])")
    end
    # Stable-sort old indices by (ball_idx > 0 ? ball_idx : ∞).
    order = sortperm(1:n; by = i -> begin
        b = rs.edges[i].ball_idx
        b == 0 ? typemax(Int) : b
    end)
    # Build old → new mapping. order[new] = old, so we invert.
    old_to_new = zeros(Int, n)
    for new in 1:n
        old_to_new[order[new]] = new
    end
    # Permute the edges array in place.
    new_edges = similar(rs.edges)
    for new in 1:n
        new_edges[new] = rs.edges[order[new]]
    end
    rs.edges = new_edges
    # Remap vertex.edges (Set of indices).
    for v in rs.vertices
        v.atom == 0 && continue
        if !isempty(v.edges)
            new_set = Set{Int}()
            for old_ei in v.edges
                new_ei = old_to_new[old_ei]
                new_ei > 0 && push!(new_set, new_ei)
            end
            v.edges = new_set
        end
    end
    # Remap face.edges (Vector of indices).
    for f in rs.faces
        f.v1 == 0 && continue
        for k in eachindex(f.edges)
            old_ei = f.edges[k]
            f.edges[k] = old_to_new[old_ei]
        end
    end
end

# BALL's getStartPosition (reducedSurface.C:1176).
function _rs_get_start_position!(cm::_RSComputer)
    _rs_find_first_face!(cm) !== 0 && return :face
    _rs_find_first_edge!(cm) !== 0 && return :edge
    _rs_find_first_vertex!(cm) !== 0 && return :vertex
    :none
end

# BALL's findFirstFace iterates 3 directions × 2 extremes.
function _rs_find_first_face!(cm::_RSComputer)
    for direction in 1:3, extreme in 1:2
        face = _rs_find_face!(cm, direction, extreme)
        face != 0 && return face
    end
    0
end

# BALL's findFace(direction, extreme).
function _rs_find_face!(cm::_RSComputer{T}, direction::Int, extreme::Int) where T
    a1 = _rs_find_first_atom(cm, direction, extreme)
    a1 == 0 && return 0
    a2 = _rs_find_second_atom(cm, a1, direction, extreme)
    a2 == 0 && return 0
    # Find a third atom giving a valid (unoccluded, UNKNOWN-status) probe.
    common = _cm_common_neighbors(cm, a1, a2)
    isempty(common) && return 0
    a3 = 0
    probe_pt = zero(Vector3{T})
    for c in common
        (c == a1 || c == a2) && continue
        cm.atom_status[c] == _RS_STATUS_UNKNOWN || continue
        pp = _cm_probe_positions(cm, a1, a2, c)
        pp.valid || continue
        for k in 1:2
            probe = k == 1 ? pp.p0 : pp.p1
            if _cm_check_probe(cm, probe, a1, a2, c)
                a3 = c
                probe_pt = probe
                break
            end
        end
        a3 != 0 && break
    end
    if a3 == 0
        # No valid third atom — mark a1, a2 INSIDE (BALL line 1307-1308).
        cm.atom_status[a1] = _RS_STATUS_INSIDE
        cm.atom_status[a2] = _RS_STATUS_INSIDE
        if haskey(ENV, "RS_TRACE")
            println(stderr, "MARK INSIDE atom=$(a1-1) src=findFace-no-third")
            println(stderr, "MARK INSIDE atom=$(a2-1) src=findFace-no-third")
        end
        return 0
    end
    # Build the seed face's 3 vertices + 3 edges (cyclic).
    v1 = _cm_insert_vertex!(cm, a1)
    v2 = _cm_insert_vertex!(cm, a2)
    v3 = _cm_insert_vertex!(cm, a3)
    fn = _face_normal_for_triple(cm, a1, a2, a3, probe_pt)
    # BALL `RSComputer::correctProbePosition` (reducedSurface.C:832):
    # `new_face->singular_ = isLess(GetDistance(probe.p, plane), probe_radius)`.
    # The face is SINGULAR if the probe center is closer than probe_radius
    # to the plane of the 3 atom centers — i.e. the probe sphere protrudes
    # through that plane, indicating 4+ atoms could touch at the singular
    # point. Without this flag, treatFirstCategory/treatSecondCategory
    # never pair up adjacent singular spheric faces and never replace
    # contact-corner vertices with 3-probe intersection vertices on their
    # boundaries — which leaves the BST receiving slightly-wrong-position
    # boundary vertices and producing AMBIG-cascade dups downstream.
    nrm_fn = sqrt(dot(fn, fn))
    is_singular = if nrm_fn > eps(T)
        pa = Vector3{T}(cm.atoms[a1].center...)
        abs(dot(probe_pt - pa, fn)) / nrm_fn < cm.rs.probe_radius - T(1e-6)
    else
        false
    end
    push!(cm.rs.faces, RSFace{T}(v1, v2, v3, Int[], probe_pt, fn, is_singular))
    fi = length(cm.rs.faces)
    for vi in (v1, v2, v3)
        push!(cm.rs.vertices[vi].faces, fi)
    end
    # Create 3 edges; attach to face and vertices.
    _rs_create_edge!(cm, fi, v1, v2, a1, a2)
    _rs_create_edge!(cm, fi, v2, v3, a2, a3)
    _rs_create_edge!(cm, fi, v3, v1, a3, a1)
    push!(cm.new_faces, fi)
    # `_cm_insert_vertex!` already pushed v1/v2/v3 to new_vertices, mirroring
    # BALL's RSComputer::insert(RSVertex*).
    fi
end

@inline function _face_normal_for_triple(cm::_RSComputer{T}, a::Int, b::Int,
                                          c::Int, probe::Vector3{T}) where T
    pa = cm.atoms[a].center; pb = cm.atoms[b].center; pc = cm.atoms[c].center
    n = plane_normal(pa, pb, pc)
    n === nothing && return Vector3{T}(0, 0, 0)
    # BALL `getFaceNormal` (reducedSurface.C:1685-1697): negate the
    # plane normal if `norm · probe < norm · atom1`, i.e. iff
    # `n · (probe - pa) < 0`. After the conditional, `n` points TOWARD
    # the probe centre (away from the inside of the molecule). My
    # previous condition was INVERTED — pointing away from the probe.
    # That flipped the sign of every thirdAtom axis and inverted every
    # rolling-angle comparison, leaving 138 atoms wrongly buried on
    # BPTI vs BALL's reference.
    dot(n, probe - Vector3{T}(pa[1], pa[2], pa[3])) < 0 ?
        Vector3{T}(-n[1], -n[2], -n[3]) : n
end

# BALL's findFirstAtom (reducedSurface.C:1357): pick the unprocessed atom
# whose SPHERE EXTREME (center ± radius) is most extreme in the given
# direction. extreme=1 = minimum (center - radius), extreme=2 = maximum
# (center + radius).
function _rs_find_first_atom(cm::_RSComputer, direction::Int, extreme::Int)
    best = 0
    best_val = extreme == 1 ? typemax(Float64) : -typemax(Float64)
    for i in 1:length(cm.atoms)
        cm.atom_status[i] == _RS_STATUS_UNKNOWN || continue
        s = cm.atoms[i]
        v = Float64(s.center[direction]) +
            (extreme == 1 ? -Float64(s.r) : Float64(s.r))
        if extreme == 1
            v < best_val && (best_val = v; best = i)
        else
            v > best_val && (best_val = v; best = i)
        end
    end
    best
end

# BALL's findSecondAtom (reducedSurface.C:1398): pick the neighbor whose
# intersection-circle (with atom a1's inflated sphere) has the most
# extreme value in the direction. For extreme=1, BALL starts at
# `first_atom.p + first_atom.radius` (= the SAME side that findFirstAtom
# picked, INVERTED). This is the "next atom toward the extreme corner".
function _rs_find_second_atom(cm::_RSComputer{T}, a1::Int, direction::Int,
                                extreme::Int) where T
    best = 0
    s1 = cm.atoms[a1]
    inflated_r = s1.r + cm.pr
    # BALL's extrem_value initial: opposite-side extreme of a1's inflated
    # sphere (line 1421-1422). For extreme=1 (minimum sought): use the
    # MAXIMUM side of a1's inflated sphere.
    best_val = Float64(s1.center[direction]) +
               (extreme == 1 ? inflated_r : -inflated_r)
    for i in cm.neighbors[a1]
        cm.atom_status[i] == _RS_STATUS_UNKNOWN || continue
        s2 = cm.atoms[i]
        inflated_r2 = s2.r + cm.pr
        # Intersection-circle extreme: project the circle onto `direction`.
        # The circle is the boundary of two inflated spheres' intersection;
        # we approximate the extreme by the circle's center ± its radius
        # along the direction.
        # circle.p = a1.center + t * (a2.center - a1.center) / d_ab where
        # t = ((R1² - R2²) + d_ab²) / (2 d_ab).
        ab = Vector3{T}(s2.center...) - Vector3{T}(s1.center...)
        d_ab2 = dot(ab, ab)
        d_ab2 < eps(T) && continue
        d_ab = sqrt(d_ab2)
        t = ((inflated_r^2 - inflated_r2^2) + d_ab2) / (2 * d_ab)
        circle_p = Vector3{T}(s1.center...) + t * (ab / d_ab)
        # circle_radius² = inflated_r² - t²; if negative, no intersection.
        rc2 = inflated_r^2 - t^2
        rc2 < 0 && continue
        rc = sqrt(rc2)
        # Circle plane normal = ab/d_ab. Direction component of circle
        # center: circle_p[direction]. Circle radius projected onto
        # direction depends on angle between direction and circle normal.
        # Compute extremes as circle_p[direction] ± rc * sin(angle).
        n_d = ab[direction] / d_ab
        # Component of circle along `direction` axis perpendicular to ab.
        proj_radius = rc * sqrt(max(zero(T), one(T) - n_d * n_d))
        cval = Float64(circle_p[direction]) +
               (extreme == 1 ? -Float64(proj_radius) : Float64(proj_radius))
        if extreme == 1
            cval < best_val && (best_val = cval; best = i)
        else
            cval > best_val && (best_val = cval; best = i)
        end
    end
    best
end

# Create an RSEdge between two RSVertices and attach to a face. Used by
# both seed-face creation and treatEdge.
function _rs_create_edge!(cm::_RSComputer{T}, fi::Int, v1::Int, v2::Int,
                          a1::Int, a2::Int) where T
    ccs = _contact_circles(cm.rs, a1, a2)
    e = RSEdge{T}(v1, v2)
    e.f1 = fi
    if ccs !== nothing
        c1, c2, c3 = ccs
        e.center_of_torus = c1.p
        e.radius_of_torus = c1.r
        e.contact_circle1 = c2
        e.contact_circle2 = c3
    end
    push!(cm.rs.edges, e)
    ei = length(cm.rs.edges)
    push!(cm.rs.vertices[v1].edges, ei)
    push!(cm.rs.vertices[v2].edges, ei)
    push!(cm.rs.faces[fi].edges, ei)
    ei
end

# BALL `RSComputer::findFirstEdge` (reducedSurface.C:1215-1230): iterate
# the same 3 directions × 2 extremes as findFirstFace, calling `findEdge`
# at each step. Returns the new edge's index (>0) on success, 0 if none.
function _rs_find_first_edge!(cm::_RSComputer)
    for direction in 1:3, extreme in 1:2
        ei = _rs_find_edge!(cm, direction, extreme)
        ei != 0 && return ei
    end
    0
end

# BALL `RSComputer::findEdge` (reducedSurface.C:1313-1346): pick extreme
# atom a1, then extreme atom a2 in the neighbours of a1, then attempt
# `createFreeEdge(a1, a2)`. On failure, remove a2 from a1's neighbours
# (and vice versa) so the same pair is not retried.
function _rs_find_edge!(cm::_RSComputer{T}, direction::Int, extreme::Int) where T
    a1 = _rs_find_first_atom(cm, direction, extreme)
    a1 == 0 && return 0
    a2 = _rs_find_second_atom(cm, a1, direction, extreme)
    a2 == 0 && return 0
    # BALL also calls neighboursOfTwoAtoms(SortedPosition2(a1,a2)) here
    # to populate the cache; my `_cm_common_neighbors` does the same on
    # demand inside `_rs_create_free_edge`.
    fe = _rs_create_free_edge(cm, a1, a2)
    if fe === nothing
        # BALL removes the failed pair from both atoms' neighbour lists
        # (reducedSurface.C:1339-1341) so the same pair isn't retried.
        idx = findfirst(==(a2), cm.neighbors[a1])
        idx !== nothing && deleteat!(cm.neighbors[a1], idx)
        idx = findfirst(==(a1), cm.neighbors[a2])
        idx !== nothing && deleteat!(cm.neighbors[a2], idx)
        return 0
    end
    v1 = _cm_insert_vertex!(cm, a1)
    v2 = _cm_insert_vertex!(cm, a2)
    fe.v1 = v1
    fe.v2 = v2
    # BALL `findEdge` (reducedSurface.C:1357-1364): `insert(edge)` immediately
    # assigns `edge->index_ = number_of_edges_++`. Free edges (no `treatEdge`
    # cycle) get their ball-index here at creation, unlike treatment-driven
    # edges which get theirs at first `treatEdge` success.
    cm.ball_edge_counter[] += 1
    fe.ball_idx = cm.ball_edge_counter[]
    push!(cm.rs.edges, fe)
    ei = length(cm.rs.edges)
    push!(cm.rs.vertices[v1].edges, ei)
    push!(cm.rs.vertices[v2].edges, ei)
    ei
end

# BALL's findFirstVertex: isolated atom (no neighbors).
function _rs_find_first_vertex!(cm::_RSComputer)
    for i in 1:length(cm.atoms)
        cm.atom_status[i] == _RS_STATUS_UNKNOWN || continue
        isempty(cm.neighbors[i]) || continue
        _cm_insert_vertex!(cm, i)
        return i
    end
    0
end

# BALL `RSComputer::getRSComponent` (reducedSurface.C:703-732): iterate
# faces by index, processing each via treatFace. NEW faces created
# during processing are appended to rs.faces and picked up by
# subsequent iterations. If treatFace returns false (= correct(atom)
# ran and the RS state changed), BALL resets the iterator to 0 to
# re-scan from the beginning (reducedSurface.C:708-714). After all
# faces are processed, calls `extendComponent()` (line 723).
function _rs_get_rs_component!(cm::_RSComputer)
    empty!(cm.new_faces)  # not used in index-based traversal
    i = 1
    while i <= length(cm.rs.faces)
        f = cm.rs.faces[i]
        if f.v1 == 0
            i += 1
            continue
        end
        if _rs_treat_face!(cm, i)
            i += 1
        else
            i = 1   # BALL: `if (!treatFace(faces_[i])) i = 0;` (0-based)
        end
    end
    # BALL `getRSComponent` ends with `extendComponent()`.
    _rs_extend_component!(cm)
end

# BALL `RSComputer::treatFace` (reducedSurface.C:732-757). Iterates the
# face's edges; for each edge with `face_[1] == NULL` (not yet rolled),
# calls treatEdge. If treatEdge returns false (= correct() was invoked,
# topology changed) treatFace returns false too, signalling
# getRSComponent to reset its iterator. Returns true on success.
function _rs_treat_face!(cm::_RSComputer{T}, fi::Int) where T
    f = cm.rs.faces[fi]
    f.v1 == 0 && return true
    for ei in copy(f.edges)
        e = cm.rs.edges[ei]
        e.f1 != 0 && e.f2 != 0 && continue
        _rs_treat_edge!(cm, ei) || return false
    end
    true
end

# BALL `RSComputer::treatEdge` (reducedSurface.C:759-917). Returns true
# on success. Returns false (BALL: caught PROBE TOUCHES FOUR ATOMS) when
# thirdAtom called `correct()` — caller resets its face iterator.
function _rs_treat_edge!(cm::_RSComputer{T}, ei::Int) where T
    e = cm.rs.edges[ei]
    v1 = e.v1; v2 = e.v2
    (v1 == 0 || v2 == 0) && return true   # deleted edge (slot-replacement sentinel)
    e.f1 == 0 && return true              # malformed
    f1_idx = e.f1
    start_face = cm.rs.faces[f1_idx]
    a1 = cm.rs.vertices[v1].atom
    a2 = cm.rs.vertices[v2].atom
    # thirdAtom: find next atom via rolling.
    local res
    try
        res = _rs_third_atom_with_status!(cm, a1, a2, start_face.center,
                                            _other_atom_in_face(cm, start_face, a1, a2),
                                            start_face.normal)
    catch err
        err isa CorrectAtomException || rethrow()
        # BALL `treatEdge` (reducedSurface.C:782-798) catches the
        # GeneralException, returns false. The correct(atom) ran
        # in-place inside thirdAtom, so the RS state has already been
        # peeled; the caller just needs to restart.
        return false
    end
    if res === nothing
        haskey(ENV, "RS_ROLL_TRACE") && println("  treatEdge $ei: thirdAtom returned nothing for (a1=$a1, a2=$a2)")
        return true
    end
    a3, probe, ang = res
    # BALL `treatEdge` (reducedSurface.C:907): `edge->angle_ = phi`.
    # Store the rolling angle returned by thirdAtom so downstream
    # consumers (cleanToricFace, triangulation) see a non-zero arc.
    e.angle = ang
    # Build candidate face (v1, v2, new vertex3 on a3).
    v3 = _cm_insert_vertex!(cm, a3)
    new_normal = _face_normal_for_triple(cm, a1, a2, a3, probe)
    # BALL `RSComputer::correctProbePosition` (reducedSurface.C:832):
    # mark face as singular if probe protrudes through the atom-triangle
    # plane (= probe center within probe_radius of the plane).
    nrm_nn = sqrt(dot(new_normal, new_normal))
    new_singular = if nrm_nn > eps(T)
        pa = Vector3{T}(cm.atoms[a1].center...)
        abs(dot(probe - pa, new_normal)) / nrm_nn < cm.rs.probe_radius - T(1e-6)
    else
        false
    end
    new_face = RSFace{T}(v1, v2, v3, Int[], probe, new_normal, new_singular)
    push!(cm.rs.faces, new_face)
    fi_new = length(cm.rs.faces)
    push!(cm.rs.vertices[v1].faces, fi_new)
    push!(cm.rs.vertices[v2].faces, fi_new)
    push!(cm.rs.vertices[v3].faces, fi_new)
    # faceExists check: does THIS face (same atom triple + same probe center)
    # already exist among v1's incident faces? Find it.
    existing = _rs_face_exists(cm, fi_new, a1)
    haskey(ENV, "RS_ROLL_TRACE") && println("  treatEdge $ei: rolled (a1=$a1, a2=$a2) → a3=$a3, faceExists=$existing")
    if haskey(ENV, "RS_FACE_EXISTS_TRACE") && existing != 0
        ex_f = cm.rs.faces[existing]
        println(stderr, "FACE_EXISTS_HIT new_fi=$fi_new (a1=$(a1-1),a2=$(a2-1),a3=$(a3-1)) probe=($(probe[1]),$(probe[2]),$(probe[3])) matched_fi=$existing matched_center=($(ex_f.center[1]),$(ex_f.center[2]),$(ex_f.center[3]))")
    end
    if existing == 0
        # New face: link edge → new face on slot 2, create 2 new edges
        # (v2, v3) and (v3, v1) attached to new_face.
        e.f2 = fi_new
        push!(new_face.edges, ei)
        _rs_create_edge!(cm, fi_new, v2, v3, a2, a3)
        _rs_create_edge!(cm, fi_new, v3, v1, a3, a1)
        push!(cm.new_faces, fi_new)
        # `_cm_insert_vertex!` already pushed v3 to new_vertices.
    else
        # Face already exists — BALL's join + edge-replace logic
        # (reducedSurface.C:836-915):
        #   1. Find test_edge = existing's edge on the SAME atom pair as ei.
        #   2. Merge v1↔test_vertex1 and v2↔test_vertex2 (the test_edge's
        #      endpoints) — delete the existing's vertices, keep new's.
        #   3. REPLACE the slot in existing.edges: test_edge → ei.
        #   4. Delete test_edge (mark deleted).
        #   5. Set ei.f2 = existing.
        existing_f = cm.rs.faces[existing]
        # Find test_edge — the existing face's edge matching the (a1, a2)
        # atom pair of ei.
        test_eid = 0
        test_slot = 0
        for (k, ex_ei) in enumerate(existing_f.edges)
            ex_ei == ei && continue   # ei itself isn't yet in existing
            ex_e = cm.rs.edges[ex_ei]
            (ex_e.v1 == 0 || ex_e.v2 == 0) && continue   # deleted edge
            ex_a1 = cm.rs.vertices[ex_e.v1].atom
            ex_a2 = cm.rs.vertices[ex_e.v2].atom
            (ex_a1 == 0 || ex_a2 == 0) && continue   # deleted vertex
            if (ex_a1 == a1 && ex_a2 == a2) || (ex_a1 == a2 && ex_a2 == a1)
                test_eid = ex_ei
                test_slot = k
                break
            end
        end
        # Identify test_vertex1 and test_vertex2 = endpoints of test_edge
        # ordered to match (a1, a2).
        ex_v1 = ex_v2 = 0
        if test_eid != 0
            te = cm.rs.edges[test_eid]
            if cm.rs.vertices[te.v1].atom == a1
                ex_v1 = te.v1; ex_v2 = te.v2
            else
                ex_v1 = te.v2; ex_v2 = te.v1
            end
        else
            # Fallback: search existing face's verts for atom matches.
            for vi in (existing_f.v1, existing_f.v2, existing_f.v3)
                ai = cm.rs.vertices[vi].atom
                ai == a1 && (ex_v1 = vi)
                ai == a2 && (ex_v2 = vi)
            end
        end
        # Merge existing's endpoint vertices into the new ones.
        if ex_v1 != 0 && ex_v1 != v1
            _rs_join_vertices!(cm, v1, ex_v1)
        end
        if ex_v2 != 0 && ex_v2 != v2
            _rs_join_vertices!(cm, v2, ex_v2)
        end
        # Pop the freshly added face and the new vertex3. Keep all index
        # bookkeeping (vertices_per_atom AND new_vertices) consistent —
        # `_cm_insert_vertex!` queued v3 in new_vertices, so we must
        # un-queue it too, otherwise extendComponent will dereference a
        # stale index.
        pop!(cm.rs.faces)
        last_vi = length(cm.rs.vertices)
        pop!(cm.rs.vertices)
        idx = findfirst(==(last_vi), cm.vertices_per_atom[a3])
        idx !== nothing && deleteat!(cm.vertices_per_atom[a3], idx)
        nv_idx = findfirst(==(last_vi), cm.new_vertices)
        nv_idx !== nothing && deleteat!(cm.new_vertices, nv_idx)
        # Remove fi_new references from v1/v2 (we just popped fi_new).
        delete!(cm.rs.vertices[v1].faces, fi_new)
        delete!(cm.rs.vertices[v2].faces, fi_new)
        # Replace test_edge slot in existing.edges with ei; delete test_edge.
        # BALL `treatEdge` (reducedSurface.C:872-887): the old edge is
        # removed from BOTH endpoint vertices' edges sets before being
        # deleted, and the face's slot is replaced via `setEdge(i, ei)`.
        # If test_edge was also referenced by another face Y (= the
        # other side of the shared (a1, a2) edge), Y's slot must be
        # SUBSTITUTED (not removed) — ei represents the same atom-pair
        # boundary and must take the slot. Mirrors BALL `setEdge`
        # invoked transitively via the shared-edge invariant.
        if test_eid != 0 && test_eid != ei
            old_e = cm.rs.edges[test_eid]
            old_v1 = old_e.v1; old_v2 = old_e.v2
            old_v1 != 0 && delete!(cm.rs.vertices[old_v1].edges, test_eid)
            old_v2 != 0 && delete!(cm.rs.vertices[old_v2].edges, test_eid)
            for other_fi in (old_e.f1, old_e.f2)
                (other_fi == 0 || other_fi == existing) && continue
                slot = findfirst(==(test_eid), cm.rs.faces[other_fi].edges)
                # Substitute (not delete) — the atom-pair side persists,
                # only the carrier edge changes.
                slot !== nothing && (cm.rs.faces[other_fi].edges[slot] = ei)
                # And register ei as a member of Y on the edge side too
                # (ei.f1/f2 will already include `existing`; we add the
                # other neighbour Y on the still-empty slot below).
                if ei != 0 && other_fi != 0 && other_fi != existing
                    if e.f1 == 0
                        e.f1 = other_fi
                    elseif e.f2 == 0 && e.f1 != other_fi
                        e.f2 = other_fi
                    end
                end
            end
            existing_f.edges[test_slot] = ei
            old_e.v1 = 0; old_e.v2 = 0; old_e.f1 = 0; old_e.f2 = 0
        end
        # Attach existing face to ei's f2 slot.
        e.f2 = existing
    end
    # BALL `treatEdge` (reducedSurface.C:895-907): after attaching the
    # new face to this edge, recompute the toric geometry and detect
    # singular crossings (probe sphere intersects the (a1, a2) line)
    # in-flight. ip1 is the intersection closer to atom1.
    _rs_update_edge_geometry!(cm, ei, a1, a2, probe)
    # BALL `treatEdge` (reducedSurface.C:914-917): if this edge has not
    # yet been registered with the RS (`edge->index_ == -1`), assign it
    # the next ball-index now. This is BALL's actual edge-index ordering
    # — created in seed face / treatment with index=-1, finalized in this
    # spot on the FIRST successful treatment that attaches the second face.
    if e.ball_idx == 0
        cm.ball_edge_counter[] += 1
        e.ball_idx = cm.ball_edge_counter[]
    end
    true
end

# BALL `treatEdge` (reducedSurface.C:889-907) end-of-function geometry
# update. Sets contact-circle / torus geometry from the (atom1, atom2)
# pair and detects whether the new probe at `probe_pos` makes the edge
# singular (probe sphere crosses the (atom1, atom2) line). When
# singular, `intersection1` is the intersection closer to atom1.
function _rs_update_edge_geometry!(cm::_RSComputer{T}, ei::Int,
                                     atom1::Int, atom2::Int,
                                     probe_pos::Vector3{T}) where T
    e = cm.rs.edges[ei]
    ccs = _contact_circles(cm.rs, atom1, atom2)
    if ccs !== nothing
        c1, c2, c3 = ccs
        e.center_of_torus = c1.p
        e.radius_of_torus = c1.r
        e.contact_circle1 = c2
        e.contact_circle2 = c3
    end
    pa = Vector3{T}(cm.atoms[atom1].center...)
    pb = Vector3{T}(cm.atoms[atom2].center...)
    pts = intersect_sphere_line(Sphere{T}(probe_pos, cm.pr), pa, pb - pa)
    if pts === nothing
        e.singular = false
        return
    end
    ip1, ip2 = pts
    # BALL: ip1 closer to vertex_[0] (= atom1) by SQUARE-distance to atom2;
    # i.e., the one MORE distant from atom2 is closer to atom1.
    d1 = ip1 - pb
    d2 = ip2 - pb
    if dot(d1, d1) < dot(d2, d2)
        # BALL line 904: if ip1.distanceSq(sphere2) < ip2.distanceSq(sphere2), swap.
        ip1, ip2 = ip2, ip1
    end
    e.intersection1 = ip1
    e.intersection2 = ip2
    e.singular = true
    nothing
end

# BALL's RSVertex::join + RSVertex::substitute (graphVertex.h:411 +
# RSVertex.C:122): transfer all edges/faces of `from_vi` into `into_vi`,
# update all references, then delete `from_vi`.
function _rs_join_vertices!(cm::_RSComputer, into_vi::Int, from_vi::Int)
    into = cm.rs.vertices[into_vi]
    from = cm.rs.vertices[from_vi]
    atom = from.atom
    if haskey(ENV, "RS_JOIN_TRACE")
        println(stderr, "JOIN_VERT into=$into_vi (atom=$(into.atom-1)) from=$from_vi (atom=$(atom-1))")
    end
    # Move edges.
    for ei in from.edges
        e = cm.rs.edges[ei]
        if e.v1 == from_vi
            e.v1 = into_vi
        elseif e.v2 == from_vi
            e.v2 = into_vi
        end
        push!(into.edges, ei)
    end
    # Move faces.
    for fi in from.faces
        f = cm.rs.faces[fi]
        if f.v1 == from_vi
            f.v1 = into_vi
        elseif f.v2 == from_vi
            f.v2 = into_vi
        elseif f.v3 == from_vi
            f.v3 = into_vi
        end
        push!(into.faces, fi)
    end
    # Mark from_vi as deleted (atom = 0), remove from vertices_per_atom.
    from.atom = 0
    empty!(from.edges)
    empty!(from.faces)
    idx = findfirst(==(from_vi), cm.vertices_per_atom[atom])
    idx !== nothing && deleteat!(cm.vertices_per_atom[atom], idx)
    push!(cm.rm_vertices, from_vi)
    nothing
end

# BALL `GraphEdge::remove(face)` (graphEdge.h:525). Removes `face_idx`
# from the edge's face slots and returns the OTHER face (or 0 if none).
# Mirrors BALL's slot-shifting semantics: if face was in f2, just clear
# f2; if in f1, move f1 ← f2 and clear f2.
@inline function _rs_edge_remove_face!(cm::_RSComputer, ei::Int, face_idx::Int)
    e = cm.rs.edges[ei]
    if e.f2 == face_idx
        e.f2 = 0
    elseif e.f1 == face_idx
        e.f1 = e.f2
        e.f2 = 0
    end
    e.f1
end

# BALL `RSFace::remove` (RSFace.C:208). Walks the face's 3 edges; for
# each edge, if the face was its sole incident face, the edge becomes
# orphaned (added to `delete_edges`, its endpoints added to
# `test_vertices`). Otherwise detach the face from the edge and add
# the edge's OTHER incident face to `treat_faces` (it needs
# re-treatment because one of its neighbours is going away).
# Also erases the face from each of its 3 vertices' incidence sets.
function _rs_face_remove!(cm::_RSComputer, fi::Int,
                           delete_edges::Set{Int},
                           test_vertices::Set{Int},
                           treat_faces::Set{Int})
    f = cm.rs.faces[fi]
    # Detach from incident vertices.
    for vi in (f.v1, f.v2, f.v3)
        vi == 0 && continue
        delete!(cm.rs.vertices[vi].faces, fi)
    end
    # Process each edge.
    for k in eachindex(f.edges)
        ei = f.edges[k]
        ei == 0 && continue
        e = cm.rs.edges[ei]
        e.v1 == 0 && continue   # already deleted edge
        if e.f1 == fi && e.f2 == 0
            # Edge had only this face → orphan. Mark for deletion.
            push!(test_vertices, e.v1)
            push!(test_vertices, e.v2)
            push!(delete_edges, ei)
            # Detach edge from its endpoints' incidence sets.
            delete!(cm.rs.vertices[e.v1].edges, ei)
            delete!(cm.rs.vertices[e.v2].edges, ei)
            f.edges[k] = 0
        else
            other = _rs_edge_remove_face!(cm, ei, fi)
            other != 0 && push!(treat_faces, other)
        end
    end
    nothing
end

# BALL `RSComputer::correctProbePosition` (reducedSurface.C:1794-1835).
# Invalidates every cached probe position whose key triple contains
# `atom` so subsequent rolling/seed lookups recompute against the
# shrunk atom. My port stores ProbePos in `cm.probe_cache` keyed by
# sorted triples; iterate and delete matching keys.
function _rs_correct_probe_position!(cm::_RSComputer, atom::Int)
    for key in collect(keys(cm.probe_cache))
        (key[1] == atom || key[2] == atom || key[3] == atom) &&
            delete!(cm.probe_cache, key)
    end
    nothing
end

# BALL `RSComputer::correct(atom)` (reducedSurface.C:919-990). For
# every lobe-vertex on `atom`, walks its incident faces, peels them
# from the RS via `_rs_face_remove!`, deletes orphan edges, removes
# now-isolated vertices, and re-queues orphaned neighbour faces so
# `getRSComponent` picks them up again. Then shrinks the atom by
# `10*EPSILON` (= 1e-3 Å during RS construction), marks it UNKNOWN,
# and invalidates the probe cache.
function _rs_correct!(cm::_RSComputer{T}, atom::Int) where T
    if haskey(ENV, "RS_CORRECT_TRACE")
        s = cm.atoms[atom]
        println("CORRECT atom=$atom pos=($(round(s.center[1]; digits=4)),$(round(s.center[2]; digits=4)),$(round(s.center[3]; digits=4))) r=$(round(s.r; digits=4))")
    end
    # Snapshot the vertex list — _rs_face_remove! mutates incidence sets.
    vlist = copy(cm.vertices_per_atom[atom])
    for vi in vlist
        v = cm.rs.vertices[vi]
        v.atom == 0 && continue
        delete_edges = Set{Int}()
        test_vertices = Set{Int}()
        treat_faces = Set{Int}()
        # Snapshot faces — _rs_face_remove! removes from v.faces.
        faces_here = collect(v.faces)
        for fi in faces_here
            cm.rs.faces[fi].v1 == 0 && continue
            _rs_face_remove!(cm, fi, delete_edges, test_vertices, treat_faces)
        end
        # Delete the faces themselves (mark v1 = 0 sentinel).
        for fi in faces_here
            f = cm.rs.faces[fi]
            f.v1 = 0; f.v2 = 0; f.v3 = 0
            empty!(f.edges)
            delete!(treat_faces, fi)
            # remove from new_faces queue if present (BALL: new_faces_.erase)
            idx = findfirst(==(fi), cm.new_faces)
            idx !== nothing && deleteat!(cm.new_faces, idx)
        end
        # BALL `correct` (reducedSurface.C:955-961): for each face in
        # `treat_faces` (= face on the OTHER side of an edge that lost
        # one of its incident faces but wasn't itself deleted), NULL the
        # OLD slot and RE-APPEND to the END of rs_.faces_ with a NEW
        # index. This forces getRSComponent to re-treat the surviving
        # face at its new position. The face index ORDER changes —
        # any references to the old face index must be remapped to the
        # new index.
        for old_fi in treat_faces
            f = cm.rs.faces[old_fi]
            f.v1 == 0 && continue   # already deleted
            # Append the SAME face object to the end (it's a mutable
            # struct — sharing the reference keeps incidence sets and
            # edge.f1/f2 self-consistent before remapping).
            push!(cm.rs.faces, f)
            new_fi = length(cm.rs.faces)
            # Remap references: vertex.faces (Set), edge.f1/f2, and the
            # new_faces queue.
            for vi in (f.v1, f.v2, f.v3)
                vi == 0 && continue
                delete!(cm.rs.vertices[vi].faces, old_fi)
                push!(cm.rs.vertices[vi].faces, new_fi)
            end
            for ei in f.edges
                ei == 0 && continue
                e = cm.rs.edges[ei]
                e.f1 == old_fi && (e.f1 = new_fi)
                e.f2 == old_fi && (e.f2 = new_fi)
            end
            nfi = findfirst(==(old_fi), cm.new_faces)
            nfi !== nothing && (cm.new_faces[nfi] = new_fi)
            # NULL out the old slot (sentinel).
            old_f = RSFace{T}(0, 0, 0, Int[],
                              zero(Vector3{T}), zero(Vector3{T}), false)
            cm.rs.faces[old_fi] = old_f
        end
        empty!(treat_faces)
        # Delete orphan edges (v1 = 0 sentinel).
        for ei in delete_edges
            e = cm.rs.edges[ei]
            e.v1 = 0; e.v2 = 0; e.f1 = 0; e.f2 = 0
        end
        # test_vertices = vertices that lost an edge. Drop the current
        # vertex from this set (it's being removed anyway). Then for
        # each remaining test vertex with no incident edges, mark it
        # deleted.
        delete!(test_vertices, vi)
        for tv in test_vertices
            tv_v = cm.rs.vertices[tv]
            tv_v.atom == 0 && continue
            if isempty(tv_v.edges)
                tv_atom = tv_v.atom
                tv_v.atom = 0
                empty!(tv_v.faces)
                idx = findfirst(==(tv), cm.vertices_per_atom[tv_atom])
                idx !== nothing && deleteat!(cm.vertices_per_atom[tv_atom], idx)
                push!(cm.rm_vertices, tv)
                # also unqueue from new_vertices to avoid stale derefs
                nv_idx = findfirst(==(tv), cm.new_vertices)
                nv_idx !== nothing && deleteat!(cm.new_vertices, nv_idx)
            end
        end
        # Remove the lobe-vertex itself.
        v.atom = 0
        empty!(v.edges)
        empty!(v.faces)
        idx = findfirst(==(vi), cm.vertices_per_atom[atom])
        idx !== nothing && deleteat!(cm.vertices_per_atom[atom], idx)
        push!(cm.rm_vertices, vi)
        nv_idx = findfirst(==(vi), cm.new_vertices)
        nv_idx !== nothing && deleteat!(cm.new_vertices, nv_idx)
    end
    # Shrink atom by 10*EPSILON (BALL uses Constants::EPSILON = 1e-4
    # during RS construction → shrink by 1e-3 Å).
    s = cm.atoms[atom]
    cm.atoms[atom] = Sphere{T}(s.center, s.r - T(1e-3))
    # Keep the RS view consistent with the computer's view.
    cm.rs.atoms[atom] = cm.atoms[atom]
    cm.atom_status[atom] = _RS_STATUS_UNKNOWN
    _rs_correct_probe_position!(cm, atom)
    nothing
end

# BALL's faceExists: scan v1's atom's vertex list and check each vertex's
# incident faces for the same atom-triple+probe-center match.
function _rs_face_exists(cm::_RSComputer{T}, new_fi::Int, a1::Int) where T
    new_f = cm.rs.faces[new_fi]
    new_atoms = sort([cm.rs.vertices[new_f.v1].atom,
                      cm.rs.vertices[new_f.v2].atom,
                      cm.rs.vertices[new_f.v3].atom])
    new_center = new_f.center
    for vi in cm.vertices_per_atom[a1]
        for fi in cm.rs.vertices[vi].faces
            fi == new_fi && continue
            f = cm.rs.faces[fi]
            f.v1 == 0 && continue
            atoms = sort([cm.rs.vertices[f.v1].atom,
                          cm.rs.vertices[f.v2].atom,
                          cm.rs.vertices[f.v3].atom])
            atoms == new_atoms || continue
            # BALL `RSFace::operator==` (RSFace.C:106-130): atom-set match
            # AND `center_ == rsface.center_` (exact `TVector3<double>`
            # equality). The same probe-center formula evaluated on the
            # same atoms gives bit-equal results in Float64. In Float32,
            # however, two evaluations of the same probe formula can
            # differ in the LSB (~1e-7 relative), which bit-exact
            # equality treats as different faces — that floods the SES
            # with spurious extra faces (saw +28 RSFaces on BPTI/F32
            # vs F64) and produces hundreds of holes downstream.
            # Use a small absolute tolerance scaled to T's precision:
            # 100·eps(T) ≈ 1.2e-5 in F32, 2.2e-14 in F64 — tighter than
            # any geometric feature.
            tol2 = T(10000) * eps(T) * eps(T)   # |Δ|² threshold
            d = f.center - new_center
            dot(d, d) <= tol2 && return fi
        end
    end
    0
end

# Helper: get the third atom of a face given two of its atoms.
function _other_atom_in_face(cm::_RSComputer, f::RSFace, a1::Int, a2::Int)
    for vi in (f.v1, f.v2, f.v3)
        atom = cm.rs.vertices[vi].atom
        atom != a1 && atom != a2 && return atom
    end
    0
end

# BALL `RSComputer::thirdAtom` (reducedSurface.C:1099-1179). Rolls the
# probe around the (atom_a, atom_b) axis starting from `start_probe`,
# finds the next atom whose probe position gives the smallest positive
# rotation angle. Losers (atoms with strictly larger angles, or strictly
# smaller angles than later winners) are marked STATUS_INSIDE.
# Degeneracies trigger `correct()` + throw `CorrectAtomException`:
#   (a) PROBE TOUCHES FOUR ATOMS — rotation angle ≈ 0 or 2π.
#   (b) MULTIPLE ATOMS TIE — two or more atoms share the smallest angle.
# Returns (a3, probe_pos, angle) on success, or `nothing` if no
# candidate (all common neighbours skipped/start-probe-matching).
function _rs_third_atom_with_status!(cm::_RSComputer{T},
                                       atom_a::Int, atom_b::Int,
                                       start_probe::Vector3{T},
                                       third_face_atom::Int,
                                       face_normal::Vector3{T}) where T
    pa = Vector3{T}(cm.atoms[atom_a].center...)
    pb = Vector3{T}(cm.atoms[atom_b].center...)
    axis = pb - pa
    third_p = Vector3{T}(cm.atoms[third_face_atom].center...)
    test_vec = cross(face_normal, axis)
    if dot(test_vec, third_p) < dot(test_vec, pa)
        axis = -axis
    end
    ax_len2 = dot(axis, axis)
    rel_start = start_probe - pa
    proj_start = rel_start - (dot(rel_start, axis) / ax_len2) * axis
    common = _cm_common_neighbors(cm, atom_a, atom_b)
    isempty(common) && return nothing
    if haskey(ENV, "RS_TRACE_THIRD")
        println(stderr, "thirdAtom enter a1=$(atom_a-1) a2=$(atom_b-1) third_face_atom=$(third_face_atom-1) start_probe=($(start_probe[1]),$(start_probe[2]),$(start_probe[3])) axis=($(axis[1]),$(axis[2]),$(axis[3]))")
    end
    # third[]: candidates at current minimum angle (BALL's third deque).
    third = Tuple{Int, Vector3{T}, T}[]
    old_angle = T(3 * π)
    for c in common
        c == atom_a && continue
        c == atom_b && continue
        # BALL thirdAtom (reducedSurface.C:1122-1160) does NOT filter
        # candidates by atom_status — INSIDE atoms can still win and
        # become ON_SURFACE (line 1177 unconditionally promotes the
        # winner). Only `findFace` (seed search) filters by UNKNOWN.
        pp = _cm_probe_positions(cm, atom_a, atom_b, c)
        pp.valid || continue
        for probe in (pp.p0, pp.p1)
            # Skip the starting position (BALL: `k->second.p != start_probe`).
            if (probe[1] ≈ start_probe[1] && probe[2] ≈ start_probe[2] &&
                probe[3] ≈ start_probe[3])
                continue
            end
            rel_p = probe - pa
            proj_p = rel_p - (dot(rel_p, axis) / ax_len2) * axis
            new_angle = oriented_angle(proj_start, proj_p, axis)
            haskey(ENV, "RS_TRACE_THIRD") &&
                println(stderr, "  cand atom=$(c-1) probe=($(probe[1]),$(probe[2]),$(probe[3])) angle=$new_angle")
            # BALL: PROBE TOUCHES FOUR ATOMS — angle ≈ 0 or 2π.
            #   call correct(c); throw GeneralException;
            #   treatEdge catches → returns false → getRSComponent resets.
            if new_angle < T(1e-4) || new_angle > 2 * T(π) - T(1e-4)
                _rs_correct!(cm, c)
                throw(CorrectAtomException(c))
            end
            # BALL (reducedSurface.C:1136-1158) uses `<=` / `<` exactly,
            # with no ε margin. The previous port allowed a 1e-4 slack
            # which let near-tied candidates accumulate in `third` and
            # then all get marked INSIDE when a new minimum arrived —
            # over-marking by ~226 atoms on BPTI.
            if new_angle <= old_angle
                if new_angle < old_angle
                    for t in third
                        _cm_mark_inside!(cm, t[1], "thirdAtom-loser-new-min")
                    end
                    empty!(third)
                    old_angle = new_angle
                end
                push!(third, (c, probe, new_angle))
            else
                _cm_mark_inside!(cm, c, "thirdAtom-strictly-larger")
            end
        end
    end
    if length(third) > 1
        # BALL MULTIPLE ATOMS TIE (reducedSurface.C:1162-1174). pop_front
        # the first winner, then call correct() on every remaining tied
        # atom and throw the GeneralException.
        first_atom = third[1][1]
        for k in 2:length(third)
            _rs_correct!(cm, third[k][1])
        end
        throw(CorrectAtomException(first_atom))
    end
    isempty(third) && return nothing
    a3, probe, ang = third[1]
    cm.atom_status[a3] = _RS_STATUS_ON_SURFACE
    haskey(ENV, "RS_TRACE") &&
        println(stderr, "MARK ON_SURFACE atom=$(a3-1) src=thirdAtom-winner")
    return (a3, probe, ang)
end

# BALL `RSComputer::extendComponent` (reducedSurface.C:993). For each
# vertex queued in `new_vertices_`, look at its neighbour atoms that
# are STATUS_UNKNOWN; call `findThirdAtom` to enumerate all candidate
# (third-atom, probe-center) pairs; the first candidate that is
# STATUS_UNKNOWN AND passes `checkProbe` becomes a new face. If
# `findThirdAtom` returns no candidates, try `createFreeEdge` between
# atom1 and atom2 to handle disconnected components reachable only
# through a free toric edge. After a face is created, recursively
# `getRSComponent()` to roll out the new component.
function _rs_extend_component!(cm::_RSComputer{T}) where T
    while !isempty(cm.new_vertices)
        vi = popfirst!(cm.new_vertices)
        vi in cm.rm_vertices && continue
        v = cm.rs.vertices[vi]
        v.atom == 0 && continue
        atom1 = v.atom
        created_face = false
        for atom2 in cm.neighbors[atom1]
            cm.atom_status[atom2] == _RS_STATUS_UNKNOWN || continue
            candidates = _rs_find_third_atom(cm, atom1, atom2)
            if isempty(candidates)
                # BALL: try `createFreeEdge(vertex1, vertex2)` — succeeds
                # only when the rolling circle on (atom1, atom2) is
                # unobstructed by any mutual neighbour (full 2π toric
                # face, no third atom can host a probe).
                fe = _rs_create_free_edge(cm, atom1, atom2)
                if fe !== nothing
                    v2 = _cm_insert_vertex!(cm, atom2)
                    # BALL `extendComponent` (reducedSurface.C:1029): `insert(edge)`
                    # immediately. Free edges get their ball-index here.
                    cm.ball_edge_counter[] += 1
                    fe.ball_idx = cm.ball_edge_counter[]
                    push!(cm.rs.edges, fe)
                    ei = length(cm.rs.edges)
                    fe.v1 = vi; fe.v2 = v2
                    push!(cm.rs.vertices[vi].edges, ei)
                    push!(cm.rs.vertices[v2].edges, ei)
                    # BALL also re-pushes vertex1 to new_vertices_ so its
                    # other neighbours get retried on the next outer
                    # iteration.
                    push!(cm.new_vertices, vi)
                    break
                end
            else
                # Iterate candidates in BALL's order; first UNKNOWN +
                # checkProbe-passing one wins.
                for (atom3, probe) in candidates
                    cm.atom_status[atom3] == _RS_STATUS_UNKNOWN || continue
                    _cm_check_probe(cm, probe, atom1, atom2, atom3) || continue
                    v2 = _cm_insert_vertex!(cm, atom2)
                    v3 = _cm_insert_vertex!(cm, atom3)
                    nrm_vec = _face_normal_for_triple(cm, atom1, atom2, atom3, probe)
                    # BALL `RSComputer::correctProbePosition` (reducedSurface.C:832):
                    # mark singular if probe center lies within probe_radius of
                    # the plane through the three atom centers.
                    nrm_len = sqrt(dot(nrm_vec, nrm_vec))
                    is_sing = if nrm_len > eps(T)
                        pa = Vector3{T}(cm.atoms[atom1].center...)
                        abs(dot(probe - pa, nrm_vec)) / nrm_len < cm.rs.probe_radius - T(1e-6)
                    else
                        false
                    end
                    push!(cm.rs.faces, RSFace{T}(vi, v2, v3, Int[], probe,
                                                  nrm_vec,
                                                  is_sing))
                    fi = length(cm.rs.faces)
                    for vk in (vi, v2, v3)
                        push!(cm.rs.vertices[vk].faces, fi)
                    end
                    _rs_create_edge!(cm, fi, vi, v2, atom1, atom2)
                    _rs_create_edge!(cm, fi, v2, v3, atom2, atom3)
                    _rs_create_edge!(cm, fi, v3, vi, atom3, atom1)
                    push!(cm.new_vertices, vi)
                    created_face = true
                    break
                end
                created_face && break
            end
        end
        created_face && _rs_get_rs_component!(cm)
    end
    empty!(cm.rm_vertices)
end

# BALL `RSComputer::findThirdAtom` (reducedSurface.C:1244). Returns the
# full list of (third-atom, probe-center) candidates that can host the
# probe sphere on the rolling axis (atom1, atom2). NO status filter,
# NO occlusion check — those happen in the outer loop. Both probe
# positions (above/below the (a,b,c) plane) are appended per candidate.
function _rs_find_third_atom(cm::_RSComputer{T}, atom1::Int,
                              atom2::Int) where T
    common = _cm_common_neighbors(cm, atom1, atom2)
    out = Tuple{Int, Vector3{T}}[]
    for c in common
        c == atom1 && continue
        c == atom2 && continue
        pp = _cm_probe_positions(cm, atom1, atom2, c)
        pp.valid || continue
        push!(out, (c, pp.p0))
        push!(out, (c, pp.p1))
    end
    out
end

# BALL `RSComputer::createFreeEdge` (reducedSurface.C:1578-1625). Builds
# a free toric edge between atoms (atom1, atom2) if (a) the inflated
# spheres of the two atoms intersect (yielding `circle1` of nonzero
# radius), (b) the rolling circle radius exceeds the probe radius, and
# (c) no mutual neighbour's inflated sphere creates an intersection
# circle in the toric plane that overlaps `circle1`. Returns the
# fully-populated RSEdge (with v1/v2 left as 0 sentinel — caller
# assigns them) or `nothing` on failure.
function _rs_create_free_edge(cm::_RSComputer{T}, atom1::Int,
                                atom2::Int) where T
    ccs = _contact_circles(cm.rs, atom1, atom2)
    ccs === nothing && return nothing
    circle1, circle2, circle3 = ccs
    circle1.r > cm.pr || return nothing      # BALL: isGreater(circle1.radius, probe_radius)

    # Plane through circle1.p with normal circle1.n.
    n_hat = circle1.n
    p0    = circle1.p

    # For each mutual neighbour, project its inflated sphere onto the
    # toric plane and check whether the resulting circle intersects
    # `circle1`. If any does, the probe cannot roll freely → return nothing.
    common = _cm_common_neighbors(cm, atom1, atom2)
    for ni in common
        ni == atom1 && continue
        ni == atom2 && continue
        sc = Vector3{T}(cm.atoms[ni].center...)
        sr = cm.atoms[ni].r + cm.pr
        # Signed distance from sphere center to plane.
        d  = dot(sc - p0, n_hat)
        abs(d) > sr && continue             # plane misses the sphere
        # Test-circle in the plane: center = sc - d·n_hat, radius = sqrt(sr² - d²).
        tc_center = sc - d * n_hat
        tc_radius = sqrt(sr * sr - d * d)
        # BALL: two circles in the same plane intersect iff
        #       (r1 - r2)^2 ≤ |c1 - c2|^2 ≤ (r1 + r2)^2
        radius_dist  = tc_radius - circle1.r
        radius_sum   = tc_radius + circle1.r
        center_dist2 = dot(tc_center - p0, tc_center - p0)
        if radius_dist * radius_dist <= center_dist2 &&
           center_dist2 <= radius_sum * radius_sum
            return nothing
        end
    end

    # Construct the free RSEdge. Mirrors BALL's
    #   new RSEdge(v1, v2, NULL, NULL, circle1.p, circle1.radius,
    #              2π, circle2, circle3, 0, 0, false, -1)
    e = RSEdge{T}(0, 0)   # caller sets v1, v2 to the actual vertex indices
    e.f1              = 0
    e.f2              = 0
    e.center_of_torus = circle1.p
    e.radius_of_torus = circle1.r
    e.angle           = 2 * T(π)
    e.contact_circle1 = circle2
    e.contact_circle2 = circle3
    e.singular        = false
    return e
end

