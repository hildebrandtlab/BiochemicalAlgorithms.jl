export
    check_rs,
    check_ses,
    clean_ses!

# ---------------------------------------------------------------------------
# Topology validation and cleanup, mirroring BALL's
# `SolventExcludedSurface::clean()`, `::check()`, `::splitSphericFaces()`,
# `::deleteSmallToricFace()`, and the post-processing pass of
# `SESSingularityCleaner`.
#
# The Julia version is best-effort because BALL's pointer graph offers more
# fine-grained edits than our index-based one. We provide:
#
#   * `check_rs(rs)` / `check_ses(ses)`: topology integrity checks returning
#     a list of human-readable issues (empty = healthy).
#   * `clean_ses!(ses; min_toric_angle, min_torus_radius)`: drops degenerate
#     toric faces below the angular/radial threshold and any spheric/contact
#     faces that became orphaned, matching BALL's `deleteSmallToricFace`
#     and `clean[Toric|Spheric|Contact]Faces` passes.
#   * `split_spheric_faces!(ses)`: splits any spheric face whose boundary
#     decomposes into multiple disjoint loops (this happens after the
#     singularity cleaner introduces new vertices). Matches BALL's
#     `splitSphericFace`.
# ---------------------------------------------------------------------------

# Returns the unordered set of incident edges/faces at each vertex, derived
# from the canonical face/edge tables. Used to cross-check that the cached
# `vertex.edges`/`vertex.faces` are consistent with the canonical data.
function _rs_incidence(rs::ReducedSurface{T}) where T
    n_v = length(rs.vertices)
    inc_e = [Set{Int}() for _ in 1:n_v]
    inc_f = [Set{Int}() for _ in 1:n_v]
    for (i, e) in enumerate(rs.edges)
        e.v1 == 0 && continue
        push!(inc_e[e.v1], i)
        push!(inc_e[e.v2], i)
    end
    for (i, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        push!(inc_f[f.v1], i)
        push!(inc_f[f.v2], i)
        push!(inc_f[f.v3], i)
    end
    inc_e, inc_f
end

"""
    check_rs(rs::ReducedSurface) -> Vector{String}

Topology integrity check for a [`ReducedSurface`](@ref). Returns a list of
human-readable issues (empty vector = healthy). Validates:

  - every face's vertex/edge slot is within bounds and non-zero;
  - every edge's vertex slot is within bounds;
  - faces sit on three *distinct* atoms;
  - the cached `vertex.edges` / `vertex.faces` sets agree with the
    canonical face/edge tables.
"""
function check_rs(rs::ReducedSurface{T}) where T
    issues = String[]
    nv = length(rs.vertices)
    ne = length(rs.edges)
    nf = length(rs.faces)
    for (i, e) in enumerate(rs.edges)
        # e.v1 == 0 marks a deleted edge (in-place pruning).
        e.v1 == 0 && continue
        if !(1 <= e.v1 <= nv) || !(1 <= e.v2 <= nv)
            push!(issues, "edge $i has out-of-range vertex slot(s)")
        end
        if e.f1 != 0 && !(1 <= e.f1 <= nf)
            push!(issues, "edge $i references missing face $(e.f1)")
        end
        if e.f2 != 0 && !(1 <= e.f2 <= nf)
            push!(issues, "edge $i references missing face $(e.f2)")
        end
    end
    for (i, f) in enumerate(rs.faces)
        # f.v1 == 0 marks a face that was deleted (in-place RS pruning via
        # STATUS_INSIDE simulation or _delete_similar_faces!). Skip it.
        f.v1 == 0 && continue
        for vs in (f.v1, f.v2, f.v3)
            (1 <= vs <= nv) || push!(issues, "face $i vertex slot out of range")
        end
        length(f.edges) < 3 && push!(issues, "face $i has fewer than 3 edges")
        for es in f.edges
            (1 <= es <= ne) ||
                push!(issues, "face $i references missing edge $(es)")
        end
        a = rs.vertices[f.v1].atom
        b = rs.vertices[f.v2].atom
        c = rs.vertices[f.v3].atom
        (a == b || b == c || a == c) && push!(issues,
            "face $i sits on non-distinct atoms ($a, $b, $c)")
    end
    inc_e, inc_f = _rs_incidence(rs)
    for (i, v) in enumerate(rs.vertices)
        v.edges != inc_e[i] && push!(issues,
            "vertex $i edge-incidence cache disagrees with canonical data")
        v.faces != inc_f[i] && push!(issues,
            "vertex $i face-incidence cache disagrees with canonical data")
    end
    issues
end

"""
    check_ses(ses::SolventExcludedSurface) -> Vector{String}

Topology integrity check for a [`SolventExcludedSurface`](@ref). Returns a
list of issues (empty vector = healthy). Validates that:

  - every SES edge references valid SES vertex indices;
  - every face's vertex / edge list refers to in-range indices;
  - every face is one of the four `SESFaceType` enums;
  - the underlying reduced surface itself passes [`check_rs`](@ref).
"""
function check_ses(ses::SolventExcludedSurface{T}) where T
    issues = check_rs(ses.reduced_surface)
    nv = length(ses.vertices)
    ne = length(ses.edges)
    for (i, e) in enumerate(ses.edges)
        for v in (e.v1, e.v2)
            v != 0 && !(1 <= v <= nv) &&
                push!(issues, "ses edge $i has out-of-range vertex $(v)")
        end
    end
    all_faces = vcat(ses.contact_faces, ses.toric_faces, ses.spheric_faces)
    for (i, f) in enumerate(all_faces)
        for v in f.vertices
            # v == 0 is BALL's NULL-vertex marker, used by full-circle
            # SES edges on free toric faces (createFreeToricFace).
            v == 0 && continue
            (1 <= v <= nv) || push!(issues,
                "ses face $i (type $(f.type)) has out-of-range vertex $(v)")
        end
        for e in f.edges
            (1 <= e <= ne) || push!(issues,
                "ses face $i (type $(f.type)) has out-of-range edge $(e)")
        end
    end
    issues
end

# ---------------------------------------------------------------------------
# Cleanup: drop degenerate toric faces and propagate the consequences.
# ---------------------------------------------------------------------------

"""
    clean_ses!(ses; min_torus_radius = 0, min_angle = 0, density = 0.0)

Drop toric faces whose `radius_of_torus` or `angle` is below the supplied
threshold. With `density > 0`, additionally drop toric faces whose
triangulation would have fewer than 0.1 segments (BALL's
`cleanToricFace` criterion: `angle * radius * sqrt(density) < 0.1`).
Then prune any spheric/contact face whose `rs_index` references a
no-longer-present element. Mirrors BALL's `deleteSmallToricFace` +
`cleanSphericFaces` + `cleanContactFaces` cascade
(solventExcludedSurface.C:249-303).

`density` defaults to `0.0`, which **disables** the segment-count
criterion (only `min_torus_radius` / `min_angle` apply). Pass the
density you'll triangulate with to enable the BALL-style sliver
removal — sub-segment sliver toric faces are the dominant source of
T-junctions in the triangulator, so enabling it directly improves
mesh manifoldness.

Returns the number of toric faces removed.
"""
function clean_ses!(ses::SolventExcludedSurface{T};
                    min_torus_radius::Real = 0,
                    min_angle::Real = 0,
                    density::Real = 0.0) where {T<:Real}
    rs = ses.reduced_surface
    kept = SESFace{T}[]
    dropped_rs_edges = Set{Int}()
    sqrt_density = sqrt(T(density))
    use_seg_criterion = density > 0
    for f in ses.toric_faces
        # Sentinel toric face for a deleted RS edge: rs_index == 0. Keep it
        # (the deleted_face placeholder) so toric_faces[ei] stays aligned
        # with rs.edges[ei]; nothing to drop.
        if f.rs_index == 0
            push!(kept, f)
            continue
        end
        e = rs.edges[f.rs_index]
        # BALL's `cleanToricFace`: drop if the toric arc would produce
        # fewer than 0.1 segments at the given density (this criterion
        # is enabled only when `density > 0` so the default call is a
        # no-op on healthy SES). Otherwise drop only if below the
        # explicit user thresholds.
        exact_segs = e.angle * e.radius_of_torus * sqrt_density
        if e.radius_of_torus < T(min_torus_radius) ||
           e.angle < T(min_angle) ||
           (use_seg_criterion && exact_segs < T(0.1))
            push!(dropped_rs_edges, f.rs_index)
        else
            push!(kept, f)
        end
    end
    n_removed = length(ses.toric_faces) - length(kept)
    ses.toric_faces = kept
    # Spheric faces referencing an RSFace whose edges were all dropped lose
    # their boundary connection; drop them too.
    if !isempty(dropped_rs_edges)
        keep_spheric = SESFace{T}[]
        for f in ses.spheric_faces
            f.rs_index == 0 && (push!(keep_spheric, f); continue)
            rsf = rs.faces[f.rs_index]
            if !all(ei -> ei in dropped_rs_edges, rsf.edges)
                push!(keep_spheric, f)
            end
        end
        ses.spheric_faces = keep_spheric
    end
    n_removed
end

# ---------------------------------------------------------------------------
# Split spheric faces with disconnected boundaries — ports
# `SolventExcludedSurface::splitSphericFace`.
#
# A spheric face's boundary may decompose into several closed loops after
# the singularity cleaner adds intersection vertices. We walk the boundary
# by following `edge.v1`/`edge.v2` adjacencies and, if not all of the face's
# edges were visited, emit the unvisited edges as a *new* spheric face on
# the same probe sphere.
# ---------------------------------------------------------------------------

"""
    split_spheric_faces!(ses) -> Int

Walk each spheric face's edge list; if the boundary decomposes into more
than one closed loop, split the face into separate spheric faces (one per
loop). Returns the number of new faces created. Idempotent.
"""
function split_spheric_faces!(ses::SolventExcludedSurface{T}) where T
    new_faces = SESFace{T}[]
    n_new = 0
    for f in ses.spheric_faces
        # BALL `splitSphericFace` (solventExcludedSurface.C:649-657) bails
        # out if any edge has vertex_[0] == NULL (= deleted edge from the
        # small-toric cleanup pass). We mirror that: if any edge has v1=0
        # or v2=0, leave the face unsplit. Without this, my loop detector
        # treated deleted-edges' (0, 0) endpoints as real vertices and
        # produced bogus tiny loops, reducing real faces from N edges
        # down to 3.
        has_deleted = false
        for eid in f.edges
            e = ses.edges[eid]
            if e.v1 == 0 || e.v2 == 0
                has_deleted = true
                break
            end
        end
        has_deleted && continue
        loops = _spheric_face_loops(ses, f)
        # Only split when there are at least 2 CLOSED loops. Open loops
        # represent boundary chains with implicit contact-circle arcs on
        # the atom (CONVEX edges that aren't in face.edges); splitting on
        # them produces fragments that can't be triangulated.
        closed_count = count(L -> L.closed, loops)
        if closed_count > 1
            sphere = f.sphere
            rs_index = f.rs_index
            # Replace the original with the first loop in-place; append the
            # rest as new faces. (We do this by overwriting the
            # vertex/edges lists of `f` and pushing further loops.)
            first_loop = loops[1]
            empty!(f.vertices); append!(f.vertices, first_loop.vertices)
            empty!(f.edges);    append!(f.edges,    first_loop.edges)
            for L in loops[2:end]
                push!(new_faces, SESFace{T}(SESFaceType.Spheric, sphere,
                                            rs_index, L.vertices, L.edges))
                n_new += 1
            end
        end
    end
    append!(ses.spheric_faces, new_faces)
    n_new
end

# ---------------------------------------------------------------------------
# Probe-probe intersection cleanup. Mirrors BALL's `SESSingularityCleaner`
# first-category logic: two adjacent probe spheres (sharing an RS edge)
# overlap when their centers are closer than 2·probe_radius. Each such pair
# yields an intersection circle on the probe sphere; the two endpoints of
# the arc-segment overlap (in the plane of the two adjacent contact points)
# are stored as additional SES vertices and noted on both spheric faces.
# ---------------------------------------------------------------------------

"""
    resolve_probe_intersections!(ses) -> Int

For every pair of adjacent RS faces (sharing an RS edge) whose probe spheres
overlap, materialise the two intersection-circle endpoints as additional
SES vertices and register them on both spheric faces' vertex lists. Returns
the number of intersecting probe pairs handled.

This is the basic ("first-category") singularity cleanup. It does *not*
subdivide the spheric or toric faces — call [`split_spheric_faces!`](@ref)
afterwards if the resulting boundary may decompose into multiple loops.
"""
function resolve_probe_intersections!(ses::SolventExcludedSurface{T}) where T
    rs = ses.reduced_surface
    pr = rs.probe_radius
    seen = Set{NTuple{2,Int}}()
    n_pairs = 0
    # Map RSFace index → SESFace (spheric) index for fast lookup.
    rsface_to_spheric = Dict{Int,Int}()
    for (i, sf) in enumerate(ses.spheric_faces)
        rsface_to_spheric[sf.rs_index] = i
    end
    for e in rs.edges
        e.f1 == 0 && continue
        e.f2 == 0 && continue
        f1_idx = e.f1; f2_idx = e.f2
        key = f1_idx < f2_idx ? (f1_idx, f2_idx) : (f2_idx, f1_idx)
        key in seen && continue
        push!(seen, key)
        c1 = rs.faces[f1_idx].center
        c2 = rs.faces[f2_idx].center
        d_vec = c2 - c1
        d2 = dot(d_vec, d_vec)
        d2 < (2 * pr)^2 - T(1e-9) || continue
        # Probe spheres overlap. Compute the intersection circle.
        d = sqrt(d2)
        a = d / 2          # by symmetric pr=pr radius
        h2 = pr^2 - a^2
        h2 < 0 && continue
        n = d_vec / d
        circle_p = c1 + a * n
        # Pick a basis perpendicular to `n` to choose two opposite endpoints.
        # Use the axis between the two RS edge atoms as a hint so the two
        # endpoints sit in the relevant plane.
        atom_a = Vector3{T}(rs.atoms[rs.vertices[e.v1].atom].center...)
        atom_b = Vector3{T}(rs.atoms[rs.vertices[e.v2].atom].center...)
        hint = atom_b - atom_a
        u = hint - dot(hint, n) * n
        nu = sqrt(dot(u, u))
        if iszero(nu)
            # Degenerate hint — pick any perpendicular.
            u = abs(n[1]) < T(0.9) ? Vector3{T}(1,0,0) : Vector3{T}(0,1,0)
            u = u - dot(u, n) * n
            nu = sqrt(dot(u, u))
        end
        u = u / nu
        h = sqrt(h2)
        p_plus  = circle_p + h * u
        p_minus = circle_p - h * u
        push!(ses.vertices, SESVertex{T}(p_plus,  0, f1_idx, 0))
        v_plus_idx = length(ses.vertices)
        push!(ses.vertices, SESVertex{T}(p_minus, 0, f2_idx, 0))
        v_minus_idx = length(ses.vertices)
        if haskey(rsface_to_spheric, f1_idx)
            sf = ses.spheric_faces[rsface_to_spheric[f1_idx]]
            push!(sf.vertices, v_plus_idx)
            push!(sf.vertices, v_minus_idx)
        end
        if haskey(rsface_to_spheric, f2_idx)
            sf = ses.spheric_faces[rsface_to_spheric[f2_idx]]
            push!(sf.vertices, v_plus_idx)
            push!(sf.vertices, v_minus_idx)
        end
        n_pairs += 1
    end
    n_pairs
end

# Returns boundary loops of a spheric face as a vector of `(edges, vertices)`
# named tuples, where consecutive edges share a vertex.
function _spheric_face_loops(ses::SolventExcludedSurface{T},
                             f::SESFace{T}) where T
    remaining = Set{Int}(f.edges)
    loops = NamedTuple{(:edges, :vertices, :closed),
                        Tuple{Vector{Int}, Vector{Int}, Bool}}[]
    while !isempty(remaining)
        start_e = first(remaining)
        e0 = ses.edges[start_e]
        loop_edges = Int[start_e]
        loop_vertices = Int[e0.v1, e0.v2]
        delete!(remaining, start_e)
        # Walk forward
        current_tip = e0.v2
        start_tail = e0.v1
        while current_tip != start_tail
            next_e = 0
            for ei in remaining
                e = ses.edges[ei]
                if e.v1 == current_tip
                    next_e = ei; current_tip = e.v2; break
                elseif e.v2 == current_tip
                    next_e = ei; current_tip = e.v1; break
                end
            end
            next_e == 0 && break    # open boundary — bail
            push!(loop_edges, next_e)
            push!(loop_vertices, current_tip)
            delete!(remaining, next_e)
        end
        is_closed = current_tip == start_tail
        # Walk backward from start_tail if loop is open, to extend the chain
        # in the opposite direction.
        if !is_closed
            current_tail = start_tail
            while true
                next_e = 0
                for ei in remaining
                    e = ses.edges[ei]
                    if e.v1 == current_tail
                        next_e = ei; current_tail = e.v2; break
                    elseif e.v2 == current_tail
                        next_e = ei; current_tail = e.v1; break
                    end
                end
                next_e == 0 && break
                pushfirst!(loop_edges, next_e)
                pushfirst!(loop_vertices, current_tail)
                delete!(remaining, next_e)
                # Now check if reverse walk closed the loop back to current_tip
                if current_tail == current_tip
                    is_closed = true
                    break
                end
            end
        end
        push!(loops, (edges=loop_edges, vertices=loop_vertices,
                       closed=is_closed))
    end
    loops
end
