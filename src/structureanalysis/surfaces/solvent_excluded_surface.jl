export
    SolventExcludedSurface,
    compute_ses

@enumx SESFaceType Contact Toric Spheric ToricSingular
@enumx SESEdgeType Convex Concave Singular

"""
    SESVertex{T}

A vertex on the solvent-excluded surface — the point at which the rolling
probe makes contact with an atom. Each `SESVertex` records the contact
location, the atom it sits on, and (for traceability) the parent RS face or
edge that produced it.

# Fields
 - `point::Vector3{T}`     — the contact point in 3D
 - `atom::Int`             — 1-based atom index (in the source `ReducedSurface.atoms`)
 - `rs_face::Int`          — source RS face index (0 if not face-derived)
 - `rs_edge::Int`          — source RS edge index (0 if not edge-derived)
"""
struct SESVertex{T<:Real}
    point::Vector3{T}
    atom::Int
    rs_face::Int
    rs_edge::Int
end

"""
    SESEdge{T}

An edge on the SES, either *convex* (an arc on an atom — the boundary
between a contact face and a toric face), *concave* (an arc on the probe
sphere — the boundary between a spheric face and a toric face), or
*singular* (a self-intersection where 4+ atoms touch the probe).

# Fields
 - `type::SESEdgeType.T`
 - `v1::Int`, `v2::Int`    — `SESVertex` indices
 - `circle::Circle3{T}`    — carrier arc
 - `f1::Int`, `f2::Int`    — adjacent SES face indices (0 if missing)
"""
mutable struct SESEdge{T<:Real}
    type::SESEdgeType.T
    v1::Int
    v2::Int
    circle::Circle3{T}
    f1::Int
    f2::Int
end

"""
    SESFace{T}

A face on the SES. Three sub-types are represented:

 - `SESFaceType.Contact`  — patch on an atom's surface visible to the probe
   (one per `RSVertex`). `sphere` holds the atom sphere.
 - `SESFaceType.Toric`    — toroidal saddle swept by the probe between two
   atoms (one per `RSEdge`). `sphere` holds the torus's enclosing sphere
   (center + major radius).
 - `SESFaceType.Spheric`  — inverted-probe patch when the probe sits on three
   atoms (one per `RSFace`). `sphere` holds the probe sphere.
 - `SESFaceType.ToricSingular` — toric face whose probe touches a fourth
   atom; treated as a degenerate case (recorded but not split).

# Fields
 - `type::SESFaceType.T`
 - `sphere::Sphere{T}`     — the surface the face lies on
 - `rs_index::Int`         — index of the source RS element (face / edge / vertex)
 - `vertices::Vector{Int}` — SES vertex indices on this face's boundary
 - `edges::Vector{Int}`    — SES edge indices bounding this face
"""
struct SESFace{T<:Real}
    type::SESFaceType.T
    sphere::Sphere{T}
    rs_index::Int
    vertices::Vector{Int}
    edges::Vector{Int}
end

"""
    SolventExcludedSurface{T}

Analytical solvent-excluded surface (SES) built from a [`ReducedSurface`](@ref).
The SES has three distinct face types — contact (atom patches), toric
(saddle regions), and spheric (inverted-probe lobes) — connected by convex
and concave arcs.

This implementation builds the topological skeleton but does *not* split
singular toric faces (probe touching 4+ atoms). Such faces are recorded with
type `ToricSingular` and left intact, which is sufficient for triangulation
of typical molecular structures.
"""
mutable struct SolventExcludedSurface{T<:Real} <: AbstractMolecularSurface{T}
    reduced_surface::ReducedSurface{T}
    vertices::Vector{SESVertex{T}}
    edges::Vector{SESEdge{T}}
    contact_faces::Vector{SESFace{T}}
    toric_faces::Vector{SESFace{T}}
    spheric_faces::Vector{SESFace{T}}
    # Singular SES edges (indices into `edges`). Each connects the two
    # axis-intersection vertices on either side of a singular toric face
    # and bounds the two adjacent CONTACT faces (not the toric face).
    # Mirrors BALL's `SolventExcludedSurface::singular_edges_`.
    singular_edges::Vector{Int}
end

# Backwards-compat constructor: empty singular_edges if not provided.
SolventExcludedSurface{T}(rs, vs, es, cf, tf, sf) where T =
    SolventExcludedSurface{T}(rs, vs, es, cf, tf, sf, Int[])

@inline nvertices(ses::SolventExcludedSurface)     = length(ses.vertices)
@inline nedges(ses::SolventExcludedSurface)        = length(ses.edges)
# Count only alive faces (rs_index != 0), mirroring the RS-side helpers
# which skip sentinel entries (`atom == 0`, `v1 == 0`).
@inline ncontact_faces(ses::SolventExcludedSurface) =
    count(f -> f.rs_index != 0, ses.contact_faces)
@inline ntoric_faces(ses::SolventExcludedSurface)   =
    count(f -> f.rs_index != 0, ses.toric_faces)
@inline nspheric_faces(ses::SolventExcludedSurface) =
    count(f -> f.rs_index != 0, ses.spheric_faces)
@inline nfaces(ses::SolventExcludedSurface) =
    ncontact_faces(ses) + ntoric_faces(ses) + nspheric_faces(ses)

# Helpers for accumulating per-RS-element incidence during construction.
struct _SESBuilder{T<:Real}
    rs::ReducedSurface{T}
    vertices::Vector{SESVertex{T}}
    edges::Vector{SESEdge{T}}
    # rs_face_to_contact_points[i][k] = SESVertex index where the probe of
    # RSFace[i] touches the kth atom of that face (k ∈ 1:3).
    rs_face_to_contact_points::Vector{Vector{Int}}
end

_SESBuilder{T}(rs::ReducedSurface{T}) where T = _SESBuilder{T}(
    rs,
    SESVertex{T}[],
    SESEdge{T}[],
    # Size by raw array length (NOT `nfaces`, which only counts alive faces)
    # because the iteration loop uses raw RS-face indices.
    [Int[] for _ in 1:length(rs.faces)],
)

# Where the probe at RSFace.center touches the given atom: along the line
# from probe→atom, at distance probe_radius from the probe center (= distance
# atom_radius from the atom center).
function _probe_contact_point(rs::ReducedSurface{T}, probe::Vector3{T}, atom_idx::Int) where T
    # BALL `SESComputer::getPoint` (solventExcludedSurface.C:1562-1573)
    # used to compute contact corner positions (line 1214: `getPoint(atom->p,
    # probe_center, atom->radius, vertex->point_)`):
    #   result = atom.center + atom.radius * (probe - atom.center) / |probe - atom.center|
    # i.e. ATOM-centered, scaled by atom.radius — not probe-centered with
    # probe_radius. The two formulas agree when probe and atom are exactly
    # touching (|d| = probe_radius + atom.radius), but differ by ~epsilon
    # when RS construction leaves residual numerical drift in the distance.
    # Match BALL exactly.
    a = rs.atoms[atom_idx]
    ac = Vector3{T}(a.center...)
    dir = probe - ac
    ℓ = norm(dir)
    iszero(ℓ) && return ac
    ac + dir * (a.r / ℓ)
end

"""
    $(TYPEDSIGNATURES)

Build the analytical solvent-excluded surface (SES) from a precomputed
[`ReducedSurface`](@ref).

Mirrors BALL `SESComputer::run` (solventExcludedSurface.C:1045-1058):
build the SES once, run the singularity cleaner, and if the cleaner
called `deleteSimilarFaces` on the underlying RS (returning false),
discard the SES and rebuild it from the now-modified RS. Loops until
the cleaner reports no further modifications.
"""
function compute_ses(rs::ReducedSurface{T}; density::Real = 1.0) where {T<:Real}
    # BALL `SESComputer::run` (solventExcludedSurface.C:1045-1058) loops
    # `while (!sessc.run())` without a max-iter cap. The 9-edge
    # deleteSimilarFaces case modifies the RS so the next iteration sees
    # a different topology and (eventually) returns true. Match BALL.
    local ses
    while true
        # BALL `SESComputer::run` shares `vertex_grid_` across all
        # singular-vertex creation paths so an intersection point seen
        # at toric-face construction (during `_compute_ses_once`) is
        # reused when `_treat_singular_toric_faces!` adds the same
        # corner via `createSingularVertex`. Rebuilt per restart since
        # the SES is rebuilt too.
        vertex_grid = Dict{NTuple{3,Int}, Vector{Int}}()
        ses = _compute_ses_once(rs; density, vertex_grid)
        _treat_singular_toric_faces!(ses, vertex_grid)
        if _ses_treat_first_category!(ses)
            # Cleaner ran clean (no 9-edge deleteSimilarFaces fired).
            _ses_treat_second_category!(ses)
            return ses
        end
        # 9-edge case fired and modified the RS. BALL: clear SES, rebuild.
        # No explicit `rs.clean()` is needed here — my `_merge_similar_rs_faces!`
        # already marks the affected faces/edges/vertices as deleted in
        # place, and `_compute_ses_once` skips deleted slots.
    end
end

# Inner builder — pure SES construction with no singularity cleaning.
# Splitting this out lets `compute_ses` invoke it in the restart loop
# that BALL's `SESComputer::run` runs around `SESSingularityCleaner`.
# The optional `vertex_grid` is the dedup hash shared with
# `_treat_singular_toric_faces!`; pass an empty Dict to skip dedup.
function _compute_ses_once(rs::ReducedSurface{T}; density::Real = 1.0,
                            vertex_grid::Dict{NTuple{3,Int}, Vector{Int}}
                                = Dict{NTuple{3,Int}, Vector{Int}}()) where {T<:Real}
    b = _SESBuilder{T}(rs)

    # ----- For each RS face, create 3 contact-point vertices and 3 concave edges
    # bounding the spheric face. -----
    # BALL `pushVertex` (solventExcludedSurface.C:1150-1178) registers every
    # SES vertex (including contact-point vertices) in the shared
    # `vertex_grid_` as it's created, so `vertexExists` in `createSingularVertex`
    # can find coincidences with existing contact points. My port had been
    # registering only singular-intersection vertices, which meant a
    # downstream `_b_create_or_find_singular_vertex!` could create a
    # duplicate vertex at a position that coincides with an existing
    # contact point — splitting topology that BALL would have shared.
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue   # face was deleted by _delete_similar_faces!
        v_idx = ntuple(k -> begin
            atom_idx = rs.vertices[(f.v1, f.v2, f.v3)[k]].atom
            point = _probe_contact_point(rs, f.center, atom_idx)
            push!(b.vertices, SESVertex{T}(point, atom_idx, fi, 0))
            vidx = length(b.vertices)
            _insert_vertex_grid!(vertex_grid, point, vidx)
            vidx
        end, 3)
        b.rs_face_to_contact_points[fi] = collect(v_idx)
    end

    # ----- Spheric faces (one per RSFace). Slot-style construction, per
    # BALL's `SESComputer::pushConcaveEdge` (solventExcludedSurface.C:1182):
    # for each atom-pair SIDE of the RS face (3 sides total), BALL calls
    # `rsface->getEdge(v1, v2, rsedge)` which returns ONE RS edge per
    # vertex pair. So BALL also creates exactly 3 concave SES edges per
    # spheric face. Stacked-probe RS edges sharing an atom-pair side
    # leave one toric face with only 1 concave + 2 convex = 3 edges; BALL
    # eliminates that scenario upstream via `correct(atom)` (not yet
    # ported). My earlier multi-edge variant overgenerates compared to
    # BALL and is reverted.
    # `spheric_faces` is INDEX-ALIGNED with `rs.faces` (size = length, not
    # nfaces): elsewhere the code indexes `ses.spheric_faces[e.f1]` directly,
    # so deleted-RS-face slots need a sentinel placeholder (rs_index=0).
    sentinel_sphere = Sphere{T}(Vector3{T}(0, 0, 0), zero(T))
    spheric_faces = SESFace{T}[
        SESFace{T}(SESFaceType.Spheric, sentinel_sphere, 0, Int[], Int[])
        for _ in 1:length(rs.faces)
    ]
    for (fi, f) in enumerate(rs.faces)
        f.v1 == 0 && continue
        v123 = b.rs_face_to_contact_points[fi]
        probe_sphere = Sphere{T}(f.center, rs.probe_radius)
        face_atoms = (rs.vertices[f.v1].atom,
                      rs.vertices[f.v2].atom,
                      rs.vertices[f.v3].atom)
        edge_idx = Int[]
        seen_pairs = Set{Tuple{Int,Int}}()
        for k in 1:3
            va_slot = k; vb_slot = mod1(k + 1, 3)
            atom_a = face_atoms[va_slot]
            atom_b = face_atoms[vb_slot]
            pair_key = atom_a < atom_b ? (atom_a, atom_b) : (atom_b, atom_a)
            pair_key in seen_pairs && continue
            push!(seen_pairs, pair_key)
            ei = 0
            for cand in f.edges
                re = rs.edges[cand]
                # BALL `face->edge_[i]` returns NULL for deleted edges
                # after correct(atom) / setEdge slot replacement. My
                # port marks deleted edges with v1=v2=0 sentinel; skip
                # them.
                re.v1 == 0 && continue
                a1 = rs.vertices[re.v1].atom
                a2 = rs.vertices[re.v2].atom
                if (a1 == atom_a && a2 == atom_b) || (a1 == atom_b && a2 == atom_a)
                    ei = cand; break
                end
            end
            ei == 0 && continue
            va = v123[va_slot]; vb = v123[vb_slot]
            pa = b.vertices[va].point
            pb = b.vertices[vb].point
            normal = cross(pa - f.center, pb - f.center)
            ℓ = norm(normal)
            normal_n = ℓ > eps(T) ? normal / ℓ : Vector3{T}(0, 0, 1)
            circle = Circle3{T}(f.center, normal_n, rs.probe_radius)
            push!(b.edges, SESEdge{T}(SESEdgeType.Concave, va, vb, circle, fi, ei))
            push!(edge_idx, length(b.edges))
        end
        spheric_faces[fi] = SESFace{T}(SESFaceType.Spheric, probe_sphere, fi,
                                        copy(v123), edge_idx)
    end

    # ----- Toric faces (one per RSEdge). Each toric face is bounded by two
    # convex edges (arcs on the two atoms) and two concave edges (already
    # created above). For a singular RSEdge the toric face is marked. -----
    toric_faces = SESFace{T}[]
    sizehint!(toric_faces, nedges(rs))
    # Sentinel face used for deleted RS edges so toric_faces[ei] stays
    # aligned with rs.edges[ei] indexing.
    deleted_face = SESFace{T}(SESFaceType.Toric,
                              Sphere{T}(Vector3{T}(0,0,0), zero(T)),
                              0, Int[], Int[])
    for (ei, e) in enumerate(rs.edges)
        if e.v1 == 0
            push!(toric_faces, deleted_face)
            continue
        end
        a = rs.vertices[e.v1].atom
        bb = rs.vertices[e.v2].atom
        type = e.singular ? SESFaceType.ToricSingular : SESFaceType.Toric
        # Torus's enclosing sphere — for triangulation we only need the
        # center (torus center) and an upper-bound radius (major + probe).
        sphere = Sphere{T}(rs.atoms[a].center, rs.atoms[a].r + rs.probe_radius)
        # Vertices on this toric face are the four contact points where the
        # two adjacent RS faces' probes touch atoms a and b.
        v_list = Int[]
        for fi in (e.f1, e.f2)
            fi == 0 && continue
            f = rs.faces[fi]
            f.v1 == 0 && continue
            for (k, vidx) in enumerate((f.v1, f.v2, f.v3))
                atom_k = rs.vertices[vidx].atom
                (atom_k == a || atom_k == bb) || continue
                push!(v_list, b.rs_face_to_contact_points[fi][k])
            end
        end
        # Singular edges: add the two line/probe-sphere intersection points
        # as additional SES vertices on this face. They mark where the toric
        # saddle "pinches" and the face splits into multiple lobes.
        if e.singular
            # BALL `createSingularVertex` semantics (solventExcludedSurface.C:1458):
            # check vertex_grid for an existing SES vertex at this point
            # (tolerance 1e-3 Å), reuse if found. Two adjacent singular
            # toric faces can share an intersection point — BALL dedupes
            # via the grid; we mirror that.
            ip1 = _b_create_or_find_singular_vertex!(b, vertex_grid,
                                                       e.intersection1, a, 0, ei)
            push!(v_list, ip1)
            ip2 = _b_create_or_find_singular_vertex!(b, vertex_grid,
                                                       e.intersection2, bb, 0, ei)
            push!(v_list, ip2)
        end
        push!(toric_faces, SESFace{T}(type, sphere, ei, v_list, Int[]))
    end

    # ----- Contact faces (one per RSVertex). For each atom that participates
    # in the surface, a contact face is the part of the atom sphere not
    # occluded by toric/spheric pieces. Vertices: all contact points on this
    # atom. Edges: convex edges connecting consecutive contact points. -----
    contact_faces = SESFace{T}[]
    sizehint!(contact_faces, length(rs.vertices))
    for (vi, rs_v) in enumerate(rs.vertices)
        atom_idx = rs_v.atom
        if atom_idx == 0
            push!(contact_faces, deleted_face)
            continue
        end
        atom = rs.atoms[atom_idx]
        sphere = Sphere{T}(atom.center, atom.r)
        # Collect SES vertices on this atom (from each incident RS face).
        v_list = Int[]
        for fi in rs_v.faces
            f = rs.faces[fi]
            f.v1 == 0 && continue
            for (k, vidx) in enumerate((f.v1, f.v2, f.v3))
                if rs.vertices[vidx].atom == atom_idx
                    push!(v_list, b.rs_face_to_contact_points[fi][k])
                end
            end
        end
        # Convex edges: one per RS edge incident to this vertex. The carrier
        # circle is the contact circle (atom_radius scaled) from the RS edge.
        # BALL `createConvexEdge` (solventExcludedSurface.C:1342) sets the
        # circle as `(circle_own.p, circle_own.p - circle_other.p, circle_own.radius)`,
        # so the normal points FROM the neighbor atom's contact circle TOWARD
        # the owning atom's contact circle. Our `re.contact_circle1/2` inherit
        # `c1.n` (atom→neighbor) from `_contact_circles`, which is the OPPOSITE
        # sign. Reconstruct the circle with BALL's sign here so sampling along
        # the arc takes the short way (otherwise `ball_partition_of_circle`
        # walks 2π - arc, sampling the long way around).
        edge_idx = Int[]
        # BALL `SESComputer::createConvexEdge` (solventExcludedSurface.C:1313)
        # is invoked from the toric-face iteration loop, which iterates
        # `toric_faces_` in INDEX order = RS-edge-INDEX order. Convex edges
        # get pushed to `contact_faces_[atom]->edge_` in RS-edge-index order.
        # `_rs_permute_edges_to_ball_order!` aligns my `rs.edges` indices
        # with BALL's, so sorting `rs_v.edges` (a Set, randomized iteration)
        # by index gives BALL's exact contact-face edge ordering — which is
        # what `firstSESEdge` (BFT seed selection) depends on.
        for ei in sort!(collect(rs_v.edges))
            re = rs.edges[ei]
            re.v1 == 0 && continue   # edge deleted by _delete_similar_faces!
            # Which contact-circle slot belongs to this atom?
            own_circle = rs.vertices[re.v1].atom == atom_idx ? re.contact_circle1 :
                                                                re.contact_circle2
            other_circle = rs.vertices[re.v1].atom == atom_idx ? re.contact_circle2 :
                                                                  re.contact_circle1
            normal = own_circle.p - other_circle.p
            ℓ = sqrt(dot(normal, normal))
            normal_n = ℓ > eps(T) ? normal / ℓ : Vector3{T}(0, 0, 1)
            target_circle = Circle3{T}(own_circle.p, normal_n, own_circle.r)
            # Boundary vertices: the two SES vertices on this atom from re's
            # two adjacent RS faces.
            two = Int[]
            for fi in (re.f1, re.f2)
                fi == 0 && continue
                f = rs.faces[fi]
                f.v1 == 0 && continue
                for (k, vidx) in enumerate((f.v1, f.v2, f.v3))
                    if rs.vertices[vidx].atom == atom_idx
                        push!(two, b.rs_face_to_contact_points[fi][k])
                    end
                end
            end
            if haskey(ENV, "SES_TRI_TRACE") && vi == 71
                println("  cf71 loop: edge $ei, two=$two, length(two)=$(length(two))")
            end
            length(two) >= 2 || continue
            # f1 = contact face index (= vi, since contact_faces[vi] is for
            # RSVertex vi). f2 = toric face index for RS edge ei.
            v0_idx, v1_idx = two[1], two[2]
            v0_p = b.vertices[v0_idx].point
            v1_p = b.vertices[v1_idx].point
            # BALL `createConvexEdge` (solventExcludedSurface.C:1352-1362):
            # check whether sampling from vertex[0] around circle.n reaches
            # vertex[1] in the same angular direction as the rsedge's
            # rolling angle. If not, REVERT (swap vertex[0]/vertex[1]).
            test_phi = oriented_angle(v0_p - target_circle.p,
                                       v1_p - target_circle.p,
                                       target_circle.n)
            if (test_phi - T(π)) * (re.angle - T(π)) < zero(T)
                v0_idx, v1_idx = v1_idx, v0_idx
            end
            push!(b.edges, SESEdge{T}(SESEdgeType.Convex, v0_idx, v1_idx,
                                       target_circle, vi, ei))
            push!(edge_idx, length(b.edges))
        end
        push!(contact_faces, SESFace{T}(SESFaceType.Contact, sphere, vi,
                                         v_list, edge_idx))
    end

    # Post-pass: populate the toric face edge lists. Each toric face is
    # bounded by 2 concave edges (one from each adjacent RSFace) and 2
    # convex edges (one per atom of the RS edge). We've stored f2 = ei
    # on those edges, so a single linear scan suffices.
    for (e_idx, se) in enumerate(b.edges)
        if se.type == SESEdgeType.Concave || se.type == SESEdgeType.Convex
            tf_idx = se.f2  # toric face index = RS edge index
            tf_idx == 0 && continue
            tf_idx <= length(toric_faces) || continue
            push!(toric_faces[tf_idx].edges, e_idx)
        end
    end
    # NOTE: a second-pass for sharing concave SES edges across stacked-
    # probe toric faces was tried (so each stacked toric has 4 edges) but
    # it caused double partitioning of the same SES edge by both toric
    # triangulators and increased overall hole count. Faithful handling
    # requires sharing the partition via `tr.edge[se_idx]` — see TODO in
    # ses_triangulator.jl. Until then we leave stacked toric faces at 3
    # edges; the triangulator skips them, leaving their associated
    # contact-side partitioning at 0.

    ses = SolventExcludedSurface{T}(rs, b.vertices, b.edges,
                                     contact_faces, toric_faces, spheric_faces,
                                     Int[])
    # BALL's `createFreeToricFace` (solventExcludedSurface.C:1525): for
    # each toric face whose RS edge has at least one missing RSFace
    # (= isolated 2-atom configuration), construct the 2 full-circle
    # convex SES edges on atoms a and b. These have v1=v2=0 (no SES
    # corner vertices — the boundary wraps fully around the rolling
    # axis). The triangulator dispatches free tori through
    # `ball_triangulate_free_toric!`.
    _create_free_toric_edges!(ses)
    # BALL's SESFace::normalize: order each face's edges (and vertices)
    # in CCW cyclic boundary order so the SES-driven triangulator can
    # walk them as `edge[0] → edge[1] → ...`. Mirrors what BALL's
    # `triangulateNonSingularToricFace` and `treatSingularToricFace`
    # both assume about the input face structure.
    _normalize_all_ses_faces!(ses)
    # BALL's `SolventExcludedSurface::clean()` (solventExcludedSurface.C:152):
    # delete small toric faces and convert one of their convex edges to a
    # singular SES edge between adjacent spheric faces. Runs AFTER
    # normalize so edges are in canonical concave-convex-concave-convex
    # order.
    # BALL `SolventExcludedSurface::compute()` does NOT call `clean()`.
    # The triangulator's preProcessing calls `ses_->clean(density)`
    # (triangulatedSES.C:188); we mirror that — `triangulate_ses` runs
    # `_ses_clean_small_toric_faces!` itself. Running it here would
    # delete ~78 toric faces on BPTI that BALL keeps in the SES
    # post-compute, breaking parity.
    # Singularity cleaning (`_treat_singular_toric_faces!`,
    # `_ses_treat_first_category!`, `_ses_treat_second_category!`) is
    # invoked by the outer `compute_ses` driver so the 9-edge
    # deleteSimilarFaces restart loop has somewhere to retry from.
    ses
end

# ---------------------------------------------------------------------------
# Small graph primitives matching BALL's `GraphFace::getEdges` /
# `GraphFace::getEdge` / `GraphEdge::other` (graphFace.h:395-446,
# graphEdge.h:437). Used by `_ses_sort_two_cuts!` and `_ses_two_cuts!`.

# BALL `GraphFace::getEdges` (graphFace.h:395): scan `face.edges` for the
# two SES edges that are incident to `vidx`. Returns (e1, e2); zero(s) if
# fewer than two are found.
function _ses_face_get_edges(ses::SolventExcludedSurface, face_edges::Vector{Int},
                              vidx::Int)
    e1 = 0; e2 = 0
    for eid in face_edges
        e = ses.edges[eid]
        if e.v1 == vidx || e.v2 == vidx
            if e1 == 0
                e1 = eid
            else
                e2 = eid
                break
            end
        end
    end
    e1, e2
end

# BALL `GraphFace::getEdge` (graphFace.h:429): scan `face.edges` for the
# edge connecting (v1, v2). Returns 0 if none.
function _ses_face_get_edge(ses::SolventExcludedSurface,
                             face_edges::Vector{Int},
                             v1::Int, v2::Int)
    for eid in face_edges
        e = ses.edges[eid]
        ((e.v1 == v1 && e.v2 == v2) || (e.v1 == v2 && e.v2 == v1)) && return eid
    end
    0
end

# BALL `GraphEdge::other` (graphEdge.h:437): return the OTHER endpoint
# of `edge` given one endpoint.
@inline _ses_edge_other(e::SESEdge, vidx::Int) =
    e.v1 == vidx ? e.v2 : e.v1

# Port of BALL `SESSingularityCleaner::sort` (solventExcludedSurface.C:1858).
# Aligns the 7 edges and 7 vertices of two spheric faces (face1, face2)
# that share two atom triple but have intersecting probe spheres. The
# canonical ordering: sesvertex1[0] / sesvertex2[0] is a shared vertex;
# walking the face boundary visits sesvertex[1..6] in order; sesedge[i]
# connects sesvertex[i] and sesvertex[i+1] (modulo 7 for sesedge[6] which
# connects sesvertex[6] and sesvertex[0]).
function _ses_sort_two_cuts!(ses::SolventExcludedSurface,
                              face1::SESFace, face2::SESFace,
                              sesedge1::Vector{Int}, sesedge2::Vector{Int},
                              sesvertex1::Vector{Int}, sesvertex2::Vector{Int})
    # find two equal vertices (BALL line 1869-1880)
    shared_v = 0
    for v1 in face1.vertices
        if v1 in face2.vertices
            shared_v = v1
            break
        end
    end
    shared_v == 0 && error("_ses_sort_two_cuts!: no shared vertex between faces")
    sesvertex1[1] = shared_v
    sesvertex2[1] = shared_v
    # find first corresponding edges (BALL line 1882-1903)
    e1a, e1b = _ses_face_get_edges(ses, face1.edges, shared_v)
    e2a, e2b = _ses_face_get_edges(ses, face2.edges, shared_v)
    sesedge1[1] = e1a
    sesedge2[1] = e2a
    # BALL: if sesedge1[0] == sesedge2[1] → sesedge1[0] = sesedge1[1]
    #       else if sesedge1[1] == sesedge2[0] → sesedge2[0] = sesedge2[1]
    #       else if sesedge1[0] == sesedge2[0] → both swap to the other
    if e1a == e2b
        sesedge1[1] = e1b
    elseif e1b == e2a
        sesedge2[1] = e2b
    elseif e1a == e2a
        sesedge1[1] = e1b
        sesedge2[1] = e2b
    end
    # find remaining edges and vertices (BALL line 1905-1922)
    sesvertex1[2] = _ses_edge_other(ses.edges[sesedge1[1]], sesvertex1[1])
    sesvertex2[2] = _ses_edge_other(ses.edges[sesedge2[1]], sesvertex2[1])
    for i in 2:6
        ea, eb = _ses_face_get_edges(ses, face1.edges, sesvertex1[i])
        sesedge1[i] = (ea != sesedge1[i-1]) ? ea : eb
        ea, eb = _ses_face_get_edges(ses, face2.edges, sesvertex2[i])
        sesedge2[i] = (ea != sesedge2[i-1]) ? ea : eb
        sesvertex1[i+1] = _ses_edge_other(ses.edges[sesedge1[i]], sesvertex1[i])
        sesvertex2[i+1] = _ses_edge_other(ses.edges[sesedge2[i]], sesvertex2[i])
    end
    sesedge1[7] = _ses_face_get_edge(ses, face1.edges, sesvertex1[1], sesvertex1[7])
    sesedge2[7] = _ses_face_get_edge(ses, face2.edges, sesvertex2[1], sesvertex2[7])
    # BALL line 1924-1942: if sesvertex1[2] != sesvertex2[2], reverse
    # the orientation of both vectors (swap positions 0↔6, 1↔5, 2↔4 in
    # vertices; 0↔5, 1↔4, 2↔3 in edges — using BALL's 0-based indices
    # 6-i with i=0..2 = 6,5,4 ↔ 0,1,2). In 1-based terms: 1↔7, 2↔6, 3↔5
    # for vertices; 1↔6, 2↔5, 3↔4 for edges.
    if sesvertex1[3] != sesvertex2[3]
        for i in 1:3
            sesvertex1[i], sesvertex1[8-i] = sesvertex1[8-i], sesvertex1[i]
            sesvertex2[i], sesvertex2[8-i] = sesvertex2[8-i], sesvertex2[i]
            sesedge1[i], sesedge1[7-i] = sesedge1[7-i], sesedge1[i]
            sesedge2[i], sesedge2[7-i] = sesedge2[7-i], sesedge2[i]
        end
    end
    nothing
end

# Port of BALL `SESSingularityCleaner::twoCuts`
# (solventExcludedSurface.C:1789-1855). The 7-edge case: face1 and face2
# share two SES vertices and several edges around them. We sort
# vertex/edge order, compute the intersection circle of the two probe
# spheres, and create two new SINGULAR SES edges connecting
# (sesvertex1[1], sesvertex1[3]) and (sesvertex1[4], sesvertex1[7]).
# Then delete the shared centre edges if they are the SAME SES edge in
# both faces (sesedge1[3] == sesedge2[3], sesedge1[7] == sesedge2[7]).
function _ses_two_cuts!(ses::SolventExcludedSurface{T},
                         si1::Int, si2::Int) where T
    face1 = ses.spheric_faces[si1]
    face2 = ses.spheric_faces[si2]
    sesedge1   = zeros(Int, 7)
    sesedge2   = zeros(Int, 7)
    sesvertex1 = zeros(Int, 7)
    sesvertex2 = zeros(Int, 7)
    _ses_sort_two_cuts!(ses, face1, face2,
                         sesedge1, sesedge2, sesvertex1, sesvertex2)
    rs = ses.reduced_surface
    pr = rs.probe_radius
    rsf1 = rs.faces[face1.rs_index]
    rsf2 = rs.faces[face2.rs_index]
    c1 = rsf1.center
    c2 = rsf2.center
    # Probe-probe intersection circle (sphere1 = (c1, pr), sphere2 = (c2, pr)).
    d  = c2 - c1
    ℓ2 = dot(d, d)
    ℓ2 > eps(T) || return
    ℓ  = sqrt(ℓ2)
    ℓ >= 2 * pr && return
    a = ℓ / 2
    r2 = pr * pr - a * a
    r2 > zero(T) || return
    r_int = sqrt(r2)
    circle_p = c1 + (a / ℓ) * d
    circle_n = d / ℓ
    # BALL: if oriented angle from sesvertex1[0] to sesvertex1[2] around
    # circle.n > π, negate circle.n.
    p0 = ses.vertices[sesvertex1[1]].point
    p2 = ses.vertices[sesvertex1[3]].point
    phi = oriented_angle(p0 - circle_p, p2 - circle_p, circle_n)
    if phi > T(π)
        circle_n = -circle_n
    end
    new_circle = Circle3{T}(circle_p, circle_n, r_int)
    # New singular edge 1 connects sesvertex1[1] and sesvertex1[3];
    # bordered by face1 (= f1 slot) and face2 (= f2 slot).
    push!(ses.edges, SESEdge{T}(SESEdgeType.Singular,
                                   sesvertex1[1], sesvertex1[3],
                                   new_circle, si1, si2))
    eid_new1 = length(ses.edges)
    push!(ses.singular_edges, eid_new1)
    push!(face1.edges, eid_new1)
    push!(face2.edges, eid_new1)
    # New singular edge 2 connects sesvertex1[4] and sesvertex1[7].
    push!(ses.edges, SESEdge{T}(SESEdgeType.Singular,
                                   sesvertex1[4], sesvertex1[7],
                                   new_circle, si1, si2))
    eid_new2 = length(ses.edges)
    push!(ses.singular_edges, eid_new2)
    push!(face1.edges, eid_new2)
    push!(face2.edges, eid_new2)
    # Delete the shared centre edges if they coincide across the two
    # faces (BALL: sesedge1[2] == sesedge2[2], sesedge1[6] == sesedge2[6]
    # in 0-based; mine = sesedge1[3] == sesedge2[3], sesedge1[7] ==
    # sesedge2[7]).
    for (k_vert_a, k_vert_b, idx_a, idx_b) in
            ((3, 4, 3, 3), (7, 1, 7, 7))
        eA = sesedge1[idx_a]
        eB = sesedge2[idx_b]
        eA == 0 && continue
        if eA == eB
            _ses_drop_edge_from_face!(face1, eA)
            _ses_drop_edge_from_face!(face2, eA)
            e = ses.edges[eA]
            # remove from endpoints' incidence (no per-vertex edge set in
            # SES; nothing to detach there)
            e.v1 = 0; e.v2 = 0; e.f1 = 0; e.f2 = 0
            # also remove from ses.singular_edges if present
            idx = findfirst(==(eA), ses.singular_edges)
            idx !== nothing && deleteat!(ses.singular_edges, idx)
        end
    end
    nothing
end

# Helper: remove first occurrence of `eid` from `face.edges`.
@inline function _ses_drop_edge_from_face!(face::SESFace, eid::Int)
    i = findfirst(==(eid), face.edges)
    i !== nothing && deleteat!(face.edges, i)
    nothing
end

# Port of BALL's `SESSingularityCleaner::treatFirstCategory` + `noCut`
# (solventExcludedSurface.C:1664-1786). Pairs up singular spheric faces
# that share the same 3-atom set; for each pair, if both faces have
# exactly 3 SES edges (the simple "noCut" case), computes the intersection
# circle of the 2 probe spheres. If the circle's centre is inside the
# molecular atom-triangle, adds a full-circle singular SES edge bounding
# the 2 spheric faces (no endpoint vertices). The triangulator's
# `partitionSingularEdges` then samples the circle and feeds the resulting
# TEdges as boundary candidates for both spheric face's BFT.
function _ses_treat_first_category!(ses::SolventExcludedSurface{T}) where T
    rs = ses.reduced_surface
    # Group singular spheric faces by sorted atom triple.
    by_atoms = Dict{NTuple{3,Int}, Vector{Int}}()
    for (si, sf) in enumerate(ses.spheric_faces)
        sf.rs_index == 0 && continue
        rsf = rs.faces[sf.rs_index]
        rsf.v1 == 0 && continue
        rsf.singular || continue
        atoms = (rs.vertices[rsf.v1].atom,
                 rs.vertices[rsf.v2].atom,
                 rs.vertices[rsf.v3].atom)
        key = NTuple{3,Int}(sort(collect(atoms)))
        push!(get!(by_atoms, key, Int[]), si)
    end
    pr = rs.probe_radius
    modified = false
    for (_, sf_list) in by_atoms
        length(sf_list) >= 2 || continue
        # Process pairs (BALL pops 2 at a time).
        i = 1
        while i + 1 <= length(sf_list)
            si1 = sf_list[i]; si2 = sf_list[i + 1]
            i += 2
            sf1 = ses.spheric_faces[si1]
            sf2 = ses.spheric_faces[si2]
            # BALL `treatFirstCategory` switches on face1's edge count
            # (solventExcludedSurface.C:1681-1692):
            #   3 → noCut          (handled below)
            #   5 → no-op
            #   7 → twoCuts
            #   9 → deleteSimilarFaces on the underlying RSFaces; set
            #       `modified=true` so the SES gets rebuilt from scratch.
            n_e = length(sf1.edges)
            if n_e == 5
                continue
            elseif n_e == 7
                _ses_two_cuts!(ses, si1, si2)
                continue   # twoCuts replaces the noCut path
            elseif n_e == 9
                _merge_similar_rs_faces!(rs, sf1.rs_index, sf2.rs_index)
                modified = true
                continue
            elseif n_e != 3
                continue
            end
            c1 = Vector3{T}(sf1.sphere.center...)
            c2 = Vector3{T}(sf2.sphere.center...)
            d = c2 - c1
            ℓ2 = dot(d, d)
            ℓ2 > eps(T) || continue
            ℓ = sqrt(ℓ2)
            ℓ >= 2 * pr && continue   # spheres don't intersect
            # Intersection circle: centre is along the line c1→c2 at
            # distance d_int from c1; radius is r_int.
            a = ℓ / 2
            d_int = a
            r_int_sq = pr * pr - d_int * d_int
            r_int_sq > zero(T) || continue
            r_int = sqrt(r_int_sq)
            circle_p = c1 + (d_int / ℓ) * d
            circle_n = d / ℓ
            # Plane test: is circle_p inside the 3-atom triangle?
            rsf1 = rs.faces[sf1.rs_index]
            atoms = (rs.vertices[rsf1.v1].atom,
                     rs.vertices[rsf1.v2].atom,
                     rs.vertices[rsf1.v3].atom)
            p1 = Vector3{T}(rs.atoms[atoms[1]].center...)
            p2 = Vector3{T}(rs.atoms[atoms[2]].center...)
            p3 = Vector3{T}(rs.atoms[atoms[3]].center...)
            nrm = rsf1.normal
            u = cross(nrm, p1 - p2)
            v = cross(nrm, p2 - p3)
            w = cross(nrm, p3 - p1)
            diff1 = p1 - circle_p
            diff2 = p2 - circle_p
            t1 = dot(u, diff1)
            t2 = dot(v, diff2)
            t3 = dot(w, diff1)
            inside = (t1 < zero(T) && t2 < zero(T) && t3 < zero(T)) ||
                     (t1 > zero(T) && t2 > zero(T) && t3 > zero(T))
            haskey(ENV, "SES_DEBUG_NOCUT") && println(stderr,
                "[NOCUT] si1=$si1 si2=$si2 inside=$inside t1=$(round(t1;digits=4)) t2=$(round(t2;digits=4)) t3=$(round(t3;digits=4))")
            inside || continue
            # Add full-circle singular SES edge.
            circle = Circle3{T}(circle_p, circle_n, r_int)
            push!(ses.edges, SESEdge{T}(SESEdgeType.Singular, 0, 0, circle,
                                          si1, si2))
            eid = length(ses.edges)
            push!(ses.singular_edges, eid)
            push!(sf1.edges, eid)
            push!(sf2.edges, eid)
        end
    end
    # BALL `SESSingularityCleaner::run` returns false when modified=true
    # (= deleteSimilarFaces was called), prompting `SESComputer::run` to
    # discard the SES and rebuild from the now-modified RS. Caller
    # interprets the return value the same way.
    !modified
end

# Port of BALL's `SESComputer::createFreeToricFace`
# (solventExcludedSurface.C:1525). For each free toric face (RS edge
# with f1==0 or f2==0), add 2 full-circle convex SES edges — one on
# each atom — to the toric face and to each atom's contact face. These
# edges have NULL endpoints (v1=v2=0) and a circle whose center is the
# atom's contact-circle center, normal pointing toward the other atom,
# and radius from the atom's contact circle.
function _create_free_toric_edges!(ses::SolventExcludedSurface{T}) where T
    rs = ses.reduced_surface
    atom_to_cf = Dict{Int, Int}()
    for (vi, rv) in enumerate(rs.vertices)
        atom_to_cf[rv.atom] = vi
    end
    for (toric_idx, tf) in enumerate(ses.toric_faces)
        tf.type == SESFaceType.Toric || continue
        ei = tf.rs_index
        ei == 0 && continue
        re = rs.edges[ei]
        (re.f1 != 0 && re.f2 != 0) && continue   # not free
        re.v1 == 0 && continue                    # deleted edge
        !isempty(tf.edges) && continue            # already has edges
        atom_a = rs.vertices[re.v1].atom
        atom_b = rs.vertices[re.v2].atom
        # Normal direction: from contact_circle1.p toward contact_circle2.p
        # (BALL: `edge->circle_.set(circle1.p, circle1.p - circle2.p, ...)`).
        normal = re.contact_circle1.p - re.contact_circle2.p
        ℓ = sqrt(dot(normal, normal))
        normal_n = ℓ > eps(T) ? normal / ℓ : Vector3{T}(0, 0, 1)
        cf_a = get(atom_to_cf, atom_a, 0)
        cf_b = get(atom_to_cf, atom_b, 0)
        # Edge on atom_a: circle from contact_circle1
        circle_a = Circle3{T}(re.contact_circle1.p, normal_n, re.contact_circle1.r)
        push!(ses.edges, SESEdge{T}(SESEdgeType.Convex, 0, 0, circle_a, cf_a, ei))
        eid_a = length(ses.edges)
        push!(tf.edges, eid_a)
        cf_a != 0 && push!(ses.contact_faces[cf_a].edges, eid_a)
        # Edge on atom_b: circle from contact_circle2, normal flipped
        circle_b = Circle3{T}(re.contact_circle2.p, -normal_n, re.contact_circle2.r)
        push!(ses.edges, SESEdge{T}(SESEdgeType.Convex, 0, 0, circle_b, cf_b, ei))
        eid_b = length(ses.edges)
        push!(tf.edges, eid_b)
        cf_b != 0 && push!(ses.contact_faces[cf_b].edges, eid_b)
    end
    ses
end

# Reorder each SES face's `edges`/`vertices` into cyclic boundary
# order (each consecutive edge shares a vertex with the next; the
# vertices list is the corresponding cycle of corner SES vertex
# indices). Mirrors BALL's `SESFace::normalize` semantics. For faces
# whose boundary edges don't form a single cycle (singular toric
# after split — 2 disconnected triangles), each connected component
# is written contiguously; the caller can detect the boundary
# discontinuity by following edge↔vertex adjacency.
function _normalize_all_ses_faces!(ses::SolventExcludedSurface{T}) where T
    # BALL's `SESFace::normalize` (SESFace.C:152) only normalizes toric
    # faces (regular + singular); contact and spheric faces are left in
    # their as-constructed edge order. We follow the same convention so
    # the seed-edge selection in `ball_build_first_triangle!` matches
    # BALL's choice.
    for f in ses.toric_faces;   _normalize_face!(ses, f); end
end

# Helper: replace every occurrence of `old_vidx` in SES face vertex/edge
# lists and SES edge endpoints with `new_vidx`. Used by the small-toric
# cleanup to "merge" endpoint vertices that collapse together when a
# small toric face is deleted (BALL: `SESVertex::substitute`).
function _ses_substitute_vertex!(ses::SolventExcludedSurface{T},
                                   old_vidx::Int, new_vidx::Int) where T
    old_vidx == new_vidx && return
    for cf in ses.contact_faces
        for i in eachindex(cf.vertices)
            cf.vertices[i] == old_vidx && (cf.vertices[i] = new_vidx)
        end
    end
    for tf in ses.toric_faces
        for i in eachindex(tf.vertices)
            tf.vertices[i] == old_vidx && (tf.vertices[i] = new_vidx)
        end
    end
    for sf in ses.spheric_faces
        for i in eachindex(sf.vertices)
            sf.vertices[i] == old_vidx && (sf.vertices[i] = new_vidx)
        end
    end
    for e in ses.edges
        e.v1 == old_vidx && (e.v1 = new_vidx)
        e.v2 == old_vidx && (e.v2 = new_vidx)
    end
end

# Remove one specific occurrence of `vidx` from `vlist` (first found).
@inline function _remove_first!(vlist::Vector{Int}, vidx::Int)
    idx = findfirst(==(vidx), vlist)
    idx !== nothing && deleteat!(vlist, idx)
end

# Port of BALL's `cleanSingularToricFace` + `deleteSmallSingularToricFace`
# (solventExcludedSurface.C:193-246, 525-626). A singular toric face has
# 6 edges and 6 vertices forming TWO probe-triangles linked by 2 concave
# bridge edges. Layout in BALL-normalized order:
#   edge[0]=convex,  edge[1]=concave, edge[2]=convex (triangle 1)
#   edge[3]=convex,  edge[4]=concave, edge[5]=convex (triangle 2)
#   vertices p0,p1,p2 (tri1), p3,p4,p5 (tri2)
# Returns `true` if the face was deleted.
function _ses_clean_singular_toric_face!(ses::SolventExcludedSurface{T},
                                         toric_idx::Int,
                                         sqrt_density::T) where T
    tf = ses.toric_faces[toric_idx]
    tf.type == SESFaceType.ToricSingular || return false
    tf.rs_index == 0 && return false
    length(tf.edges)    == 6 || return false
    length(tf.vertices) == 6 || return false
    rs = ses.reduced_surface
    re = rs.edges[tf.rs_index]
    # BALL-normalized layout (matches `_normalize_singular_toric_face!`,
    # ported from SESFace.C:250 + findTriangle_):
    #   triangle 1: edges[1]=convex (atom arc), [2]=concave (atom_a→apex),
    #               [3]=concave (apex→atom_b)
    #               vertices[1]=atom_a, [2]=apex, [3]=atom_b
    #   triangle 2: edges[4]=convex, [5]=concave, [6]=concave
    #               vertices[4]=atom_a', [5]=apex', [6]=atom_b'
    # BALL's `cleanSingularToricFace` (solventExcludedSurface.C:193)
    # checks `v0 == v2` (triangle 1's atom arc collapses) or `v3 == v5`
    # (triangle 2's atom arc collapses).
    e0_id, e1_id, e2_id = tf.edges[1], tf.edges[2], tf.edges[3]
    e3_id, e4_id, e5_id = tf.edges[4], tf.edges[5], tf.edges[6]
    v0 = tf.vertices[1];                      v2 = tf.vertices[3]
    v3 = tf.vertices[4];                      v5 = tf.vertices[6]
    e0 = ses.edges[e0_id]; e3 = ses.edges[e3_id]
    # Sanity (mirror BALL's assumption that the layout is correct).
    e0.type == SESEdgeType.Convex || return false
    e3.type == SESEdgeType.Convex || return false
    coll1 = (v0 == v2)
    coll2 = (v3 == v5)
    should_delete = false
    if coll1
        if re.angle < T(π)
            should_delete = true
        else
            re.angle = 2 * T(π)
            return false
        end
    elseif coll2
        if re.angle < T(π)
            should_delete = true
        else
            re.angle = 2 * T(π)
            return false
        end
    else
        # BALL: `exact_number_of_segments = rsedge.angle * edge3.radius * sqrt_density`
        exact = re.angle * e3.circle.r * sqrt_density
        should_delete = exact < T(0.1)
    end
    should_delete || return false
    if haskey(ENV, "SES_DEBUG_DEL")
        println(stderr, "[DEL_SINGT] toric idx=$toric_idx rsedge_idx=$(tf.rs_index) (angle=$(round(re.angle; digits=6)) radius=$(round(e3.circle.r; digits=6)))")
    end
    # Delete: BALL's `deleteSmallSingularToricFace` (solventExcludedSurface.C:525).
    # Keep edges [1] and [4] (= concave edge1 of each triangle) and convert
    # them to singular edges. Delete the other 4 (edges 0, 2, 3, 5).
    e_cv1a = e0_id   # = convex of triangle 1 (BALL edge[0])
    e_cv1b = e2_id   # = closing concave of triangle 1 (BALL edge[2])
    e_cv2a = e3_id   # = convex of triangle 2 (BALL edge[3])
    e_cv2b = e5_id   # = closing concave of triangle 2 (BALL edge[5])
    e_conc1 = e1_id  # = concave to-be-singular for triangle 1 (BALL edge[1])
    e_conc2 = e4_id  # = concave to-be-singular for triangle 2 (BALL edge[4])
    for e_cv in (e_cv1a, e_cv1b, e_cv2a, e_cv2b)
        # remove edge from neighbour contact face
        e = ses.edges[e_cv]
        if e.f1 > 0 && e.f1 <= length(ses.contact_faces)
            _remove_first!(ses.contact_faces[e.f1].edges, e_cv)
        end
        e.v1 = 0; e.v2 = 0
    end
    # BALL's `deleteSmallSingularToricFace` (solventExcludedSurface.C:552-563)
    # ALWAYS merges p[2]→p[0] and p[5]→p[3], regardless of whether this is
    # the collapse case or the small-size-delete case. Each merge replaces
    # one of the triangle's atom-contact vertices with the other; the
    # adjacent contact face has both in its boundary so the merge closes
    # the boundary chain into a pentagon. Without this, the contact face's
    # edge list ends up missing the closing edge but keeps two distinct
    # endpoint vertices → BFT walks an open chain instead of a pentagon →
    # 32 tris emitted vs BALL's 3.
    # In my port: v0 = vertices[1] = p[0], v2 = vertices[3] = p[2],
    #             v3 = vertices[4] = p[3], v5 = vertices[6] = p[5].
    if v0 != v2
        # Remove v2 from neighbour0's vertex list (= contact face on the
        # OTHER side of edge[0]); then substitute v2 → v0 everywhere.
        n0 = e0.f1   # convex edge bounds a contact face
        if n0 > 0 && n0 <= length(ses.contact_faces)
            _remove_first!(ses.contact_faces[n0].vertices, v2)
        end
        _ses_substitute_vertex!(ses, v2, v0)
    end
    if v3 != v5
        n3 = e3.f1   # convex edge bounds a contact face
        if n3 > 0 && n3 <= length(ses.contact_faces)
            _remove_first!(ses.contact_faces[n3].vertices, v5)
        end
        _ses_substitute_vertex!(ses, v5, v3)
    end
    # Convert the concave edges to singular and add to singular list.
    for e_conc in (e_conc1, e_conc2)
        e = ses.edges[e_conc]
        e.type = SESEdgeType.Singular
        # Find the OTHER spheric face = the one across the probe axis.
        # In BALL, this is `neighbour2 = edge[2]->other(face)`. In my
        # data model, both concave edges' `f1` is the spheric face on
        # the THIS-side of the singular toric face. The other spheric
        # face is reachable via the convex edge's f1? No — convex edges
        # bound contact faces, not spheric.
        # Pragmatic: the singular toric face links TWO spheric faces
        # (the two probes' interior caps). Their indices are
        # rs.faces[re.f1] and rs.faces[re.f2] → spheric face indices
        # equal those (since ses.spheric_faces is index-aligned with
        # rs.faces). The original e.f1 == one of them; the OTHER is f2.
        if e.f1 == re.f1
            e.f2 = re.f2
        elseif e.f1 == re.f2
            e.f2 = re.f1
        else
            # fallback: pick the unused one
            e.f2 = (e.f1 == re.f1) ? re.f2 : re.f1
        end
        if e.v1 != 0 && e.v2 != 0
            v0p = ses.vertices[e.v1].point
            v1p = ses.vertices[e.v2].point
            phi = oriented_angle(v0p - e.circle.p, v1p - e.circle.p, e.circle.n)
            if phi > T(π)
                e.circle = Circle3{T}(e.circle.p, -e.circle.n, e.circle.r)
            end
        end
        push!(ses.singular_edges, e_conc)
    end
    # BALL `deleteSmallSingularToricFace` (solventExcludedSurface.C:571-572):
    # the two CONCAVE edges that bordered the OTHER spheric face have just
    # been nulled (e_cv1b = e2 and e_cv2b = e5 above). In their slots in
    # that spheric face's edge list, BALL substitutes the kept-and-now-
    # SINGULAR edges (e_conc1 = e1 and e_conc2 = e4). Without this
    # substitution, the spheric face's boundary keeps two null-endpoint
    # edges and the BFT walks them twice → duplicate triangles.
    other_sph_idx = (e_conc1 != 0) ? begin
        ec1 = ses.edges[e_conc1]
        # ses.spheric_faces is index-aligned with rs.faces; ec1.f2 was just
        # set to re.f1 or re.f2 (the other probe). That index targets the
        # spheric face whose boundary list needs the substitution.
        ec1.f2 > 0 ? ec1.f2 : 0
    end : 0
    if other_sph_idx > 0 && other_sph_idx <= length(ses.spheric_faces)
        sph_other = ses.spheric_faces[other_sph_idx]
        idx = findfirst(==(e_cv1b), sph_other.edges)
        idx !== nothing && (sph_other.edges[idx] = e_conc1)
        idx = findfirst(==(e_cv2b), sph_other.edges)
        idx !== nothing && (sph_other.edges[idx] = e_conc2)
    end
    # Mark the singular toric face as deleted.
    ses.toric_faces[toric_idx] = SESFace{T}(SESFaceType.ToricSingular,
                                            Sphere{T}(Vector3{T}(0,0,0), zero(T)),
                                            0, Int[], Int[])
    return true
end

# Port of BALL's `SolventExcludedSurface::clean()` small-toric pass +
# `deleteSmallToricFace` (solventExcludedSurface.C:152, 441-522). For
# each non-singular, non-free toric face whose
# `re.angle * edge3.circle.r * √density < 0.1` (= negligibly small
# saddle), delete the face and convert one of its concave edges into a
# SINGULAR SES edge bounding the 2 adjacent spheric faces directly.
# This removes the spurious tiny triangles that the regular toric
# triangulator would otherwise produce in dense clusters (e.g. on
# BPTI's sulfur surroundings).
function _ses_clean_small_toric_faces!(ses::SolventExcludedSurface{T};
                                        density::T = one(T)) where T
    rs = ses.reduced_surface
    # BALL: `TriangulatedSES::create` calls `ses_->clean(density)` with the
    # USER density; the small-toric criterion `angle*radius*sqrt(density) <
    # 0.1` is correspondingly tighter for higher densities (BPTI at d=6.5
    # → sqrt_density ≈ 2.55, so 2.55× more toric faces are kept than the
    # density=1 default).
    sqrt_density = sqrt(density)
    n_deleted = 0
    # BALL's `clean()` outer loop: keep sweeping until no face is deleted
    # in a full pass. A deletion may substitute vertices, which can make
    # subsequent faces' v0==v3 / v1==v2 conditions true and trigger
    # further deletions.
    done = false
    while !done
        done = true
        for toric_idx in eachindex(ses.toric_faces)
        tf = ses.toric_faces[toric_idx]
        # Handle SINGULAR toric face cleanup
        if tf.type == SESFaceType.ToricSingular
            if _ses_clean_singular_toric_face!(ses, toric_idx, sqrt_density)
                n_deleted += 1
                done = false
            end
            continue
        end
        tf.type == SESFaceType.Toric || continue
        tf.rs_index == 0 && continue
        re = rs.edges[tf.rs_index]
        (re.f1 == 0 || re.f2 == 0) && continue
        length(tf.edges) == 4 || continue
        length(tf.vertices) == 4 || continue
        # BALL's `cleanToricFace` (solventExcludedSurface.C:249-303):
        # delete OR mark-as-free the toric face if it has collapsed
        # corners (v0==v3 or v1==v2) or is too small.
        v0 = tf.vertices[1]; v1 = tf.vertices[2]
        v2 = tf.vertices[3]; v3 = tf.vertices[4]
        edge1 = ses.edges[tf.edges[2]]
        edge3 = ses.edges[tf.edges[4]]
        should_delete = false
        if v0 == v3
            if re.angle < T(π)
                should_delete = true
            else
                re.angle = 2 * T(π)   # mark as free; skip deletion
                continue
            end
        elseif v1 == v2
            if re.angle < T(π)
                should_delete = true
            else
                re.angle = 2 * T(π)
                continue
            end
        else
            exact = re.angle * edge3.circle.r * sqrt_density
            should_delete = exact < T(0.1)
        end
        should_delete || continue
        if haskey(ENV, "SES_DEBUG_SMALLT")
            println(stderr, "[DEL_SMALL_TORIC] tf_idx=$(toric_idx-1) rsedge_idx=$(tf.rs_index-1) v0=$(v0-1) v1=$(v1-1) v2=$(v2-1) v3=$(v3-1) angle=$(round(re.angle; digits=6)) radius=$(round(ses.edges[tf.edges[4]].circle.r; digits=6))")
        end
        # Identify edges in BALL's normalized order:
        #   [0] = first concave, [1] = first convex (= contact-side edge,
        #   bordering neighbour1 = contact face),
        #   [2] = second concave (= bordering neighbour2 = spheric face),
        #   [3] = second convex (= bordering neighbour3 = contact face).
        eid_0, eid_1, eid_2, eid_3 = tf.edges[1], tf.edges[2], tf.edges[3], tf.edges[4]
        e0 = ses.edges[eid_0]
        e1 = ses.edges[eid_1]
        e2 = ses.edges[eid_2]
        e3 = ses.edges[eid_3]
        # Sanity: e0/e2 should be Concave (= my f1 = RSFace index = spheric);
        # e1/e3 should be Convex (= my f1 = RSVertex index = contact face).
        e0.type == SESEdgeType.Concave || continue
        e2.type == SESEdgeType.Concave || continue
        e1.type == SESEdgeType.Convex  || continue
        e3.type == SESEdgeType.Convex  || continue
        # Neighbour indices (= the OTHER face on each edge, opposite from
        # this toric face). For concave edges, that's the spheric face
        # (= f1). For convex edges, the contact face (= f1).
        neighbour1_idx = e1.f1   # contact face for atom on edge[1]
        neighbour2_idx = e2.f1   # spheric face on probe for edge[2]
        neighbour3_idx = e3.f1   # contact face for atom on edge[3]
        # Merge p[0] with p[3]: replace p[3] with p[0] everywhere except
        # neighbour3 (where we just REMOVE p[3] from the vertex list —
        # mirrors BALL's `neighbour3->vertex_.remove(p[3]); p[3]->substitute(p[0])`).
        # BALL `deleteSmallToricFace` (solventExcludedSurface.C:465-476)
        # ALWAYS substitutes when v0 != v3 (or v1 != v2); the neighbour-side
        # vertex-list remove is just a list maintenance step. Don't gate the
        # substitute on neighbour_idx > 0 — that would leave an orphan
        # vertex referenced by other edges/faces and break boundary closure.
        if v0 != v3
            if neighbour3_idx > 0
                _remove_first!(ses.contact_faces[neighbour3_idx].vertices, v3)
            end
            _ses_substitute_vertex!(ses, v3, v0)
        end
        if v1 != v2
            if neighbour1_idx > 0
                _remove_first!(ses.contact_faces[neighbour1_idx].vertices, v2)
            end
            _ses_substitute_vertex!(ses, v2, v1)
        end
        # BALL `deleteSmallToricFace` (solventExcludedSurface.C:480-503):
        # when v[2] == v[1] (already-collapsed degenerate corner), null
        # the adjacent contact face entirely instead of just removing the
        # toric-side edge. Same for v[3] == v[0]. Without this, the
        # degenerate contact face survives into triangulation and emits
        # spurious icosphere-based triangles for a zero-area patch.
        if v1 == v2 && neighbour1_idx > 0
            ses.contact_faces[neighbour1_idx] = SESFace{T}(SESFaceType.Contact,
                                                            Sphere{T}(Vector3{T}(0,0,0), zero(T)),
                                                            0, Int[], Int[])
        elseif neighbour1_idx > 0
            _remove_first!(ses.contact_faces[neighbour1_idx].edges, eid_1)
        end
        if v0 == v3 && neighbour3_idx > 0
            ses.contact_faces[neighbour3_idx] = SESFace{T}(SESFaceType.Contact,
                                                            Sphere{T}(Vector3{T}(0,0,0), zero(T)),
                                                            0, Int[], Int[])
        elseif neighbour3_idx > 0
            _remove_first!(ses.contact_faces[neighbour3_idx].edges, eid_3)
        end
        # Substitute edge[2] -> edge[0] in neighbour2.edges (BALL:
        # `neighbour2->substitute(edge[2], edge[0])`).
        if neighbour2_idx > 0
            n2_edges = ses.spheric_faces[neighbour2_idx].edges
            idx = findfirst(==(eid_2), n2_edges)
            if idx !== nothing
                n2_edges[idx] = eid_0
            end
        end
        # Mark edges 1/2/3 as deleted: set their endpoints to 0 so any
        # downstream code that walks them treats as deleted.
        e1.v1 = 0; e1.v2 = 0
        e2.v1 = 0; e2.v2 = 0
        e3.v1 = 0; e3.v2 = 0
        # Convert e0 to a SINGULAR edge bounding the 2 spheric faces.
        # In my port the concave edge's f1 was the spheric face on e0's
        # side; the new f2 becomes neighbour2 (= the OTHER spheric face).
        e0.type = SESEdgeType.Singular
        e0.f2 = neighbour2_idx
        # BALL: if oriented angle from v0 to v1 around circle.n > π,
        # negate the normal so sampling takes the short way.
        if e0.v1 != 0 && e0.v2 != 0
            v0p = ses.vertices[e0.v1].point
            v1p = ses.vertices[e0.v2].point
            phi = oriented_angle(v0p - e0.circle.p, v1p - e0.circle.p,
                                  e0.circle.n)
            if phi > T(π)
                e0.circle = Circle3{T}(e0.circle.p, -e0.circle.n, e0.circle.r)
            end
        end
        push!(ses.singular_edges, eid_0)
        # Mark the toric face as deleted by replacing with sentinel.
        ses.toric_faces[toric_idx] = SESFace{T}(SESFaceType.Toric,
                                                  Sphere{T}(Vector3{T}(0,0,0), zero(T)),
                                                  0, Int[], Int[])
        n_deleted += 1
        done = false
        end
    end
    ses
end

function _normalize_face!(ses::SolventExcludedSurface{T}, f::SESFace{T}) where T
    n = length(f.edges)
    n <= 1 && return
    if f.type == SESFaceType.ToricSingular
        _normalize_singular_toric_face!(ses, f)
    else
        _normalize_nonsingular_toric_face!(ses, f)
    end
end

# Port of BALL `SESFace::normalizeNonSingularToricFace_` (SESFace.C:175).
# Picks the FIRST concave edge in the existing order as edge[0], then the
# SECOND concave edge as edge[2], and infers the two convex edges as the
# ones connecting them. The downstream triangulator reads edge[3] (= the
# closing convex edge in this layout) to extract the rolling-axis normal.
function _normalize_nonsingular_toric_face!(ses::SolventExcludedSurface{T}, f::SESFace{T}) where T
    edges_in = copy(f.edges)
    # Build vertex → list of incident edges (within this face).
    by_vert = Dict{Int, Vector{Int}}()
    for eid in edges_in
        e = ses.edges[eid]
        push!(get!(by_vert, e.v1, Int[]), eid)
        push!(get!(by_vert, e.v2, Int[]), eid)
    end
    remaining = Set{Int}(edges_in)
    ordered_edges = Int[]
    ordered_vertices = Int[]
    while !isempty(remaining)
        # First concave edge in the EXISTING order — match BALL's
        # `while ((*e)->type_ != TYPE_CONCAVE) e++;`.
        seed = 0
        for eid in edges_in
            if eid in remaining && ses.edges[eid].type == SESEdgeType.Concave
                seed = eid
                break
            end
        end
        seed == 0 && (seed = first(remaining))
        e0 = ses.edges[seed]
        cur_eid = seed
        next_v = e0.v2
        push!(ordered_edges, cur_eid)
        push!(ordered_vertices, e0.v1)
        delete!(remaining, cur_eid)
        guard = 0
        max_iter = length(edges_in) + 2
        while guard < max_iter
            guard += 1
            cand = 0
            for eid in get(by_vert, next_v, Int[])
                eid != cur_eid && eid in remaining && (cand = eid; break)
            end
            cand == 0 && break
            push!(ordered_edges, cand)
            push!(ordered_vertices, next_v)
            delete!(remaining, cand)
            e_next = ses.edges[cand]
            cur_eid = cand
            next_v = e_next.v1 == next_v ? e_next.v2 : e_next.v1
        end
    end
    empty!(f.edges); append!(f.edges, ordered_edges)
    empty!(f.vertices); append!(f.vertices, ordered_vertices)
end

# Port of BALL `SESFace::normalizeSingularToricFace_` (SESFace.C:250).
# Per BALL: each triangle of the 6-edge face is normalized by starting
# from a CONVEX edge (`findTriangle_`); triangle 1 starts at the FIRST
# convex edge in the existing edge list, triangle 2 starts at the LAST
# convex edge in the existing edge list. The two intermediate slots
# (edges 1 and 2 of each triangle) are filled by walking the triangle.
function _normalize_singular_toric_face!(ses::SolventExcludedSurface{T}, f::SESFace{T}) where T
    edges_in = copy(f.edges)
    by_vert = Dict{Int, Vector{Int}}()
    for eid in edges_in
        e = ses.edges[eid]
        push!(get!(by_vert, e.v1, Int[]), eid)
        push!(get!(by_vert, e.v2, Int[]), eid)
    end
    # findTriangle_ in BALL: starting from a convex `edge0` with vertices
    # (v0, v2), find edge1 = the OTHER edge incident to v0 (in this face),
    # then edge2 = the closing edge between v1 (= edge1's other endpoint)
    # and v2. Returns (edge0, edge1, edge2, v0, v1, v2).
    function _find_triangle(seed_convex_eid::Int, available::Set{Int})
        e0 = ses.edges[seed_convex_eid]
        v0 = e0.v1; v2 = e0.v2
        # edge1 = OTHER edge incident to v0 within `available`
        edge1 = 0
        for eid in get(by_vert, v0, Int[])
            (eid in available && eid != seed_convex_eid) && (edge1 = eid; break)
        end
        edge1 == 0 && return nothing
        e1 = ses.edges[edge1]
        v1 = e1.v1 == v0 ? e1.v2 : e1.v1
        # edge2 = closing edge between v1 and v2 within `available`
        edge2 = 0
        for eid in get(by_vert, v1, Int[])
            (eid == seed_convex_eid || eid == edge1) && continue
            eid in available || continue
            e2 = ses.edges[eid]
            if (e2.v1 == v1 && e2.v2 == v2) || (e2.v1 == v2 && e2.v2 == v1)
                edge2 = eid; break
            end
        end
        edge2 == 0 && return nothing
        return (seed_convex_eid, edge1, edge2, v0, v1, v2)
    end
    # Find FIRST convex edge (BALL's `findTriangle_(true)`)
    available = Set{Int}(edges_in)
    first_convex_eid = 0
    for eid in edges_in
        if ses.edges[eid].type == SESEdgeType.Convex
            first_convex_eid = eid; break
        end
    end
    first_convex_eid == 0 && return  # not a valid singular toric layout
    tri1 = _find_triangle(first_convex_eid, available)
    tri1 === nothing && return
    e0, e1, e2, p0, p1, p2 = tri1
    delete!(available, e0); delete!(available, e1); delete!(available, e2)
    # Find LAST convex edge (BALL's `findTriangle_(false)`)
    last_convex_eid = 0
    for eid in Iterators.reverse(edges_in)
        if eid in available && ses.edges[eid].type == SESEdgeType.Convex
            last_convex_eid = eid; break
        end
    end
    last_convex_eid == 0 && return
    tri2 = _find_triangle(last_convex_eid, available)
    tri2 === nothing && return
    e3, e4, e5, p3, p4, p5 = tri2
    # BALL's normalize swaps edge4 ↔ edge5 (and p5 ↔ p3) if
    # `edge1->circle_ != edge4->circle_`. This pairs the concave-edge of
    # cycle 1 with the matching concave-edge of cycle 2 (same probe
    # carrier circle).
    if ses.edges[e1].circle.p != ses.edges[e4].circle.p ||
       ses.edges[e1].circle.n != ses.edges[e4].circle.n
        e4, e5 = e5, e4
        p3, p5 = p5, p3
    end
    empty!(f.edges);    append!(f.edges,    Int[e0, e1, e2, e3, e4, e5])
    empty!(f.vertices); append!(f.vertices, Int[p0, p1, p2, p3, p4, p5])
end

# Port of BALL's `SESComputer::treatSingularToricFace`
# (solventExcludedSurface.C:1373). For each singular toric face, builds
# the singular SES edge that connects the two axis-intersection
# vertices already present on the face, and registers it as a boundary
# edge of the two adjacent contact faces. The toric face's own edge
# list is *not* augmented with the singular edge — per BALL's topology,
# the singular edge bounds the two contact faces (not the toric face).
# Port of BALL's `SESComputer::treatSingularToricFace`
# (solventExcludedSurface.C:1373). For each toric face whose underlying
# RSEdge is singular (= probe touches a 4th atom at the rolling-axis
# intersection), this routine:
#
#   1. Creates 2 NEW SES vertices at the rolling-axis intersection
#      points of the two probes (`e.intersection1` / `e.intersection2`),
#      assigning each one of the two atoms.
#   2. Creates 1 NEW *singular* SES edge connecting them.
#   3. Splits each of the 2 CONCAVE edges of the toric face into 2 halves
#      meeting at one of the new vertices. (Convex edges are NOT split.)
#   4. Updates the toric face to its 6-vertex/6-edge form (2 disconnected
#      triangular cycles, one per atom). Each cycle has 2 concave halves
#      + 1 convex edge.
#   5. Adds the singular edge + the 2 new concave halves + the 2 new
#      vertices to BOTH SPHERIC faces (= the probe faces) that share the
#      rolling axis. Contact faces are NOT modified.
#   6. Sets `tf.type = ToricSingular`.
#
# After this runs, the toric face is in BALL's normalized singular form
# and `ball_triangulate_singular_toric!` (ses_triangulator.jl) can build
# the two cycle patches. `partitionSingularEdge` (also ses_triangulator.jl)
# then samples the singular edge so the spheric-face mesh can weld with
# the singular-toric mesh at the shared boundary.
function _treat_singular_toric_faces!(ses::SolventExcludedSurface{T},
                                         vertex_grid::Dict{NTuple{3,Int}, Vector{Int}}
                                            = Dict{NTuple{3,Int}, Vector{Int}}()) where T
    rs = ses.reduced_surface
    # Seed the vertex grid with any pre-existing SES vertices so dedup
    # works on the first call too (BALL `SESComputer::run` shares
    # `vertex_grid_` across all singular-vertex creation paths).
    if isempty(vertex_grid)
        for (vidx, v) in enumerate(ses.vertices)
            _insert_vertex_grid!(vertex_grid, v.point, vidx)
        end
    end
    for (toric_idx, tf) in enumerate(ses.toric_faces)
        tf.type == SESFaceType.ToricSingular || continue
        ei = tf.rs_index
        ei == 0 && continue   # sentinel for deleted RS edge
        e = rs.edges[ei]
        e.v1 == 0 && continue # sentinel for deleted edge
        (e.f1 == 0 || e.f2 == 0) && continue
        # Build the singular-intersection circle. BALL constructs it via
        # `GetIntersection(probe1, probe2, intersection_circle)`
        # (solventExcludedSurface.C:1399) where probe1 = neighbour0
        # (the spheric face on edge[0] side = f1), probe2 = neighbour2
        # (= f2). `GetIntersection` returns normal = (probe2.p - probe1.p)/dist
        # = from f1 toward f2. We match that sign here.
        pc1 = rs.faces[e.f1].center
        pc2 = rs.faces[e.f2].center
        midp = (pc1 + pc2) / T(2)
        axis = pc2 - pc1   # BALL convention: probe2 - probe1 = f2 - f1
        ℓ = sqrt(dot(axis, axis))
        axis_n = ℓ > eps(T) ? axis / ℓ : Vector3{T}(0, 0, 1)
        half_d2 = (ℓ / 2)^2
        r2 = rs.probe_radius^2 - half_d2
        r2 < T(1e-10) && continue
        singular_circle = Circle3{T}(midp, axis_n, sqrt(r2))
        # Identify the 2 concave edges of this toric face and which RSFace
        # each lives on (= which probe).
        concave_on_f1 = 0   # concave SESEdge index with f1 == e.f1
        concave_on_f2 = 0   # concave SESEdge index with f1 == e.f2
        convex_eids = Int[]
        for se_idx in tf.edges
            se = ses.edges[se_idx]
            if se.type == SESEdgeType.Concave
                if se.f1 == e.f1
                    concave_on_f1 = se_idx
                elseif se.f1 == e.f2
                    concave_on_f2 = se_idx
                end
            elseif se.type == SESEdgeType.Convex
                push!(convex_eids, se_idx)
            end
        end
        (concave_on_f1 == 0 || concave_on_f2 == 0) && continue
        length(convex_eids) == 2 || continue
        # Atom identity of each concave-edge endpoint (a concave edge on
        # probe X connects atom_A's corner on probe X to atom_B's corner
        # on probe X). We use rs.vertices[e.v1].atom and e.v2.atom for the
        # 2 atoms participating in the rolling axis.
        atom_a = rs.vertices[e.v1].atom
        atom_b = rs.vertices[e.v2].atom
        # Create the 2 new SES vertices, one per axis intersection. BALL
        # assigns each to one of the two atoms (`face->rsedge_->getVertex(ip)`).
        # We follow the same convention: intersection1 → atom of e.v1 (= atom_a);
        # intersection2 → atom of e.v2 (= atom_b).
        # BALL `createSingularVertex` (solventExcludedSurface.C:1458):
        # check vertex_grid for an existing SES vertex at this point
        # (tolerance 1e-3 Å); reuse if found, allocate otherwise.
        ip1 = _ses_create_or_find_singular_vertex!(ses, vertex_grid,
                                                    e.intersection1, atom_a, 0, ei)
        ip2 = _ses_create_or_find_singular_vertex!(ses, vertex_grid,
                                                    e.intersection2, atom_b, 0, ei)
        # Create the singular edge between the 2 new vertices.
        push!(ses.edges, SESEdge{T}(SESEdgeType.Singular, ip1, ip2,
                                     singular_circle, e.f1, e.f2))
        singular_eid = length(ses.edges)
        push!(ses.singular_edges, singular_eid)
        # BALL's post-construction normal-flip check
        # (solventExcludedSurface.C:1427-1434): compute oriented angle phi
        # from ip1 to ip2 around the circle normal; if (rsedge.angle - π) and
        # (phi - π) have opposite signs, negate the singular circle's normal
        # so the rotation direction is consistent with the toric face's
        # rolling angle convention.
        let se = ses.edges[singular_eid]
            phi = oriented_angle(e.intersection1 - se.circle.p,
                                  e.intersection2 - se.circle.p,
                                  se.circle.n)
            if (e.angle - T(π)) * (phi - T(π)) < zero(T)
                se.circle = Circle3{T}(se.circle.p, -se.circle.n, se.circle.r)
            end
        end
        # Split the concave edges. For each, find which endpoint is on
        # atom_a vs atom_b, and put one half on each cycle.
        # Concave-on-f1: connects atom_a's corner on probe-f1 to atom_b's
        # corner on probe-f1. After split:
        #   * Original edge → atom_a-corner ↔ ip1 (atom_a cycle)
        #   * New half     → atom_b-corner ↔ ip2 (atom_b cycle)
        # Same for concave-on-f2.
        new_halves = Tuple{Int,Int}[]   # (new_eid, atom_label) where atom_label = 1 for atom_a, 2 for atom_b
        function split_concave!(ce_idx::Int, probe_f::Int)
            ce = ses.edges[ce_idx]
            v1_sv = ses.vertices[ce.v1]
            v2_sv = ses.vertices[ce.v2]
            corner_a = v1_sv.atom == atom_a ? ce.v1 :
                       v2_sv.atom == atom_a ? ce.v2 : 0
            corner_b = v1_sv.atom == atom_b ? ce.v1 :
                       v2_sv.atom == atom_b ? ce.v2 : 0
            (corner_a == 0 || corner_b == 0) && return (0, 0, 0)
            # Mutate ce in place to be the atom_a half: corner_a ↔ ip1
            ce.v1 = corner_a
            ce.v2 = ip1
            # Create new concave half: ip2 ↔ corner_b, same circle, same probe.
            # BALL `treatSingularToricFace` (solventExcludedSurface.C:1479)
            # passes `(p[1], new_point3, true)` to `updateEdge`; since the
            # source clone's `vertex_[0]` was `p[0]` (= NOT `p[1]`),
            # updateEdge takes the else branch and sets
            # `vertex_[0] = new_point3`, `vertex_[1] = p[1]`. So the new
            # edge's v1 = apex (ip2), v2 = corner_b.
            # My earlier order (corner_b, ip2) gave the COMPLEMENT arc
            # angle (~6 rad instead of ~0.4 rad), inflating n_seg from 1
            # to 13 on sph 278's edge 4231 → +12 TEdges on the boundary,
            # which the BFT then over-triangulates.
            push!(ses.edges, SESEdge{T}(SESEdgeType.Concave, ip2, corner_b,
                                         ce.circle, probe_f, ei))
            new_eid = length(ses.edges)
            return (ce_idx, new_eid, corner_b)
        end
        ok = true
        (orig_f1, new_f1, corner_b_f1) = split_concave!(concave_on_f1, e.f1)
        orig_f1 == 0 && (ok = false)
        (orig_f2, new_f2, corner_b_f2) = split_concave!(concave_on_f2, e.f2)
        orig_f2 == 0 && (ok = false)
        ok || continue
        # Find the convex edges and which cycle each belongs to (= which
        # atom). A convex edge on atom_X connects two corners on atom_X.
        convex_a = 0; convex_b = 0
        for ce_idx in convex_eids
            ce = ses.edges[ce_idx]
            v1_sv = ses.vertices[ce.v1]
            if v1_sv.atom == atom_a
                convex_a = ce_idx
            elseif v1_sv.atom == atom_b
                convex_b = ce_idx
            end
        end
        (convex_a == 0 || convex_b == 0) && continue
        # Rebuild tf.edges into 6 edges arranged as 2 contiguous cycles:
        #   cycle A (atom_a): [convex_a, concave_a_half_on_f1, concave_a_half_on_f2]
        #   cycle B (atom_b): [convex_b, concave_b_half_on_f1, concave_b_half_on_f2]
        # The 2 cycles meet at ip1 (cycle A's apex) and ip2 (cycle B's apex)
        # but share no edges — they are disconnected as required by
        # ball_triangulate_singular_toric!.
        empty!(tf.edges)
        empty!(tf.vertices)
        # Cycle A
        ce_a_on_f1_atom_a = orig_f1   # mutated to corner_a_f1 ↔ ip1
        ce_a_on_f2_atom_a = orig_f2   # mutated to corner_a_f2 ↔ ip1
        push!(tf.edges, convex_a)
        push!(tf.edges, ce_a_on_f1_atom_a)
        push!(tf.edges, ce_a_on_f2_atom_a)
        ce_b_on_f1_atom_b = new_f1    # new: corner_b_f1 ↔ ip2
        ce_b_on_f2_atom_b = new_f2    # new: corner_b_f2 ↔ ip2
        push!(tf.edges, convex_b)
        push!(tf.edges, ce_b_on_f1_atom_b)
        push!(tf.edges, ce_b_on_f2_atom_b)
        # Vertex list (3 per cycle, contiguous). BALL's `SESFace::normalize`
        # for singular toric faces produces order [atom_a, apex, atom_b]
        # per cycle (= the order `cleanSingularToricFace` reads as v[0..5]).
        # The convex edge connects atom_a to atom_b on the SAME atom
        # surface; the apex is the axis-intersection vertex shared by the
        # two concave halves of that cycle. We don't run normalize after
        # `_treat_singular_toric_faces!` (because the face type changes
        # mid-construction), so emit in the normalized order directly.
        ce = ses.edges[convex_a]
        push!(tf.vertices, ce.v1)   # atom_a (= one endpoint of convex arc)
        push!(tf.vertices, ip1)     # apex (= axis-intersection vertex)
        push!(tf.vertices, ce.v2)   # atom_b (= other endpoint of convex arc)
        ce = ses.edges[convex_b]
        push!(tf.vertices, ce.v1)   # atom_a'
        push!(tf.vertices, ip2)     # apex'
        push!(tf.vertices, ce.v2)   # atom_b'
        # Update the two SPHERIC faces (probe faces) bordering this
        # rolling axis: add the singular edge + the 2 new concave halves
        # + the 2 new vertices to each spheric face's boundary.
        # (BALL: neighbour0->edge_.push_back(new_edge0); neighbour0->edge_.push_back(new_edge); etc.)
        # spheric_faces indexed by RSFace; we add to both faces e.f1 and e.f2.
        sph_f1 = ses.spheric_faces[e.f1]
        sph_f2 = ses.spheric_faces[e.f2]
        push!(sph_f1.edges, new_f1)
        push!(sph_f1.edges, singular_eid)
        push!(sph_f1.vertices, ip1)
        push!(sph_f1.vertices, ip2)
        push!(sph_f2.edges, new_f2)
        push!(sph_f2.edges, singular_eid)
        push!(sph_f2.vertices, ip1)
        push!(sph_f2.vertices, ip2)
    end
    ses
end

"""
    $(TYPEDSIGNATURES)

Compute the analytical solvent-excluded surface directly from an atom
container. Equivalent to `compute_ses(compute_reduced_surface(ac; probe_radius))`.
"""
function compute_ses(ac::AbstractAtomContainer{T};
                     probe_radius::T = T(1.5),
                     density::Real = 1.0) where T
    compute_ses(compute_reduced_surface(ac; probe_radius); density = density)
end

# ===========================================================================
# Port of BALL's `SESSingularityCleaner::treatSecondCategory` + helpers
# (solventExcludedSurface.C:1952-2616). Given the singular SES edges that
# already exist (from singular toric face processing, small-toric cleanup,
# and `treatFirstCategory`'s noCut), look for places where each singular
# edge is geometrically intersected by a NEARBY probe sphere. At each
# such intersection point, split the original singular edge into shorter
# sub-edges. This adds vertices and arcs on the spheric face boundaries
# that close stacked-probe topology in BALL's data model.
# ===========================================================================

# Port of BALL's 3-sphere intersection (analyticalGeometry.h:929). Returns
# the 2 intersection points (or `nothing` if they don't exist). All three
# spheres are assumed to have the SAME radius (BALL's call site does so).
function _three_sphere_intersection(c1::Vector3{T}, c2::Vector3{T},
                                     c3::Vector3{T}, r::T) where T
    r² = r * r
    p1² = dot(c1, c1)
    p2² = dot(c2, c2)
    p3² = dot(c3, c3)
    # Perpendicular bisector planes (equidistant from each sphere pair):
    # Plane 1: (c2 - c1) · X = (|c2|² - |c1|²) / 2
    # Plane 2: (c3 - c1) · X = (|c3|² - |c1|²) / 2
    u = (p2² - p1²) / 2
    v = (p3² - p1²) / 2
    n1 = c2 - c1; n2 = c3 - c1
    # Line of intersection: direction = n1 × n2; need a point on the line.
    dir = cross(n1, n2)
    ℓ2 = dot(dir, dir)
    ℓ2 > eps(T) || return nothing
    # Find a point P0 on both planes by solving a 3×3 system: stack
    # n1·X = u, n2·X = v, dir·X = 0 (anchor on perpendicular plane).
    A = [n1[1] n1[2] n1[3];
         n2[1] n2[2] n2[3];
         dir[1] dir[2] dir[3]]
    detA = A[1,1]*(A[2,2]*A[3,3] - A[2,3]*A[3,2]) -
           A[1,2]*(A[2,1]*A[3,3] - A[2,3]*A[3,1]) +
           A[1,3]*(A[2,1]*A[3,2] - A[2,2]*A[3,1])
    abs(detA) > eps(T) || return nothing
    bvec = (u, v, zero(T))
    # Cramer's rule
    function _det3(c1::Int)
        cols = ntuple(j -> j == c1 ? bvec : (A[1,j], A[2,j], A[3,j]), 3)
        cols[1][1] * (cols[2][2]*cols[3][3] - cols[2][3]*cols[3][2]) -
        cols[1][2] * (cols[2][1]*cols[3][3] - cols[2][3]*cols[3][1]) +
        cols[1][3] * (cols[2][1]*cols[3][2] - cols[2][2]*cols[3][1])
    end
    P0 = Vector3{T}(_det3(1)/detA, _det3(2)/detA, _det3(3)/detA)
    # Substitute X = P0 + t·dir into |X - c1|² = r²  →  quadratic in t.
    diff = c1 - P0
    a = dot(dir, dir)
    b = -2 * dot(diff, dir)
    c = dot(diff, diff) - r²
    disc = b * b - 4 * a * c
    disc >= zero(T) || return nothing
    sd = sqrt(disc)
    t1 = (-b + sd) / (2 * a)
    t2 = (-b - sd) / (2 * a)
    return P0 + t1 * dir, P0 + t2 * dir
end

# Port of `SESSingularityCleaner::probeIntersection` with memoization.
# Cache key: sorted triple of spheric face indices. Returns
# (point1, point2) or `nothing`. The cache uses `Bool` for cached misses.
function _probe_intersection!(cache::Dict{NTuple{3,Int}, Any},
                                ses::SolventExcludedSurface{T},
                                f1::Int, f2::Int, f3::Int, pr::T) where T
    key = NTuple{3,Int}(sort([f1, f2, f3]))
    haskey(cache, key) && return cache[key]
    sf1 = ses.spheric_faces[f1]
    sf2 = ses.spheric_faces[f2]
    sf3 = ses.spheric_faces[f3]
    c1 = Vector3{T}(sf1.sphere.center...)
    c2 = Vector3{T}(sf2.sphere.center...)
    c3 = Vector3{T}(sf3.sphere.center...)
    result = _three_sphere_intersection(c1, c2, c3, pr)
    cache[key] = result
    return result
end

# Port of `SESSingularityCleaner::getIntersectionPointsAndAngles`.
# Given a circle and a reference point `p_ref` on the circle, compute
# the 3-sphere intersection of probes (f1, f2, f3) and return their
# oriented angles around the circle relative to `p_ref`, swapped so
# phi1 ≤ phi2. Returns (phi1, point1, phi2, point2) or `nothing`.
function _intersection_points_and_angles!(cache::Dict{NTuple{3,Int}, Any},
                                            circle::Circle3{T},
                                            p_ref::Vector3{T},
                                            ses::SolventExcludedSurface{T},
                                            f1::Int, f2::Int, f3::Int, pr::T) where T
    res = _probe_intersection!(cache, ses, f1, f2, f3, pr)
    res === nothing && return nothing
    p1, p2 = res
    phi1 = oriented_angle(p_ref - circle.p, p1 - circle.p, circle.n)
    phi2 = oriented_angle(p_ref - circle.p, p2 - circle.p, circle.n)
    # Wrap angles equal to 2π to 0 (BALL: with epsilon=0.001).
    abs(phi1 - 2 * T(π)) < T(0.001) && (phi1 = zero(T))
    abs(phi2 - 2 * T(π)) < T(0.001) && (phi2 = zero(T))
    if phi2 < phi1
        phi1, phi2 = phi2, phi1
        p1, p2 = p2, p1
    end
    return phi1, p1, phi2, p2
end

# Port of `SESSingularityCleaner::isIntersection`. Checks if the
# intersection segment [min_phi, max_phi] is on the singular edge's
# arc (covering [0, phi]).
function _is_intersection(min_phi::T, max_phi::T, phi::T,
                            middle::Vector3{T}, probe_p::Vector3{T},
                            probe_r::T) where T
    if max_phi > phi
        return false
    elseif !isapprox(min_phi, zero(T); atol = T(0.001)) || max_phi < phi
        return true
    else
        # max_phi == phi (within tolerance) && min_phi == 0 — degenerate.
        # Test: is `middle` INSIDE the probe (i.e., we're crossing INTO it)?
        d² = dot(middle - probe_p, middle - probe_p)
        return !(d² > probe_r * probe_r + T(1e-6))
    end
end

# Port of `SESSingularityCleaner::getIntersectionsOfSingularEdge`.
# Returns a vector of `((phi, candidate_face_idx), point)` tuples.
function _intersections_of_singular_edge(edge::SESEdge{T},
                                            phi::T,
                                            candidates::Vector{Int},
                                            ses::SolventExcludedSurface{T},
                                            cache::Dict{NTuple{3,Int}, Any},
                                            pr::T,
                                            face1::Int, face2::Int) where T
    # Midpoint of the singular arc (at phi/2).
    p0 = ses.vertices[edge.v1].point - edge.circle.p
    # Rodrigues rotation by phi/2 around edge.circle.n.
    function _rotate(v::Vector3{T}, axis::Vector3{T}, angle::T)
        c = cos(angle); s = sin(angle)
        v * c + cross(axis, v) * s + axis * dot(axis, v) * (one(T) - c)
    end
    mid_rel = _rotate(p0, edge.circle.n, phi / 2)
    middle = mid_rel + edge.circle.p

    intersections = Tuple{T, Int, Vector3{T}}[]
    p_ref = ses.vertices[edge.v1].point
    for fc in candidates
        (fc == face1 || fc == face2) && continue
        res = _intersection_points_and_angles!(cache, edge.circle, p_ref,
                                                 ses, face1, face2, fc, pr)
        res === nothing && continue
        phi1, point1, phi2, point2 = res
        sf = ses.spheric_faces[fc]
        probe_p = Vector3{T}(sf.sphere.center...)
        if _is_intersection(phi1, phi2, phi, middle, probe_p, pr)
            push!(intersections, (phi1, fc, point1))
            push!(intersections, (phi2, fc, point2))
        end
    end
    intersections
end

# Port of `SESSingularityCleaner::getExtrema`. Returns (min_list, max_list)
# where each is the set of intersections at the global min/max phi.
function _intersection_extrema(intersections::Vector{Tuple{T, Int, Vector3{T}}}) where T
    isempty(intersections) && return Tuple{T, Int, Vector3{T}}[],
                                       Tuple{T, Int, Vector3{T}}[]
    min_phi = T(2 * π); max_phi = zero(T)
    min_list = Tuple{T, Int, Vector3{T}}[]
    max_list = Tuple{T, Int, Vector3{T}}[]
    eps_ang = T(0.0001)
    for it in intersections
        phi = it[1]
        if phi <= min_phi + eps_ang
            if phi < min_phi - eps_ang
                empty!(min_list)
                min_phi = phi
            end
            push!(min_list, it)
        end
        if phi >= max_phi - eps_ang
            if phi > max_phi + eps_ang
                empty!(max_list)
                max_phi = phi
            end
            push!(max_list, it)
        end
    end
    return min_list, max_list
end

# Port of `SESSingularityCleaner::vertexExists`. Returns the index of an
# existing SES vertex at `point` (within 1e-3 Å), or 0 if none.
function _ses_vertex_exists(point::Vector3{T},
                              vertex_grid::Dict{NTuple{3,Int}, Vector{Int}},
                              ses::SolventExcludedSurface{T}) where T
    cell = T(0.05)
    k = (round(Int, point[1] / cell),
         round(Int, point[2] / cell),
         round(Int, point[3] / cell))
    # Check 3×3×3 cells around k.
    for dx in -1:1, dy in -1:1, dz in -1:1
        ck = (k[1] + dx, k[2] + dy, k[3] + dz)
        haskey(vertex_grid, ck) || continue
        for vidx in vertex_grid[ck]
            d = ses.vertices[vidx].point - point
            dot(d, d) < T(1e-6) && return vidx
        end
    end
    return 0
end

@inline function _insert_vertex_grid!(vertex_grid::Dict{NTuple{3,Int}, Vector{Int}},
                                       point::Vector3{T}, vidx::Int) where T
    cell = T(0.05)
    k = (round(Int, point[1] / cell),
         round(Int, point[2] / cell),
         round(Int, point[3] / cell))
    push!(get!(vertex_grid, k, Int[]), vidx)
end

# Port of BALL `SESComputer::createSingularVertex`
# (solventExcludedSurface.C:1458-1496). If a SES vertex already exists
# at `point` (within 1e-3 Å, matching BALL's temporary EPSILON), return
# its index. Otherwise create a new SESVertex and register it in the
# vertex grid. Mirrors BALL's `vertexExists`-then-allocate pattern.
function _ses_create_or_find_singular_vertex!(ses::SolventExcludedSurface{T},
                                                vertex_grid::Dict{NTuple{3,Int}, Vector{Int}},
                                                point::Vector3{T}, atom::Int,
                                                f1::Int, rs_index::Int) where T
    test = _ses_vertex_exists(point, vertex_grid, ses)
    test != 0 && return test
    push!(ses.vertices, SESVertex{T}(point, atom, f1, rs_index))
    vidx = length(ses.vertices)
    _insert_vertex_grid!(vertex_grid, point, vidx)
    vidx
end

# Builder-side variant of `_ses_create_or_find_singular_vertex!` used
# during `_compute_ses_once`, before the SolventExcludedSurface struct
# exists. Operates on `b.vertices` (which becomes `ses.vertices` after
# construction, so the indices remain valid).
function _b_create_or_find_singular_vertex!(b::_SESBuilder{T},
                                              vertex_grid::Dict{NTuple{3,Int}, Vector{Int}},
                                              point::Vector3{T}, atom::Int,
                                              f1::Int, rs_index::Int) where T
    # Inline the grid lookup against `b.vertices`.
    cell = T(0.05)
    k = (round(Int, point[1] / cell),
         round(Int, point[2] / cell),
         round(Int, point[3] / cell))
    for dx in -1:1, dy in -1:1, dz in -1:1
        ck = (k[1] + dx, k[2] + dy, k[3] + dz)
        haskey(vertex_grid, ck) || continue
        for vidx in vertex_grid[ck]
            d = b.vertices[vidx].point - point
            dot(d, d) < T(1e-6) && return vidx
        end
    end
    push!(b.vertices, SESVertex{T}(point, atom, f1, rs_index))
    vidx = length(b.vertices)
    push!(get!(vertex_grid, k, Int[]), vidx)
    vidx
end

# Port of `SESSingularityCleaner::buildEndEdge`. For the min OR max
# extremum group, find/create a vertex at the absolute extremum point
# and (if it differs from the current edge endpoint) create a new
# singular edge from that endpoint to the new vertex.
# Returns (vertex_idx, actual_extremum_face_idx, new_edge_idx_or_0).
function _ses_build_end_edge!(ses::SolventExcludedSurface{T},
                                edge_idx::Int,
                                extrema::Vector{Tuple{T, Int, Vector3{T}}},
                                vertex_grid::Dict{NTuple{3,Int}, Vector{Int}},
                                use_min::Bool,
                                skip_set::Union{Set{Int}, Nothing} = nothing) where T
    edge = ses.edges[edge_idx]
    vertex_idx = 0
    actual_extremum = 0
    # First, see if any extremum point coincides with an existing vertex.
    for it in extrema
        test = _ses_vertex_exists(it[3], vertex_grid, ses)
        if test != 0
            vertex_idx = test
            actual_extremum = it[2]
        end
    end
    if vertex_idx == 0
        # Pick the absolute extremum.
        abs_ext = extrema[1]
        for it in extrema
            if use_min
                it[1] < abs_ext[1] && (abs_ext = it)
            else
                it[1] > abs_ext[1] && (abs_ext = it)
            end
        end
        actual_extremum = abs_ext[2]
        # Create new SES vertex.
        outward = edge.circle.p - abs_ext[3]
        ℓ = sqrt(dot(outward, outward))
        normal = ℓ > eps(T) ? outward / ℓ : Vector3{T}(0, 0, 1)
        push!(ses.vertices, SESVertex{T}(abs_ext[3], 0, 0, 0))
        # SESVertex has (point, atom, rs_face, rs_edge) — use atom=0 (mark
        # this vertex as a probe-probe intersection point), rs_face/rs_edge=0.
        vertex_idx = length(ses.vertices)
        _insert_vertex_grid!(vertex_grid, abs_ext[3], vertex_idx)
    end
    # If the chosen vertex differs from the current edge endpoint, split.
    v_slot = use_min ? 1 : 2
    endpoint = use_min ? edge.v1 : edge.v2
    new_edge_idx = 0
    if vertex_idx != endpoint
        # Create a clone of `edge` with one endpoint replaced. BALL clones
        # the circle but uses NULL rsedge_; for us, set rs_index=0.
        # BALL `buildEndEdge` (solventExcludedSurface.C:2426-2428):
        #   new_edge = new SESEdge(*edge, true);   // clone
        #   new_edge->vertex_[1-v] = vertex;        // v=0 for min, 1 for max
        # So for min (v=0): vertex_[0] stays = edge.v1, vertex_[1] = new
        #    for max (v=1): vertex_[0] = new, vertex_[1] stays = edge.v2
        # My earlier ordering had this REVERSED for the min case, which
        # gave the COMPLEMENT arc angle in `partitionNonFreeSingularEdge`
        # and inflated the per-edge tedge count for the new sub-edges.
        new_v1 = use_min ? edge.v1 : vertex_idx
        new_v2 = use_min ? vertex_idx : edge.v2
        push!(ses.edges, SESEdge{T}(SESEdgeType.Singular,
                                      new_v1, new_v2, edge.circle,
                                      edge.f1, edge.f2))
        new_edge_idx = length(ses.edges)
        push!(ses.singular_edges, new_edge_idx)
        # BALL `buildEndEdge` (solventExcludedSurface.C:2433) uses
        # `singular_edges_.push_front(new_edge)`, so the new shortened
        # end-edge is NOT visited during the current `treatSecondCategory`
        # pass. My port pushes back; record the new eid in `skip_set` so
        # the outer loop skips it. Same effect, no index-shift bookkeeping.
        skip_set === nothing || push!(skip_set, new_edge_idx)
        if haskey(ENV, "SES_DEBUG_NEW")
            cp = edge.circle.p
            println(stderr, "[NEW_SING_END] eid=$new_edge_idx f1=$(edge.f1) f2=$(edge.f2) v1=$new_v1 v2=$new_v2 circle.p=($(round(cp[1];digits=4)),$(round(cp[2];digits=4)),$(round(cp[3];digits=4))) r=$(round(edge.circle.r;digits=4))")
        end
        # Register on both adjacent spheric faces.
        push!(ses.spheric_faces[edge.f1].edges, new_edge_idx)
        push!(ses.spheric_faces[edge.f2].edges, new_edge_idx)
        # Add the new vertex to both adjacent spheric faces' vertex lists.
        let v = use_min ? new_v1 : new_v2
            v in ses.spheric_faces[edge.f1].vertices || push!(ses.spheric_faces[edge.f1].vertices, v)
            v in ses.spheric_faces[edge.f2].vertices || push!(ses.spheric_faces[edge.f2].vertices, v)
        end
    end
    return vertex_idx, actual_extremum, new_edge_idx
end

# Port of `SESSingularityCleaner::buildEdge`. Builds a singular SES edge
# between spheric_face[face1] and spheric_face[face2] from the current
# vertex toward `end_face_idx`. Updates `vertex` and `face2` for chaining.
# Returns `(new_vertex_idx, new_face2_idx)` for the next iteration; if
# they signal "stop", new_vertex_idx == 0.
function _ses_build_edge!(ses::SolventExcludedSurface{T},
                            cache::Dict{NTuple{3,Int}, Any},
                            edge_idx::Int,
                            face1::Int, face2::Int, end_face::Int,
                            vertex_idx::Int,
                            indices::Set{Int},
                            vertex_grid::Dict{NTuple{3,Int}, Vector{Int}},
                            use_min::Bool) where T
    sf1 = ses.spheric_faces[face1]
    sf2 = ses.spheric_faces[face2]
    c1 = Vector3{T}(sf1.sphere.center...)
    c2 = Vector3{T}(sf2.sphere.center...)
    pr = ses.reduced_surface.probe_radius
    # Intersection circle of probes face1 and face2.
    d = c2 - c1
    ℓ = sqrt(dot(d, d))
    ℓ < 2 * pr || return 0, face2
    d_int = ℓ / 2
    r_int² = pr * pr - d_int * d_int
    r_int² > zero(T) || return 0, face2
    circle_p = c1 + (d_int / ℓ) * d
    circle_n = d / ℓ
    circle = Circle3{T}(circle_p, circle_n, sqrt(r_int²))
    # Orient circle.n per BALL's sign rule. Need the original singular edge.
    orig_edge = ses.edges[edge_idx]
    sign_factor = use_min ? -one(T) : one(T)
    s1 = dot(c1 - orig_edge.circle.p, orig_edge.circle.n)
    s2 = dot(c1 - circle.p, circle.n)
    if s1 * s2 * sign_factor > zero(T)
        circle = Circle3{T}(circle.p, -circle.n, circle.r)
    end
    point = ses.vertices[vertex_idx].point
    # Find the candidate face among `indices` whose intersection with
    # (face1, face2) has the SMALLEST positive phi from `point`.
    min_phi = T(2 * π)
    min_candidates = Tuple{Vector3{T}, Int}[]
    for fc in indices
        (fc == face1 || fc == face2) && continue
        res = _intersection_points_and_angles!(cache, circle, point,
                                                 ses, face1, face2, fc, pr)
        res === nothing && continue
        phi1, p1, phi2, p2 = res
        for (phi_x, p_x) in ((phi1, p1), (phi2, p2))
            if phi_x > T(0.001) && phi_x <= min_phi + T(0.001)
                if phi_x < min_phi - T(0.001)
                    empty!(min_candidates)
                    min_phi = phi_x
                end
                push!(min_candidates, (p_x, fc))
            end
        end
    end
    isempty(min_candidates) && return 0, face2
    # Prefer the candidate that matches `end_face`.
    new_face2 = face2
    new_point = min_candidates[1][1]
    found_end = false
    for c in min_candidates
        if c[2] == end_face
            new_point = c[1]
            new_face2 = end_face
            found_end = true
            break
        end
    end
    if !found_end
        # BALL `buildEdge` fallback (solventExcludedSurface.C:2451-2464):
        # iterate min_list from start; each iter SETS point2/face2 and
        # checks vertexExists. Exit early on found. If none found, point2
        # and face2 retain the LAST iterated values (not the first).
        # My earlier port defaulted to the FIRST candidate on no-match,
        # which differs from BALL's last-iterated fallback.
        for c in min_candidates
            new_point = c[1]
            new_face2 = c[2]
            test = _ses_vertex_exists(c[1], vertex_grid, ses)
            test != 0 && break
        end
    end
    # Bail out if the two spheric faces are already neighbors via the
    # singular edges list. We approximate "is_neighbored_to" by checking
    # if face1 and face2 already share a SES edge.
    sf1_now = ses.spheric_faces[face1]
    sf2_now = ses.spheric_faces[face2]
    already_neighbored = !isempty(intersect(Set(sf1_now.edges), Set(sf2_now.edges)))
    if already_neighbored
        return 0, new_face2
    end
    # Create or find new vertex.
    new_vertex_idx = _ses_vertex_exists(new_point, vertex_grid, ses)
    if new_vertex_idx == 0
        outward = circle.p - new_point
        ℓo = sqrt(dot(outward, outward))
        push!(ses.vertices, SESVertex{T}(new_point, 0, 0, 0))
        new_vertex_idx = length(ses.vertices)
        _insert_vertex_grid!(vertex_grid, new_point, new_vertex_idx)
    end
    # Create new singular SES edge between vertex_idx and new_vertex_idx,
    # bounding face1 and face2.
    push!(ses.edges, SESEdge{T}(SESEdgeType.Singular,
                                  vertex_idx, new_vertex_idx, circle,
                                  face1, face2))
    new_edge_idx = length(ses.edges)
    push!(ses.singular_edges, new_edge_idx)
    if haskey(ENV, "SES_DEBUG_NEW")
        cp = circle.p
        println(stderr, "[NEW_SING] eid=$new_edge_idx f1=$face1 f2=$face2 v1=$vertex_idx v2=$new_vertex_idx circle.p=($(round(cp[1];digits=4)),$(round(cp[2];digits=4)),$(round(cp[3];digits=4))) r=$(round(circle.r;digits=4))")
    end
    push!(sf1_now.edges, new_edge_idx)
    push!(sf2_now.edges, new_edge_idx)
    vertex_idx in sf1_now.vertices || push!(sf1_now.vertices, vertex_idx)
    vertex_idx in sf2_now.vertices || push!(sf2_now.vertices, vertex_idx)
    new_vertex_idx in sf1_now.vertices || push!(sf1_now.vertices, new_vertex_idx)
    new_vertex_idx in sf2_now.vertices || push!(sf2_now.vertices, new_vertex_idx)
    return new_vertex_idx, new_face2
end

# Port of `SESSingularityCleaner::treatSingularEdge`. Splits a singular
# SES edge at probe-probe intersection points and creates the resulting
# sub-edges. Returns `true` if the edge should be deleted.
function _ses_treat_singular_edge!(ses::SolventExcludedSurface{T},
                                     edge_idx::Int,
                                     candidates_per_pos::Vector{Vector{Int}},
                                     cache::Dict{NTuple{3,Int}, Any},
                                     vertex_grid::Dict{NTuple{3,Int}, Vector{Int}},
                                     pr::T,
                                     skip_set::Union{Set{Int}, Nothing} = nothing) where T
    edge = ses.edges[edge_idx]
    edge.type == SESEdgeType.Singular || return false
    edge.v1 == 0 && return false  # full-circle variant — not handled here
    # phi = oriented angle from v0 to v1 around circle.n
    rel0 = ses.vertices[edge.v1].point - edge.circle.p
    rel1 = ses.vertices[edge.v2].point - edge.circle.p
    phi = oriented_angle(rel0, rel1, edge.circle.n)
    # Get candidates near circle.p — for simplicity here, use the global
    # spheric face list (small N for BPTI). Faster lookup can be added.
    candidates = collect(1:length(ses.spheric_faces))
    intersections = _intersections_of_singular_edge(edge, phi, candidates,
                                                      ses, cache, pr,
                                                      edge.f1, edge.f2)
    haskey(ENV, "SES_DEBUG_TSE") && println(stderr, "[TSE] eid=$edge_idx n_intersections=$(length(intersections))")
    isempty(intersections) && return false
    min_list, max_list = _intersection_extrema(intersections)
    indices = Set{Int}()
    for it in min_list; push!(indices, it[2]); end
    for it in max_list; push!(indices, it[2]); end
    push!(indices, edge.f1); push!(indices, edge.f2)
    face1 = edge.f1; face2 = edge.f2
    v1_idx, actual_min, _ = _ses_build_end_edge!(ses, edge_idx, min_list,
                                                    vertex_grid, true, skip_set)
    v2_idx, actual_max, _ = _ses_build_end_edge!(ses, edge_idx, max_list,
                                                    vertex_grid, false, skip_set)
    # BALL `treatSingularEdge` (solventExcludedSurface.C:2049-2078) does
    # TWO passes — once walking from face1 toward face2 (= one side of
    # the singular edge), once with face1/face2 swapped (= other side).
    # Each pass has a fallback from actual_min to actual_max if the first
    # didn't reach face2. The second pass is what adds edges to the
    # OPPOSITE spheric face's boundary. My earlier single-pass port
    # left the second spheric face under-edged, which polluted its BST.
    function walk_one_pass(start_face::Int, target_face::Int)
        next_face = actual_min
        vertex = v1_idx
        iter_guard = 0
        while next_face != target_face && vertex != 0 && iter_guard < 16
            iter_guard += 1
            vertex, next_face = _ses_build_edge!(ses, cache, edge_idx,
                                                    start_face, next_face, target_face,
                                                    vertex, indices, vertex_grid,
                                                    true)
        end
        if next_face != target_face
            next_face = actual_max
            vertex = v2_idx
            iter_guard = 0
            while next_face != target_face && vertex != 0 && iter_guard < 16
                iter_guard += 1
                vertex, next_face = _ses_build_edge!(ses, cache, edge_idx,
                                                        start_face, next_face, target_face,
                                                        vertex, indices, vertex_grid,
                                                        false)
            end
        end
    end
    walk_one_pass(face1, face2)
    walk_one_pass(face2, face1)
    return true   # mark for deletion
end

# Top-level treatSecondCategory. Iterates singular SES edges, splits each
# at probe-probe intersections, and marks the original for deletion.
function _ses_treat_second_category!(ses::SolventExcludedSurface{T}) where T
    isempty(ses.singular_edges) && return ses
    pr = ses.reduced_surface.probe_radius
    before_count = length(ses.singular_edges)
    cache = Dict{NTuple{3,Int}, Any}()
    vertex_grid = Dict{NTuple{3,Int}, Vector{Int}}()
    # Seed vertex grid with existing SES vertices.
    for (vi, v) in enumerate(ses.vertices)
        _insert_vertex_grid!(vertex_grid, v.point, vi)
    end
    # BALL `treatSecondCategory` (solventExcludedSurface.C:2059-2064)
    # iterates `singular_edges_` LIVE — `buildEdge`'s `push_back(new_edge)`
    # appends to the same list during iteration, so newly-created singular
    # edges get re-treated in the same pass. The recursion terminates
    # because each new edge eventually finds no more intersections.
    deletable = Int[]
    # BALL `buildEndEdge` (solventExcludedSurface.C:2433) uses `push_front`
    # so the new shortened end-edges are NOT re-visited in this pass. My
    # port pushes back; track those eids and skip them. Mirrors BALL's
    # `treatSecondCategory` semantics: only ORIGINAL singular edges are
    # split; the resulting end-edges are kept live but not re-treated.
    skip_set = Set{Int}()
    i = 1
    while i <= length(ses.singular_edges)
        eid = ses.singular_edges[i]
        if eid in skip_set
            i += 1
            continue
        end
        edge = ses.edges[eid]
        if edge.type == SESEdgeType.Singular && edge.v1 != 0
            if _ses_treat_singular_edge!(ses, eid, Vector{Int}[], cache,
                                           vertex_grid, pr, skip_set)
                push!(deletable, eid)
            end
        end
        i += 1
    end
    after_treat = length(ses.singular_edges)
    # Mark deletable edges by clearing their endpoints (don't actually
    # remove from arrays to keep indices stable).
    for eid in deletable
        e = ses.edges[eid]
        e.v1 = 0; e.v2 = 0
        # Remove from singular_edges list.
        idx = findfirst(==(eid), ses.singular_edges)
        idx !== nothing && deleteat!(ses.singular_edges, idx)
        # Remove from spheric face edges.
        for sf in ses.spheric_faces
            i = findfirst(==(eid), sf.edges)
            i !== nothing && deleteat!(sf.edges, i)
        end
    end
    if haskey(ENV, "SES_DEBUG_TSC")
        println(stderr, "[TSC] singular before=$before_count after=$after_treat added=$(after_treat - before_count) deletable=$(length(deletable))")
    end
    ses
end
