# ===========================================================================
# Port of BALL's `SESTriangulator` class (triangulatedSES.C). This file
# replaces the per-RS-element pipeline in `advancing_front.jl` + parts of
# `triangulated_ses.jl` with a BALL-faithful traversal driven by the SES
# face topology that `compute_ses` builds (edges in cyclic boundary order,
# edge↔face cross-references, singular SES edges from
# treatSingularToricFace).
#
# BALL's flow (SESTriangulator::run):
#   1. preProcessing: clean() and splitSphericFaces() on the SES; create
#      one TrianglePoint per SES vertex.
#   2. triangulateToricFaces: dispatch per toric face — free / singular /
#      non-singular via different buildTriangles calls. Sets EPSILON = 1e-4
#      during this phase (looser comparison for toric arc snapping).
#   3. partitionSingularEdges: sample each singular SES edge into a chain
#      of TriangleEdges between contact faces.
#   4. triangulateContactFaces: for each contact face, copy + cut the
#      template icosphere against the face's bounding edges, then run
#      buildSphericTriangles (advancing-front, convex = true).
#   5. triangulateSphericFaces: for each spheric face, run
#      buildSphericTriangles (no icosphere — boundary-only, convex = false).
#
# The whole flow operates on `SESFace::edge_` lists (cyclic, normalized).
# Our `ses.contact_faces[i].edges` / `.toric_faces[i].edges` /
# `.spheric_faces[i].edges` provide the equivalent — populated by
# `compute_ses` and ordered by `_normalize_all_ses_faces!`.
# ===========================================================================

# Mirrors BALL's TrianglePoint — a vertex on the working surface with
# `.point`, `.normal`, and incidence sets to edges / faces. `global_id`
# is the index in the final mesh's vertex array.
mutable struct TPoint{T<:Real}
    point::Vector3{T}
    normal::Vector3{T}
    edges::Set{Int}
    faces::Set{Int}
    global_id::Int
end

TPoint{T}(p::Vector3{T}, n::Vector3{T}, gid::Int) where T =
    TPoint{T}(p, n, Set{Int}(), Set{Int}(), gid)

# Mirrors BALL's TriangleEdge — directed (vertex_[0] → vertex_[1]).
mutable struct TEdge
    v0::Int
    v1::Int
    f0::Int     # incident triangles (0 = none)
    f1::Int
end

# Mirrors BALL's Triangle — vertices in CCW outward order, edges in
# convention (edge_[0] between vertex_[0]/[1], edge_[1] between [1]/[2],
# edge_[2] between [2]/[0]).
mutable struct TTri
    v0::Int
    v1::Int
    v2::Int
    e01::Int
    e12::Int
    e20::Int
end

# Container for BALL's `TriangulatedSES` (== `part` inside each face
# triangulator). Mirrors BALL's per-class state.
mutable struct TSES{T<:Real}
    points::Vector{TPoint{T}}
    edges::Vector{TEdge}
    tris::Vector{TTri}
end
TSES{T}() where T = TSES{T}(TPoint{T}[], TEdge[], TTri[])

# Mirrors BALL's `SESTriangulator`. Holds:
#   * `ses` — the SES being triangulated.
#   * `density` (and `sqrt_density`) — triangulation density parameter.
#   * `tses` — the global triangulated mesh accumulating all face parts.
#   * `point` — array indexed by SES vertex index → its TPoint in `tses`.
#     Mirrors BALL's `point_[]`. Populated by preProcessing.
#   * `edge` — per-SES-edge list of TriangleEdges; mirrors BALL's
#     `edge_[]`. Populated as toric faces produce boundary edges, then
#     read by buildSphericTriangles to seed the candidate set.
#   * `template_spheres` — icosphere refinement templates indexed by
#     refinement level (0..3). Mirrors BALL's `template_spheres_`.
mutable struct SESTriangulator{T<:Real}
    ses::SolventExcludedSurface{T}
    density::T
    sqrt_density::T
    tses::TSES{T}
    # SES-vertex-index → TPoint index in tses.points.
    point::Vector{Int}
    # SES-edge-index → list of TEdge indices in tses.edges that sample
    # this SES edge (created by toric face triangulation).
    edge::Vector{Vector{Int}}
    # Icosphere templates, lazily extended on demand. `templates[k+1]` is the
    # icosphere refined `k` times. BALL `numberOfRefinements`
    # (triangulatedSES.C:1516) is unbounded; for high-density / large-radius
    # inputs the refinement level can exceed 3, so we extend dynamically.
    templates::Vector{Any}
end

# Convert density to BALL's refinement-level integer. Mirrors
# `SESTriangulator::numberOfRefinements` (triangulatedSES.C:1350).
function ball_num_refinements(density::T, radius::T) where T
    EPS = T(1e-6)
    test0 = (T(4) * density * T(π) * radius^2 - T(12)) / T(30)
    n = 0
    # BALL `Maths::isGreaterOrEqual(test0, 0)`: test0 >= -EPSILON. Bail
    # only when strictly below tolerance.
    test0 < -EPS && return 0
    test1 = one(T); test2 = one(T)
    # BALL `Maths::isLess(test2, test0)`: test2 < test0 - EPSILON.
    while test2 - test0 < -EPS
        test1 = test2
        test2 *= T(4)
        n += 1
    end
    # BALL `Maths::isLess(test2-test0, test0-test1)`.
    if (test2 - test0) - (test0 - test1) < -EPS
        n += 1
    end
    # BALL `numberOfRefinements` (triangulatedSES.C:1369-1372) clamps to
    # max 3. Earlier comment incorrectly said BALL is unbounded. Match
    # BALL exactly.
    n > 3 ? 3 : n
end

# Look up (or lazily build) the icosphere template at refinement level
# `n_ref` (0-indexed). Mirrors BALL's `template_spheres_` HashMap lookup
# in triangulatedSES.C:845-859 — BALL pre-builds templates in
# `buildTemplateSpheres` but lazily-extending here gives the same effect.
@inline function _ball_template_at(tr::SESTriangulator{T}, n_ref::Int) where T
    while length(tr.templates) <= n_ref
        push!(tr.templates, icosphere(T, length(tr.templates)))
    end
    tr.templates[n_ref + 1]
end

# Port of BALL's `SESTriangulator::partitionOfCircle`
# (triangulatedSES.C:1317). Samples `number_of_segments + 1` points
# along a circle starting at `p0` (relative to `circle.p` if
# `on_surface=true`, otherwise picks a default starting direction
# perpendicular to `circle.n`). Each step rotates by `phi` around
# `circle.n` via Rodrigues' formula.
#
# Returns a Vector{Vector3{T}} of length n_seg+1 — the first element
# is at p0, the last after n_seg rotations of psi.
function ball_partition_of_circle(circle::Circle3{T}, p0::Vector3{T},
                                   phi::T, number_of_segments::Int;
                                   on_surface::Bool = true) where T
    if on_surface
        p = p0 - circle.p
    else
        # Pick any direction perpendicular to circle.n.
        p = Vector3{T}(circle.n[2], -circle.n[1], zero(T))
        if dot(p, p) < eps(T)
            p = Vector3{T}(circle.n[3], zero(T), -circle.n[1])
        end
        ℓ = sqrt(dot(p, p))
        p = (p / ℓ) * circle.radius
    end
    out = Vector{Vector3{T}}(undef, number_of_segments + 1)
    out[1] = circle.p + p
    # Rodrigues' rotation: v' = v cos θ + (k × v) sin θ + k (k · v)(1 - cos θ)
    # requires `k` to be a UNIT vector. BALL's quaternion rotation
    # (TQuaternion::fromAxisAngle, quaternion.h:fromAxisAngle) normalizes
    # the axis internally; match here so callers don't need to ensure
    # `circle.n` is normalized.
    cosφ = cos(phi); sinφ = sin(phi)
    nk = sqrt(dot(circle.n, circle.n))
    k = nk > eps(T) ? circle.n / nk : circle.n
    for i in 1:number_of_segments
        kxp = cross(k, p)
        kdotp = dot(k, p)
        p_new = p * cosφ + kxp * sinφ + k * (kdotp * (1 - cosφ))
        p = p_new
        out[i + 1] = circle.p + p
    end
    out
end

# Sample an SES edge into `n_seg + 1` points along its carrier arc.
# Uses `ball_partition_of_circle` to walk from vertex[0] to vertex[1]
# in the direction of the carrier circle's normal. The endpoint is
# snapped to vertex[1].point to handle numerical drift (mirrors
# BALL's `edge_segments[number_of_segments] = edge->vertex_[1]->point_`
# in triangulateNonSingularToricFace, triangulatedSES.C:311 etc.).
function ball_sample_ses_edge(ses::SolventExcludedSurface{T},
                               se_idx::Int, n_seg::Int) where T
    se = ses.edges[se_idx]
    v0 = ses.vertices[se.v1].point
    v1 = ses.vertices[se.v2].point
    # Compute the oriented angle from v0 to v1 around circle.n.
    rel0 = v0 - se.circle.p
    rel1 = v1 - se.circle.p
    # Build an orthonormal basis (u, v) in the circle's plane.
    u = rel0
    ℓu = sqrt(dot(u, u))
    u = ℓu > eps(T) ? u / ℓu : Vector3{T}(1, 0, 0)
    v = cross(se.circle.n, u)
    ψ_e = atan(dot(rel1, v), dot(rel1, u))
    Δ = ψ_e
    # Wrap into the shorter arc consistent with the circle's normal
    # direction (BALL relies on caller-provided orientation via
    # circle.n, so we don't pick "short vs long" here).
    if Δ < zero(T); Δ += 2 * T(π); end
    psi = Δ / n_seg
    out = ball_partition_of_circle(se.circle, v0, psi, n_seg)
    out[end] = v1   # snap last sample
    out
end

# Helper: add a new TEdge with v0/v1 to `tr.tses`, link to both endpoints'
# incidence sets, and optionally record on a per-SES-edge edge list
# (`tr.edge[se_idx]`). Mirrors BALL's pattern of `new TriangleEdge` +
# `tses_->edges_.push_back(...)` + `vertex_[*]->edges_.insert(...)` +
# `edge_[se_idx].push_back(...)`.
@inline function _ball_new_edge!(tr::SESTriangulator{T}, v0::Int, v1::Int;
                                 ses_edge::Int = 0) where T
    push!(tr.tses.edges, TEdge(v0, v1, 0, 0))
    eid = length(tr.tses.edges)
    push!(tr.tses.points[v0].edges, eid)
    push!(tr.tses.points[v1].edges, eid)
    ses_edge > 0 && push!(tr.edge[ses_edge], eid)
    eid
end

# Helper: add a new TTri with vertices v0/v1/v2 and edges e01/e12/e20.
# Sets face_[0] = this triangle on all three edges (mirrors BALL's
# `t->edge_[i]->face_[0] = t` after construction) and inserts the
# triangle into all three vertices' incidence sets.
@inline function _ball_new_tri!(tr::SESTriangulator{T},
                                v0::Int, v1::Int, v2::Int,
                                e01::Int, e12::Int, e20::Int) where T
    # Note: BALL `SESTriangulator::buildTriangles` doesn't guard against
    # degenerate triangles. A previous guard here (refuse if any two
    # corners coincide) was a band-aid over upstream geometric
    # divergences — removed per the BALL-faithful port directive.
    push!(tr.tses.tris, TTri(v0, v1, v2, e01, e12, e20))
    tid = length(tr.tses.tris)
    push!(tr.tses.points[v0].faces, tid)
    push!(tr.tses.points[v1].faces, tid)
    push!(tr.tses.points[v2].faces, tid)
    e = tr.tses.edges
    e[e01].f0 = e[e01].f0 == 0 ? tid : (e[e01].f1 = tid; e[e01].f0)
    e[e12].f0 = e[e12].f0 == 0 ? tid : (e[e12].f1 = tid; e[e12].f0)
    e[e20].f0 = e[e20].f0 == 0 ? tid : (e[e20].f1 = tid; e[e20].f0)
    tid
end

# Port of BALL's `SESTriangulator::buildTriangles` (triangulatedSES.C:1410)
# for the regular (non-singular, non-free) case. Builds a
# (n_seg+1) × (n_tri+1) grid of TrianglePoints between edge1 and edge3,
# with `centers[i]` being the rolling axis position at step i.
#
# Arguments:
#   * edge0, edge1, edge2, edge3 — the 4 SES edges of the toric face
#     in BALL's normalize order (0 = first convex / contact arc,
#     1 = concave arc on one side, 2 = second convex, 3 = concave on
#     other side). For our SES topology each is an SES edge index.
#   * centers — n_seg+1 Vector3{T} points sampled along the rolling
#     axis (centers of the cross-section probe sphere at each step).
#   * edge1_segments, edge3_segments — n_seg+1 sample points along the
#     two concave boundary arcs.
#   * probe_radius — the rolling probe radius.
#
# Emits n_seg × n_tri × 2 triangles into `tr.tses`.
function ball_build_triangles_regular!(tr::SESTriangulator{T},
                                        edge0::Int, edge1::Int,
                                        edge2::Int, edge3::Int,
                                        centers::Vector{Vector3{T}},
                                        edge1_segments::Vector{Vector3{T}},
                                        edge3_segments::Vector{Vector3{T}},
                                        probe_radius::T) where T
    # Compute n_tri from the angle between (edge1[0]-centers[0]) and
    # (edge3[0]-centers[0]) seen from centers[0].
    a = edge1_segments[1] - centers[1]
    b = edge3_segments[1] - centers[1]
    cosθ = dot(a, b) / (sqrt(dot(a, a)) * sqrt(dot(b, b)))
    cosθ = clamp(cosθ, -one(T), one(T))
    ψ = acos(cosθ)
    n_tri = max(1, round(Int, ψ * probe_radius * tr.sqrt_density, RoundNearestTiesAway))
    phi = ψ / n_tri
    n_seg = length(centers) - 1
    # Grid indexing: grid[i+1, j+1] holds the TPoint index for cross-
    # section step `j` (rolling along centers) and circle step `i`
    # (across the strip). 1-based throughout — corresponds to BALL's
    # `points[i*end + j]` with `end = n_tri+1`.
    grid = Matrix{Int}(undef, n_seg + 1, n_tri + 1)
    e0 = tr.ses.edges[edge0]
    e2 = tr.ses.edges[edge2]
    for i in 0:n_seg
        # Cross-section circle at centers[i+1] in the plane spanned by
        # (edge1_points[i] - centers[i]) and (edge3_points[i] - centers[i]).
        u = edge1_segments[i + 1] - centers[i + 1]
        v = edge3_segments[i + 1] - centers[i + 1]
        n_circle = cross(u, v)
        ℓ = sqrt(dot(n_circle, n_circle))
        n_circle = ℓ > eps(T) ? n_circle / ℓ : Vector3{T}(0, 0, 1)
        circle = Circle3{T}(centers[i + 1], n_circle, probe_radius)
        line = ball_partition_of_circle(circle, edge1_segments[i + 1], phi, n_tri)
        for j in 0:n_tri
            # Reuse an existing SES corner TPoint at the 4 grid corners.
            corner_pid = 0
            if i == 0 && j == 0
                d1 = line[1] - tr.ses.vertices[e0.v1].point
                d2 = line[1] - tr.ses.vertices[e0.v2].point
                corner_pid = tr.point[dot(d1, d1) < dot(d2, d2) ? e0.v1 : e0.v2]
            elseif i == 0 && j == n_tri
                d1 = line[n_tri + 1] - tr.ses.vertices[e0.v1].point
                d2 = line[n_tri + 1] - tr.ses.vertices[e0.v2].point
                corner_pid = tr.point[dot(d1, d1) < dot(d2, d2) ? e0.v1 : e0.v2]
            elseif i == n_seg && j == 0
                d1 = line[1] - tr.ses.vertices[e2.v1].point
                d2 = line[1] - tr.ses.vertices[e2.v2].point
                corner_pid = tr.point[dot(d1, d1) < dot(d2, d2) ? e2.v1 : e2.v2]
            elseif i == n_seg && j == n_tri
                d1 = line[n_tri + 1] - tr.ses.vertices[e2.v1].point
                d2 = line[n_tri + 1] - tr.ses.vertices[e2.v2].point
                corner_pid = tr.point[dot(d1, d1) < dot(d2, d2) ? e2.v1 : e2.v2]
            end
            pid = if corner_pid == 0
                # Create a fresh interior point.
                p = line[j + 1]
                push!(tr.tses.points, TPoint{T}(p, centers[i + 1] - p,
                                                 0))  # global_id = 0 (interior)
                length(tr.tses.points)
            else
                corner_pid
            end
            grid[i + 1, j + 1] = pid
        end
    end
    # Build vertical edges: connect grid[j, i] ↔ grid[j, i+1] for i ∈ 0..n_tri-1.
    # In BALL these are the "edges[offset2 - j]" entries.
    vedges = Matrix{Int}(undef, n_seg + 1, n_tri)
    for j in 0:n_seg, i in 0:(n_tri - 1)
        se = 0
        j == 0    && (se = edge0)
        j == n_seg && (se = edge2)
        vedges[j + 1, i + 1] = _ball_new_edge!(tr,
            grid[j + 1, i + 1], grid[j + 1, i + 2]; ses_edge = se)
    end
    # Horizontal edges: grid[j, i] ↔ grid[j+1, i] for j ∈ 0..n_seg-1.
    hedges = Matrix{Int}(undef, n_seg, n_tri + 1)
    for j in 0:(n_seg - 1), i in 0:n_tri
        se = 0
        i == 0     && (se = edge1)
        i == n_tri && (se = edge3)
        hedges[j + 1, i + 1] = _ball_new_edge!(tr,
            grid[j + 1, i + 1], grid[j + 2, i + 1]; ses_edge = se)
    end
    # Diagonal edges: grid[j, i] ↔ grid[j+1, i+1] for j ∈ 0..n_seg-1, i ∈ 0..n_tri-1.
    dedges = Matrix{Int}(undef, n_seg, n_tri)
    for j in 0:(n_seg - 1), i in 0:(n_tri - 1)
        dedges[j + 1, i + 1] = _ball_new_edge!(tr,
            grid[j + 1, i + 1], grid[j + 2, i + 2])
    end
    # Build 2 triangles per cell (j, i).
    for i in 0:(n_tri - 1), j in 0:(n_seg - 1)
        # t1: grid[j, i], grid[j, i+1], grid[j+1, i+1]
        _ball_new_tri!(tr,
            grid[j + 1, i + 1], grid[j + 1, i + 2], grid[j + 2, i + 2],
            vedges[j + 1, i + 1], hedges[j + 1, i + 2], dedges[j + 1, i + 1])
        # t2: grid[j, i], grid[j+1, i+1], grid[j+1, i]
        _ball_new_tri!(tr,
            grid[j + 1, i + 1], grid[j + 2, i + 2], grid[j + 2, i + 1],
            dedges[j + 1, i + 1], vedges[j + 2, i + 1], hedges[j + 1, i + 1])
    end
    nothing
end

# Mirrors BALL's `SESTriangulator::preProcessing` (triangulatedSES.C:144).
# Note: BALL also runs `ses.clean(density)` and `ses.splitSphericFaces()`
# here; in our pipeline those happen earlier (in `triangulate_ses`).
# Creates one TPoint per SES vertex and populates `tr.point` to map
# SES vertex indices to TPoint indices. Also initializes `tr.edge` with
# one empty list per SES edge (populated by toric face triangulation).
function ball_pre_processing!(tr::SESTriangulator{T}) where T
    rs = tr.ses.reduced_surface
    n_v = length(tr.ses.vertices)
    resize!(tr.point, n_v)
    for i in 1:n_v
        sv = tr.ses.vertices[i]
        # BALL stores SES vertex normals as UNNORMALIZED vectors
        # (SESComputer::createContactVertex etc., e.g. solventExcludedSurface.C:1215:
        # `vertex->normal_.set(probe_center - vertex->point_)`). The
        # commented-out post-processing loop in `SESTriangulator::run`
        # (triangulatedSES.C:179-184) would normalize but is disabled.
        # Compute unnormalized outward direction for parity.
        normal = if sv.atom != 0 && sv.rs_face != 0
            # Contact corner: outward from atom (= toward probe).
            sv.point - Vector3{T}(rs.atoms[sv.atom].center...)
        elseif sv.rs_edge != 0
            # Singular axis-intersection vertex: BALL's
            # createSingularVertex stores `intersection_circle.p - point`
            # (= toward the rolling axis, inward in the saddle).
            re = rs.edges[sv.rs_edge]
            if re.f1 != 0 && re.f2 != 0
                midp = (rs.faces[re.f1].center + rs.faces[re.f2].center) / T(2)
                midp - sv.point
            else
                Vector3{T}(0, 0, 1)
            end
        else
            Vector3{T}(0, 0, 1)
        end
        push!(tr.tses.points, TPoint{T}(sv.point, normal, i))
        tr.point[i] = length(tr.tses.points)
    end
    n_e = length(tr.ses.edges)
    resize!(tr.edge, n_e)
    for i in 1:n_e
        tr.edge[i] = Int[]
    end
end

# Port of BALL's `SESTriangulator::partitionSingularEdges`
# (triangulatedSES.C:174) — samples each singular SES edge into a chain
# of `n+1` TPoints and `n` TEdges along its intersection-circle arc,
# stored in `tr.edge[singular_eid]`. The two endpoint TPoints are the
# already-existing TPoints at the singular edge's two SES vertices
# (axis-intersection points created by `_treat_singular_toric_faces!`);
# the n-1 interior TPoints are new.
#
# After this runs, the spheric face triangulator sees these TEdges via
# `tr.edge[singular_eid]` (because `_treat_singular_toric_faces!` added
# the singular edge to each spheric face's `.edges` list), giving it
# boundary samples to weld with the singular-toric mesh.
#
# Mirrors `partitionNonFreeSingularEdge` (triangulatedSES.C:617). Free
# variant (vertex[0] == NULL — for isolated 2-atom systems) is not yet
# ported.
function ball_partition_singular_edges!(tr::SESTriangulator{T}) where T
    for se_idx in tr.ses.singular_edges
        se = tr.ses.edges[se_idx]
        se.type == SESEdgeType.Singular || continue
        if se.v1 == 0 || se.v2 == 0
            # Full-circle singular edge (from `noCut`): sample n points
            # around 2π and connect as a closed ring of n TEdges.
            # BALL `partitionFreeSingularEdge` (triangulatedSES.C:691-698):
            # n = round(2π·r·√density); if n==0, n++ (= 1). No min-3 clamp.
            n = max(1, round(Int, 2 * T(π) * se.circle.r * tr.sqrt_density, RoundNearestTiesAway))
            # Pick a unit vector perpendicular to circle.n. BALL's
            # `partitionOfCircle` (triangulatedSES.C:1494-1500) with
            # `on_surface=false`: `orth = (circle.n.y, -circle.n.x, 0)`, fall
            # back to `(circle.n.z, 0, -circle.n.x)` only if first is exactly
            # zero (= circle.n is along z axis). Match BALL exactly.
            cn = se.circle.n
            ortho = Vector3{T}(cn[2], -cn[1], zero(T))
            if ortho == Vector3{T}(0, 0, 0)
                ortho = Vector3{T}(cn[3], zero(T), -cn[1])
            end
            ℓo = sqrt(dot(ortho, ortho))
            ortho = ℓo > eps(T) ? ortho / ℓo : Vector3{T}(1, 0, 0)
            p0 = se.circle.p + ortho * se.circle.r
            pts = ball_partition_of_circle(se.circle, p0, 2 * T(π) / n, n)
            # First n samples form the ring (last sample equals first).
            ring_pids = Vector{Int}(undef, n)
            for k in 1:n
                # BALL: `p->normal_ = singular_edge->circle.p - p->point_`
                # (triangulatedSES.C:711, 776), unnormalized.
                normal = se.circle.p - pts[k]
                push!(tr.tses.points, TPoint{T}(pts[k], normal, 0))
                ring_pids[k] = length(tr.tses.points)
            end
            for k in 1:n
                k_next = mod1(k + 1, n)
                _ball_new_edge!(tr, ring_pids[k], ring_pids[k_next];
                                ses_edge = se_idx)
            end
            continue
        end
        v0_p = tr.ses.vertices[se.v1].point
        v1_p = tr.ses.vertices[se.v2].point
        # Oriented angle from v0 to v1 around circle.n.
        rel0 = v0_p - se.circle.p
        rel1 = v1_p - se.circle.p
        u = rel0
        ℓu = sqrt(dot(u, u))
        u = ℓu > eps(T) ? u / ℓu : Vector3{T}(1, 0, 0)
        vbasis = cross(se.circle.n, u)
        ψ = atan(dot(rel1, vbasis), dot(rel1, u))
        if ψ < zero(T); ψ += 2 * T(π); end
        n_seg = max(1, round(Int, ψ * se.circle.r * tr.sqrt_density, RoundNearestTiesAway))
        psi = ψ / n_seg
        # Sample n+1 points starting from v0_p, snap last to v1_p.
        pts = ball_partition_of_circle(se.circle, v0_p, psi, n_seg)
        pts[end] = v1_p
        # The 2 endpoint TPoints already exist; create n-1 interior
        # TPoints and n TEdges between consecutive samples. Update
        # endpoint normals to point toward the singular circle center
        # (BALL: p1->normal_ = singular_edge->circle_.p - p1->point_).
        # BALL `partitionNonFreeSingularEdge` (triangulatedSES.C:770-790):
        # endpoint TPoints get `normal_ = singular_edge->circle.p - point`
        # (unnormalized). Interior TPoints same formula.
        prev_pid = tr.point[se.v1]
        let tp = tr.tses.points[prev_pid]
            tp.normal = se.circle.p - tp.point
        end
        for k in 2:(n_seg)
            normal = se.circle.p - pts[k]
            push!(tr.tses.points, TPoint{T}(pts[k], normal, 0))
            mid_pid = length(tr.tses.points)
            _ball_new_edge!(tr, prev_pid, mid_pid; ses_edge = se_idx)
            prev_pid = mid_pid
        end
        # Final TEdge: prev_pid → endpoint TPoint at se.v2.
        end_pid = tr.point[se.v2]
        let tp = tr.tses.points[end_pid]
            tp.normal = se.circle.p - tp.point
        end
        _ball_new_edge!(tr, prev_pid, end_pid; ses_edge = se_idx)
    end
end

# Port of BALL's `buildTriangles` variant for singular toric halves
# (edge3 == NULL in BALL). Each half-strip is a triangular patch in
# topology, bounded by:
#   * `edge0_se` — concave arc on a probe sphere (n_seg+1 samples).
#   * `edge1_se` — split contact arc on one atom (n_seg+1 samples,
#     last snapped to the apex).
#   * `edge2_se` — split contact arc on the other atom (n_seg+1 samples).
# The "edge3" collapses to a single apex vertex (vertex[1+offset] in
# BALL). Each cell at the top row degenerates: t1 fans from
# (j, n_tri-1) to apex, t2 is zero-area (apex == apex).
function ball_build_triangles_singular!(tr::SESTriangulator{T},
                                          edge0_se::Int, edge1_se::Int,
                                          edge2_se::Int,
                                          centers::Vector{Vector3{T}},
                                          edge1_segments::Vector{Vector3{T}},
                                          apex_point::Vector3{T},
                                          apex_pid::Int,
                                          probe_radius::T,
                                          apex_sv::Int = 0) where T
    n_seg = length(centers) - 1
    # n_tri from the angle between (edge1[0] - centers[0]) and
    # (apex - centers[0]) seen from centers[0].
    a = edge1_segments[1] - centers[1]
    b = apex_point - centers[1]
    cosθ = dot(a, b) / (sqrt(dot(a, a)) * sqrt(dot(b, b)))
    cosθ = clamp(cosθ, -one(T), one(T))
    ψ = acos(cosθ)
    n_tri = max(1, round(Int, ψ * probe_radius * tr.sqrt_density, RoundNearestTiesAway))
    phi = ψ / n_tri
    # Grid: n_seg+1 columns × n_tri rows + 1 apex point.
    grid = Matrix{Int}(undef, n_seg + 1, n_tri)
    # BALL `buildTriangles` (triangulatedSES.C:1722-1751) corner check uses
    # `edge0` (= function parameter) for (i==0, j==0) and `edge2` for
    # (i==n_seg, j==0). For singular toric, edge0 = CONCAVE_0 and edge2 =
    # CONCAVE_2 (passed through edge0_se / edge2_se from process_cycle).
    # The check finds the closer of those edges' endpoints to line[0] —
    # which is the corner shared with the CONVEX edge (since we now sample
    # CONVEX starting at the shared corner).
    e0 = tr.ses.edges[edge0_se]
    e2 = tr.ses.edges[edge2_se]
    for i in 0:n_seg
        u = edge1_segments[i + 1] - centers[i + 1]
        v = apex_point - centers[i + 1]
        n_circle = cross(u, v)
        ℓ = sqrt(dot(n_circle, n_circle))
        n_circle = ℓ > eps(T) ? n_circle / ℓ : Vector3{T}(0, 0, 1)
        circle = Circle3{T}(centers[i + 1], n_circle, probe_radius)
        line = ball_partition_of_circle(circle, edge1_segments[i + 1], phi, n_tri)
        for j in 0:(n_tri - 1)
            corner_pid = 0
            if i == 0 && j == 0
                d1 = line[1] - tr.ses.vertices[e0.v1].point
                d2 = line[1] - tr.ses.vertices[e0.v2].point
                corner_sv = dot(d1, d1) < dot(d2, d2) ? e0.v1 : e0.v2
                corner_pid = tr.point[corner_sv]
            elseif i == n_seg && j == 0
                d1 = line[1] - tr.ses.vertices[e2.v1].point
                d2 = line[1] - tr.ses.vertices[e2.v2].point
                corner_sv = dot(d1, d1) < dot(d2, d2) ? e2.v1 : e2.v2
                corner_pid = tr.point[corner_sv]
            end
            pid = if corner_pid == 0
                p = line[j + 1]
                push!(tr.tses.points, TPoint{T}(p, centers[i + 1] - p, 0))
                length(tr.tses.points)
            else
                corner_pid
            end
            grid[i + 1, j + 1] = pid
        end
    end
    # Edge recording follows BALL's `buildTriangles` (triangulatedSES.C:1410):
    # in our (rolling-axis × cross-section) indexing, j varies along the
    # rolling axis (= "j" in BALL) and i along the cross-section
    # (= "i" in BALL). Vertical edges connect same-j adjacent-i grid
    # points and at j=0 / j=n_seg lie on edge0_se / edge2_se (the two
    # concave arcs of the cycle). Horizontal edges connect adjacent-j
    # same-i and at i=0 lie on edge1_se (the convex arc, = the "j=0"
    # edges in BALL terms — wait no, in BALL it's i=0 too; the patch's
    # cross-section step 0 traces edge1's circle).
    vedges = Matrix{Int}(undef, n_seg + 1, n_tri - 1)
    for j in 0:n_seg, i in 0:(n_tri - 2)
        se = 0
        j == 0      && (se = edge0_se)   # bottom concave (probe 1)
        j == n_seg  && (se = edge2_se)   # top concave (probe 2)
        vedges[j + 1, i + 1] = _ball_new_edge!(tr,
            grid[j + 1, i + 1], grid[j + 1, i + 2]; ses_edge = se)
    end
    hedges = Matrix{Int}(undef, n_seg, n_tri)
    for j in 0:(n_seg - 1), i in 0:(n_tri - 1)
        se = 0
        i == 0 && (se = edge1_se)   # convex arc (atom side, bottom row)
        hedges[j + 1, i + 1] = _ball_new_edge!(tr,
            grid[j + 1, i + 1], grid[j + 2, i + 1]; ses_edge = se)
    end
    dedges = Matrix{Int}(undef, n_seg, n_tri - 1)
    for j in 0:(n_seg - 1), i in 0:(n_tri - 2)
        dedges[j + 1, i + 1] = _ball_new_edge!(tr,
            grid[j + 1, i + 1], grid[j + 2, i + 2])
    end
    # Build 2 triangles per cell for the regular rows (i = 0..n_tri-2).
    for i in 0:(n_tri - 2), j in 0:(n_seg - 1)
        _ball_new_tri!(tr,
            grid[j + 1, i + 1], grid[j + 1, i + 2], grid[j + 2, i + 2],
            vedges[j + 1, i + 1], hedges[j + 1, i + 2], dedges[j + 1, i + 1])
        _ball_new_tri!(tr,
            grid[j + 1, i + 1], grid[j + 2, i + 2], grid[j + 2, i + 1],
            dedges[j + 1, i + 1], vedges[j + 2, i + 1], hedges[j + 1, i + 1])
    end
    # Apex fan — mirrors BALL's `buildTriangles` singular branch
    # (triangulatedSES.C:1714-1757). Each cell at i = n_tri-1 produces
    # ONE triangle:  t = (grid[j, n_tri-1], grid[j+1, n_tri-1], apex).
    # Edge attribution to `tr.edge[ses_edge]`:
    #   * the FIRST apex edge (from grid[0, n_tri-1] to apex) → edge0_se
    #     (the concave arc on probe 1, = BALL's edge0 parameter)
    #   * the LAST apex edge (from grid[n_seg, n_tri-1] to apex) → edge2_se
    #     (the concave arc on probe 2)
    #   * all middle apex edges → no SES edge (purely internal to the
    #     singular toric patch's apex fan).
    apex_edges = Vector{Int}(undef, n_seg + 1)
    for j in 0:n_seg
        se_for_apex = if j == 0
            edge0_se
        elseif j == n_seg
            edge2_se
        else
            0
        end
        apex_edges[j + 1] = _ball_new_edge!(tr, grid[j + 1, n_tri], apex_pid;
                                             ses_edge = se_for_apex)
    end
    # BALL's apex-fan triangle order (triangulatedSES.C:1729-1755):
    #   t.vertex = (grid[i, n_tri-1], apex, grid[i+1, n_tri-1])
    #   t.edge[0] = apex_edge[i] (= grid[i] ↔ apex)
    #   t.edge[1] = horizontal (= grid[i] ↔ grid[i+1])
    #   t.edge[2] = apex_edge[i+1] (= grid[i+1] ↔ apex)
    # The (left, apex, right) ordering gives a right-hand normal
    # consistent with the regular-cell triangles below.
    for j in 0:(n_seg - 1)
        _ball_new_tri!(tr,
            grid[j + 1, n_tri], apex_pid, grid[j + 2, n_tri],
            apex_edges[j + 1], apex_edges[j + 2], hedges[j + 1, n_tri])
    end
    nothing
end

# Port of BALL's `triangulateSingularToricFace` (triangulatedSES.C:433).
# For a singular toric face (6 vertices, 6 edges arranged in 2
# disconnected triangular cycles after `treat_singular_toric_faces!`):
# build each half as a triangular patch fanning to its axis-intersection
# apex vertex.
function ball_triangulate_singular_toric!(tr::SESTriangulator{T},
                                            tf::SESFace{T}) where T
    length(tf.edges) == 6 || return   # treat_singular_toric_faces didn't fully process
    rs = tr.ses.reduced_surface
    re = rs.edges[tf.rs_index]
    (re.f1 == 0 || re.f2 == 0) && return
    # The 6 edges of the singular toric face form 2 disconnected
    # triangular cycles. Walk by vertex adjacency to split.
    edges_list = tf.edges
    verts_list = tf.vertices
    # Cycle 1: starts at vertex[0] (= verts_list[1]). Walk 3 edges.
    # Cycle 2: starts at vertex[3] (= verts_list[4]). Walk 3 edges.
    # In our normalize output, the cycles are contiguous in edges_list /
    # verts_list.
    # Each cycle (in `_treat_singular_toric_faces!`'s output order) has:
    #   - 1 convex edge (atom arc, base of the triangle patch)
    #   - 2 concave edges (probe arcs, the two sides meeting at the apex)
    # The apex is the SES vertex shared by both concave edges (= one of
    # the new axis-intersection vertices).
    function process_cycle(cyc_edges::Vector{Int}, cyc_verts::Vector{Int})
        length(cyc_edges) == 3 || return
        # `_treat_singular_toric_faces!` stores cycle edges as
        # [convex, concave_on_probe_f1, concave_on_probe_f2]. BALL's
        # normalize (SESFace.C findTriangle_, line 327-343) instead puts
        # the concave INCIDENT to convex.vertex_[0] first. Reorder so
        # concave_eids[1] is the concave that shares an endpoint with
        # the convex edge's v1 (= BALL's vertex_[0]).
        convex_idx = cyc_edges[1]
        cv_tmp = tr.ses.edges[convex_idx]
        c_a = tr.ses.edges[cyc_edges[2]]
        c_b = tr.ses.edges[cyc_edges[3]]
        # Find the concave that shares an endpoint with cv_tmp.v1.
        c_a_incident = (c_a.v1 == cv_tmp.v1) || (c_a.v2 == cv_tmp.v1)
        concave_eids = c_a_incident ?
            Int[cyc_edges[2], cyc_edges[3]] :
            Int[cyc_edges[3], cyc_edges[2]]
        # Apex = the vertex shared by both concave edges (not on the convex).
        ce0 = tr.ses.edges[concave_eids[1]]
        ce1 = tr.ses.edges[concave_eids[2]]
        apex_sv = 0
        for v in (ce0.v1, ce0.v2)
            if v == ce1.v1 || v == ce1.v2
                apex_sv = v
                break
            end
        end
        apex_sv == 0 && return
        apex_pid = tr.point[apex_sv]
        apex_point = tr.ses.vertices[apex_sv].point
        # n_seg from rsedge.angle × edge[1]_radius × √density. BALL
        # `triangulateSingularToricFace` (triangulatedSES.C:593-595) uses
        # `edge[1]->circle_.radius` where edge[1] is the SECOND edge in
        # the normalized cycle layout (= concave edge, on the probe sphere).
        # See SESFace.C:265 (`findTriangle_`) — after normalize, the cycle's
        # edges go [convex, concave_1, concave_2], so edge[1] is concave.
        # NOTE: an earlier "fix" briefly used cv.circle.r (= contact circle
        # radius ≈ slightly > probe), which is wrong per BALL but happened
        # to fix downstream BFT failures. Revert to the BALL formula and
        # diagnose the BFT issue separately.
        cc = tr.ses.edges[concave_eids[1]]
        n_seg = max(1, round(Int, re.angle * cc.circle.r * tr.sqrt_density, RoundNearestTiesAway))
        psi = re.angle / n_seg
        # Sample the convex edge (the contact arc) into n_seg+1 points.
        # BALL `triangulateSingularToricFace` (triangulatedSES.C:621-626) samples
        # from `edge[0+offset]->vertex_[0]` along the convex circle. BALL's
        # `findTriangle_` (SESFace.C:323-325) sets `vertex0 = edge0.vertex_[0]`
        # as the corner SHARED with concave_0 (= edge1 in BALL's normalize).
        # We must start sampling at THAT corner — NOT at cv.v1 (which may be
        # the other corner depending on how the convex edge was originally
        # constructed). Identify the shared corner as the non-apex endpoint
        # of concave_eids[1].
        cv = tr.ses.edges[convex_idx]
        ce0 = tr.ses.edges[concave_eids[1]]
        shared_corner = ce0.v1 == apex_sv ? ce0.v2 : ce0.v1
        if cv.v1 == shared_corner
            start_v, end_v = cv.v1, cv.v2
        else
            start_v, end_v = cv.v2, cv.v1
        end
        edge1_segments = ball_partition_of_circle(cv.circle,
            tr.ses.vertices[start_v].point, psi, n_seg)
        edge1_segments[end] = tr.ses.vertices[end_v].point
        # Center circle: axis = the convex edge's normal (the "rolling
        # axis" direction in this cycle's local frame); centered on the
        # torus center, radius = torus major radius. Mirrors BALL's
        # `TCircle3 center_circle(rsedge->center_of_torus, axis, rsedge->radius_of_torus)`.
        axis = cv.circle.n
        # BALL `triangulateSingularToricFace` (triangulatedSES.C:614-625):
        # `if (edge[0+offset]->vertex_[0] != vertex[0+offset]) axis.negate();`
        # — checks if the sampling-start vertex is the same as the cycle's
        # vertex[0+offset]. After my reorientation above, start_v IS the
        # vertex shared with concave_0, which equals BALL's vertex[0+offset].
        # The flip is therefore unnecessary unless cv.circle's normal is
        # oriented opposite to the rolling direction.
        if start_v != cyc_verts[1]
            axis = -axis
        end
        center_circle = Circle3{T}(re.center_of_torus, axis,
                                    re.radius_of_torus)
        # Sample center circle starting from one probe center (the first
        # concave's circle center).
        centers = ball_partition_of_circle(center_circle,
            tr.ses.edges[concave_eids[1]].circle.p, psi, n_seg)
        centers[end] = tr.ses.edges[concave_eids[2]].circle.p
        # Build the triangle patch. BALL's signature:
        # buildTriangles(edge0=concave1, edge1=convex, edge2=concave2, edge3=NULL).
        # Our `ball_build_triangles_singular!` matches this argument order.
        tri_lo = length(tr.tses.tris) + 1
        ball_build_triangles_singular!(tr, concave_eids[1], convex_idx,
            concave_eids[2],
            centers, edge1_segments, apex_point, apex_pid, rs.probe_radius,
            apex_sv)
        # BALL's `triangulateSingularToricFace` (triangulatedSES.C:525-540):
        # if `orth · (v0 - centers[0]) > 0`, flip all triangles emitted
        # by this buildTriangles call.
        if tri_lo <= length(tr.tses.tris)
            test_t = tr.tses.tris[tri_lo]
            v0p = tr.tses.points[test_t.v0].point
            v1p = tr.tses.points[test_t.v1].point
            v2p = tr.tses.points[test_t.v2].point
            orth = cross(v1p - v0p, v2p - v0p)
            # BALL `triangulateSingularToricFace` (triangulatedSES.C:498-512)
            # uses `Maths::isGreater(orth * v0, orth * centers[0])` which with
            # EPSILON=1e-4 (set at line 164) means `(a - b) >= 1e-4`.
            if dot(orth, v0p - centers[1]) >= T(1e-4)
                # BALL only swaps vertex_[0]/vertex_[1] — edge slots are not
                # touched. Match that exactly.
                for ti in tri_lo:length(tr.tses.tris)
                    t = tr.tses.tris[ti]
                    t.v0, t.v1 = t.v1, t.v0
                end
            end
        end
    end
    process_cycle(edges_list[1:3], verts_list[1:3])
    process_cycle(edges_list[4:6], verts_list[4:6])
end

# Port of BALL's `SESTriangulator::triangulateNonSingularToricFace`
# (triangulatedSES.C:272). For a non-free, non-singular toric SES face,
# extract its 4 boundary edges (in normalized cyclic order: edge[0] is
# a convex/contact arc, edge[1] is a concave arc on one probe, edge[2]
# is the other contact arc, edge[3] is the concave arc on the other
# probe), sample edge[1] and edge[3] into `n_seg+1` points each, sample
# the centers along the rolling axis, then delegate to buildTriangles.
#
# NOTE on `length(tf.edges) == 4` gate: BALL's RS construction never
# produces non-singular toric faces with fewer than 4 SES edges (BALL
# crashes on `e++` for fewer). My Julia RS construction generates
# `Vector{Int}` edge lists that, for stacked-probe topologies (3+ probes
# sharing an atom-pair axis from different triples), can yield 3-edge
# toric faces (only 1 concave because the second concave is "owned" by
# a parallel RS edge sharing the same atom pair). These 3-edge tori are
# skipped here, leaving their adjacent convex SES edges unpartitioned in
# `tr.edge[se_idx]` — which becomes the dominant source of BPTI's
# remaining boundary edges (158/189). The faithful fix requires either
# (a) BALL's SESSingularityCleaner::treatSecondCategory to split the
# stacked topology with intersection vertices, or (b) sharing the
# arc partition across multiple toric faces in this triangulator.
function ball_triangulate_non_singular_toric!(tr::SESTriangulator{T},
                                                tf::SESFace{T}) where T
    rs = tr.ses.reduced_surface
    re = rs.edges[tf.rs_index]
    # BALL's `RSEdge::isFree()` (RSEdge.C:294) returns true iff face_[0]==NULL.
    # If we reach here from the dispatcher, neither slot is 0 — but be
    # defensive: handle free/partially-free by returning so we don't
    # index a NULL/sentinel RS face below.
    (re.f1 == 0 || re.f2 == 0) && return
    # BALL `triangulateNonSingularToricFace` (triangulatedSES.C:402-412)
    # grabs the first 4 SES edges unconditionally. NO check for
    # `face.edges.size() == 4` — if a face has 6 edges (stacked-probe
    # toric), BALL still grabs the first 4 and processes. Match BALL.
    length(tf.edges) >= 4 || return
    e0_idx, e1_idx, e2_idx, e3_idx = tf.edges[1], tf.edges[2], tf.edges[3], tf.edges[4]
    edge1 = tr.ses.edges[e1_idx]
    edge3 = tr.ses.edges[e3_idx]
    # n_seg: BALL uses rsedge.angle * edge3.circle.radius * sqrt(density).
    n_seg = max(1, round(Int, re.angle * edge3.circle.r * tr.sqrt_density, RoundNearestTiesAway))
    psi = re.angle / n_seg
    # Sample edge3 into n_seg+1 points. Use BALL's partitionOfCircle and
    # snap the endpoint.
    p0 = tr.ses.vertices[tf.vertices[1]].point  # first vertex of the normalized face
    # Start the rolling-axis normal from edge3.circle.n. BALL negates it
    # together with reversing edge3_segments when edge3->vertex_[0] != p0
    # (triangulatedSES.C:313-323).
    normal = edge3.circle.n
    edge3_segments = ball_partition_of_circle(edge3.circle,
        tr.ses.vertices[edge3.v1].point, psi, n_seg)
    edge3_segments[end] = tr.ses.vertices[edge3.v2].point
    if tr.ses.vertices[edge3.v1].point != p0
        reverse!(edge3_segments)
        normal = -normal
    end
    edge1_segments = ball_partition_of_circle(edge1.circle,
        tr.ses.vertices[edge1.v1].point, psi, n_seg)
    edge1_segments[end] = tr.ses.vertices[edge1.v2].point
    p1 = tr.ses.vertices[tf.vertices[2]].point
    if tr.ses.vertices[edge1.v1].point != p1
        reverse!(edge1_segments)
    end
    # Build the rolling-axis center circle and sample its n_seg+1 centers.
    center_circle = Circle3{T}(re.center_of_torus, normal,
                                re.radius_of_torus)
    edge0 = tr.ses.edges[e0_idx]
    centers = ball_partition_of_circle(center_circle,
        edge0.circle.p, psi, n_seg)
    centers[end] = tr.ses.edges[e2_idx].circle.p
    # Track the triangle list size so we can apply BALL's
    # post-buildTriangles orientation swap below.
    tri_lo = length(tr.tses.tris) + 1
    ball_build_triangles_regular!(tr, e0_idx, e1_idx, e2_idx, e3_idx,
        centers, edge1_segments, edge3_segments, rs.probe_radius)
    # BALL's `triangulateNonSingularToricFace` (triangulatedSES.C:373-386)
    # tests the first emitted triangle's right-hand normal against
    # centers[0]: if `orth · v0 > orth · centers[0]`, the triangle normal
    # points AWAY from the rolling-axis center (= wrong direction for an
    # outward-pointing toric mesh) and ALL triangles emitted by this
    # buildTriangles call are flipped by swapping vertex 0 and 1.
    # BALL temporarily sets `Constants::EPSILON = 1e-4` during
    # `triangulateToricFaces` (triangulatedSES.C:163-170), so the
    # `isGreater` comparison here uses 1e-4 (not the default 1e-6).
    if tri_lo <= length(tr.tses.tris)
        test_t = tr.tses.tris[tri_lo]
        v0p = tr.tses.points[test_t.v0].point
        v1p = tr.tses.points[test_t.v1].point
        v2p = tr.tses.points[test_t.v2].point
        orth = cross(v1p - v0p, v2p - v0p)
        # BALL `Maths::isGreater(a, b) = (a - b >= EPSILON)`; during the
        # toric phase BALL temporarily sets EPSILON = 1e-4 (triangulatedSES.C:163-170).
        # So the comparison is `(dot(orth, v0p) - dot(orth, centers[0])) >= 1e-4`.
        toric_eps = T(1e-4)
        if dot(orth, v0p) - dot(orth, centers[1]) >= toric_eps
            # BALL `triangulateNonSingularToricFace` (triangulatedSES.C:498-512):
            # swaps vertex_[0]/vertex_[1] only. Edge slots are NOT swapped,
            # which makes the per-triangle edge_[i]↔vertex_[i] correspondence
            # inconsistent after the flip. Match BALL exactly.
            for ti in tri_lo:length(tr.tses.tris)
                t = tr.tses.tris[ti]
                t.v0, t.v1 = t.v1, t.v0
            end
        end
    end
end

# Port of BALL's `SESTriangulator::triangulateFreeToricFace`
# (triangulatedSES.C:390). A "free" toric face arises from an RSEdge
# with both endpoints free (no probe-on-both-sides) — typically an
# isolated 2-atom system where the rolling probe sweeps a full 2π saddle
# around the rolling axis. The face has 2 convex SES edges (full circles
# on atom A and atom B with no corner SES vertices) and the boundary
# wraps around.
function ball_triangulate_free_toric!(tr::SESTriangulator{T},
                                       tf::SESFace{T}) where T
    rs = tr.ses.reduced_surface
    re = rs.edges[tf.rs_index]
    e_front_idx = tf.edges[1]
    e_back_idx  = tf.edges[2]
    e_front = tr.ses.edges[e_front_idx]
    e_back  = tr.ses.edges[e_back_idx]
    # n_points around the full 2π circle. BALL uses circle1.radius
    # (front convex circle).
    n = max(1, round(Int, 2 * T(π) * e_front.circle.r * tr.sqrt_density, RoundNearestTiesAway))
    phi = 2 * T(π) / n
    # Pick a perpendicular direction to `normal`. BALL `triangulateFreeToricFace`
    # (triangulatedSES.C:528-533): orth = (normal.y, -normal.x, 0). If that is
    # the zero vector (= normal is along z), fall back to (normal.z, 0, -normal.x).
    # Match BALL exactly — check for exact zero, not a magnitude threshold.
    normal = e_front.circle.n
    orth = Vector3{T}(normal[2], -normal[1], zero(T))
    if orth == Vector3{T}(0, 0, 0)
        orth = Vector3{T}(normal[3], zero(T), -normal[1])
    end
    ℓ = sqrt(dot(orth, orth))
    orth = ℓ > eps(T) ? orth / ℓ : Vector3{T}(1, 0, 0)
    p1 = e_front.circle.p + orth * e_front.circle.r
    p2 = e_back.circle.p  + orth * e_back.circle.r
    center_circle = Circle3{T}(re.center_of_torus, normal, re.radius_of_torus)
    p3 = re.center_of_torus + orth * re.radius_of_torus
    # Sample each circle into n+1 points (BALL pops_back the last for wrap;
    # we instead use n entries directly and let the cell-builder wrap with
    # modular indexing).
    # IMPORTANT: BALL forces circle2 (back edge) to use the SAME normal as
    # circle1 (front edge) — `TCircle3(back.p, normal=front.n, back.r)` —
    # so the back-arc sampling rotates in the same direction as the front
    # (triangulatedSES.C:396). Our back edge's circle.n may have the
    # opposite sign (since the two convex edges sit on opposite atoms);
    # rebuilding the back circle here with the front normal forces
    # consistent rotation.
    back_circle_forced = Circle3{T}(e_back.circle.p, normal, e_back.circle.r)
    pts1 = ball_partition_of_circle(e_front.circle, p1, phi, n)
    pts2 = ball_partition_of_circle(back_circle_forced, p2, phi, n)
    centers = ball_partition_of_circle(center_circle, p3, phi, n)
    # Drop the duplicate last entry (which equals the first by wrap).
    pop!(pts1); pop!(pts2); pop!(centers)
    ball_build_triangles_free!(tr, e_front_idx, e_back_idx,
        centers, pts1, pts2, rs.probe_radius)
end

# Free-toric grid builder. Like `ball_build_triangles_regular!` but with
# the rolling-axis direction WRAPPING (column n_seg-1 connects back to
# column 0). edge0/edge2 are NULL (no concave boundaries), so vedges at
# the wrap seam are not recorded on any SES edge.
function ball_build_triangles_free!(tr::SESTriangulator{T},
                                     edge1::Int, edge3::Int,
                                     centers::Vector{Vector3{T}},
                                     edge1_segments::Vector{Vector3{T}},
                                     edge3_segments::Vector{Vector3{T}},
                                     probe_radius::T) where T
    n_seg = length(centers)   # wraps: column n_seg ≡ column 0
    # Cross-section refinement from the angle between edge1[0]-centers[0]
    # and edge3[0]-centers[0].
    a = edge1_segments[1] - centers[1]
    b = edge3_segments[1] - centers[1]
    cosθ = dot(a, b) / (sqrt(dot(a, a)) * sqrt(dot(b, b)))
    cosθ = clamp(cosθ, -one(T), one(T))
    ψ = acos(cosθ)
    n_tri = max(1, round(Int, ψ * probe_radius * tr.sqrt_density, RoundNearestTiesAway))
    phi = ψ / n_tri
    # Grid: n_seg × (n_tri+1) — first index wraps mod n_seg.
    grid = Matrix{Int}(undef, n_seg, n_tri + 1)
    for j in 1:n_seg
        u = edge1_segments[j] - centers[j]
        v = edge3_segments[j] - centers[j]
        n_circle = cross(u, v)
        ℓ = sqrt(dot(n_circle, n_circle))
        n_circle = ℓ > eps(T) ? n_circle / ℓ : Vector3{T}(0, 0, 1)
        circle = Circle3{T}(centers[j], n_circle, probe_radius)
        line = ball_partition_of_circle(circle, edge1_segments[j], phi, n_tri)
        for i in 1:(n_tri + 1)
            p = line[i]
            # BALL `triangulatedSES.C:1643` stores `normal_ = centers[i] -
            # line[j]` *unnormalized*; we mirror that so the mesh's normal
            # field exactly matches BALL's TSES output. Renderers that need
            # unit vectors should normalize at consumption time.
            push!(tr.tses.points, TPoint{T}(p, centers[j] - p, 0))
            grid[j, i] = length(tr.tses.points)
        end
    end
    @inline wrap(j) = mod1(j, n_seg)
    # Vertical edges: connect grid[j, i] ↔ grid[j, i+1] (cross-section
    # direction, fixed j). No SES edge recording (edge0/edge2 = NULL).
    vedges = Matrix{Int}(undef, n_seg, n_tri)
    for j in 1:n_seg, i in 1:n_tri
        vedges[j, i] = _ball_new_edge!(tr, grid[j, i], grid[j, i + 1])
    end
    # Horizontal edges: grid[j, i] ↔ grid[j+1 mod n_seg, i] (rolling-axis
    # direction, with wrap). At i=1 the boundary lies on edge1 (front
    # convex), at i=n_tri+1 on edge3 (back convex).
    hedges = Matrix{Int}(undef, n_seg, n_tri + 1)
    for j in 1:n_seg, i in 1:(n_tri + 1)
        se = 0
        i == 1            && (se = edge1)
        i == n_tri + 1    && (se = edge3)
        hedges[j, i] = _ball_new_edge!(tr,
            grid[j, i], grid[wrap(j + 1), i]; ses_edge = se)
    end
    # Diagonal edges: grid[j, i] ↔ grid[j+1 mod n_seg, i+1].
    dedges = Matrix{Int}(undef, n_seg, n_tri)
    for j in 1:n_seg, i in 1:n_tri
        dedges[j, i] = _ball_new_edge!(tr,
            grid[j, i], grid[wrap(j + 1), i + 1])
    end
    # Two triangles per cell (j, i), with wrap on j.
    for j in 1:n_seg, i in 1:n_tri
        j2 = wrap(j + 1)
        _ball_new_tri!(tr,
            grid[j, i], grid[j, i + 1], grid[j2, i + 1],
            vedges[j, i], hedges[j, i + 1], dedges[j, i])
        _ball_new_tri!(tr,
            grid[j, i], grid[j2, i + 1], grid[j2, i],
            dedges[j, i], vedges[j2, i], hedges[j, i])
    end
    nothing
end

# Helpers mirroring BALL's TrianglePoint::has / Triangle::third /
# Triangle::getRelativeIndex.

# Look for an existing TEdge between TPoints v_a and v_b. Returns edge
# index (>0) or 0 if none.
@inline function _ball_find_edge(tr::SESTriangulator{T}, v_a::Int, v_b::Int) where T
    for eid in tr.tses.points[v_a].edges
        e = tr.tses.edges[eid]
        if (e.v0 == v_a && e.v1 == v_b) || (e.v0 == v_b && e.v1 == v_a)
            return eid
        end
    end
    0
end

# Triangle.third(v1, v2): returns the third vertex of a triangle given
# 2 of its vertices.
@inline function _ball_tri_third(t::TTri, v1::Int, v2::Int)
    if (t.v0 == v1) || (t.v0 == v2)
        if (t.v1 == v1) || (t.v1 == v2)
            return t.v2
        else
            return t.v1
        end
    else
        return t.v0
    end
end

# getRelativeIndex(v): which slot (0/1/2 in BALL, 0/1/2 here too) is v in.
# Returns -1 if not found.
@inline function _ball_rel_index(t::TTri, v::Int)
    t.v0 == v && return 0
    t.v1 == v && return 1
    t.v2 == v && return 2
    -1
end

# Port of BALL's `createTriangleAndEdges` (triangulatedSES.C:1250).
# Given a popped edge and a third point, create the new triangle and
# its 2 new edges (reusing existing edges where possible). The geometric
# outward-orientation swap at the end ensures the triangle's vertices
# are in CCW outward order on the sphere.
function _ball_create_triangle_and_edges!(tr::SESTriangulator{T},
                                           edge_idx::Int, pid::Int,
                                           sphere_p::Vector3{T},
                                           convex::Bool) where T
    e = tr.tses.edges[edge_idx]
    # BALL `createTriangleAndEdges` (triangulatedSES.C:1416-1480):
    # edge1->vertex_[0] = edge->vertex_[0]; edge1->vertex_[1] = point
    # edge2->vertex_[0] = point;            edge2->vertex_[1] = edge->vertex_[1]
    # ... triangle vertex_[0] = edge->vertex_[1], vertex_[1] = edge->vertex_[0],
    # vertex_[2] = point. Then conditional swap of vertex_[0]/vertex_[1] only.
    # Match BALL's exact edge1/edge2 endpoints so subsequent border-push
    # order in `_ball_link_new_tri!` and edge identity in `_ball_find_edge`
    # match BALL's iteration.
    edge1_id = _ball_find_edge(tr, e.v0, pid)   # BALL edge1: (edge.v[0], point)
    old1 = edge1_id != 0
    if !old1
        push!(tr.tses.edges, TEdge(e.v0, pid, 0, 0))
        edge1_id = length(tr.tses.edges)
    end
    edge2_id = _ball_find_edge(tr, pid, e.v1)   # BALL edge2: (point, edge.v[1])
    old2 = edge2_id != 0
    if !old2
        push!(tr.tses.edges, TEdge(pid, e.v1, 0, 0))
        edge2_id = length(tr.tses.edges)
    end
    # Triangle layout: (e.v1, e.v0, point) initially; edges (orig, edge1, edge2).
    # Edge convention: edge_[i] between vertex_[i] and vertex_[(i+1)%3].
    v0_final = e.v1
    v1_final = e.v0
    v2 = pid
    p0 = tr.tses.points[v0_final].point
    p1 = tr.tses.points[v1_final].point
    p2 = tr.tses.points[v2].point
    test_vec = cross(p1 - p0, p2 - p0)
    test_val = dot(test_vec, sphere_p - p0)
    # BALL `createTriangleAndEdges` (triangulatedSES.C:1473-1474) uses
    # `Maths::isGreater(test_value, 0)` / `Maths::isLess(test_value, 0)`,
    # which compare against `+EPSILON` / `-EPSILON` respectively
    # (1e-6 here), not against 0 strictly.
    EPS_ORIENT = T(1e-6)
    if (convex && test_val >= EPS_ORIENT) || (!convex && -test_val >= EPS_ORIENT)
        # BALL swaps vertex_[0]/vertex_[1] only — edges stay in place.
        v0_final, v1_final = v1_final, v0_final
    end
    if haskey(ENV, "SES_DUP_DETECT")
        # Sanity check: does the new triangle have the SAME 3 vertices as any
        # existing one? If so, the BST is creating a true duplicate.
        key = sort([v0_final, v1_final, v2])
        ctx = get(ENV, "SES_TRI_CTX", "?")
        for (j, ot) in enumerate(tr.tses.tris)
            if sort([ot.v0, ot.v1, ot.v2]) == key
                println(stderr, "DUP_TRI ctx=$ctx new=(($v0_final,$v1_final,$v2)) matches existing tri $j=(($(ot.v0),$(ot.v1),$(ot.v2)))")
                break
            end
        end
    end
    push!(tr.tses.tris,
          TTri(v0_final, v1_final, v2, edge_idx, edge1_id, edge2_id))
    tid = length(tr.tses.tris)
    return tid, edge1_id, old1, edge2_id, old2
end

# Port of BALL's `buildSphericTriangles` (triangulatedSES.C:777).
# Runs the advancing-front on a spheric/contact face.
#
# Arguments:
#   * `face` — the SES face being triangulated.
#   * `sphere_p` — center of the patch sphere (atom centre for contact,
#     probe centre for spheric). Used by the geometric outward swap in
#     `_ball_create_triangle_and_edges!`.
#   * `convex` — true for contact faces (atom surface, outward = away
#     from sphere_p), false for spheric faces (probe surface, outward
#     = toward sphere_p).
#   * `extra_points` — additional candidate TPoint indices to consider
#     (= the icosphere interior points for contact faces; empty for
#     spheric faces).
#
# Walks the face's SES edges, collecting their per-edge TEdges from
# `tr.edge[se_idx]` (populated by toric face triangulation) onto the
# border. For each border edge, runs BALL's iterative quickhull on
# the candidate set and emits a new triangle.
function ball_build_spheric_triangles!(tr::SESTriangulator{T},
                                        face::SESFace{T},
                                        sphere_p::Vector3{T},
                                        convex::Bool;
                                        extra_points::Vector{Int} = Int[]) where T
    # 1. Build the candidate set: extra_points (icosphere interior) +
    #    all TPoints incident to TEdges sampling this face's SES edges.
    #    Mirrors BALL's `buildSphericTriangles` (triangulatedSES.C:777-796).
    candidates = copy(extra_points)
    seen_pts = Set{Int}(extra_points)
    for se_idx in face.edges
        se_idx == 0 && continue
        for te_idx in tr.edge[se_idx]
            te = tr.tses.edges[te_idx]
            te.v0 in seen_pts || (push!(candidates, te.v0); push!(seen_pts, te.v0))
            te.v1 in seen_pts || (push!(candidates, te.v1); push!(seen_pts, te.v1))
        end
    end
    # Deterministic iteration order in candidates can lead to order-sensitive
    # behavior in the BST main loop's "best" selection (each new "best" updates
    # the comparison plane normal). BALL's `HashSet<TrianglePoint*>` iterates
    # in pointer-hash order — random across runs but deterministic for one
    # invocation. Sort by TPoint index to give a stable, well-defined order
    # that exercises the BST consistently. Diagnostic only — toggle via env.
    haskey(ENV, "SES_BST_SORT_CANDS") && sort!(candidates)
    # 2. The border starts EMPTY (BALL triangulatedSES.C:797). It gets
    #    populated by `buildFirstTriangle` which seeds the advancing
    #    front from a single triangle on a non-singular SES edge. The
    #    main loop then expands this single front outward until it
    #    naturally reaches all the other boundary SES edges (whose
    #    toric-side TEdges are already in `tr.edge[se_idx]` and get
    #    attached on the second pass via `_ball_find_edge`).
    border = Int[]
    seeded = ball_build_first_triangle!(tr, face, border, candidates, sphere_p, convex)
    seeded || return false
    EPS = T(1e-6)
    iter = 0
    do_trace = haskey(ENV, "SES_TRI_DEEP") && get(ENV, "SES_TRI_CTX", "") == get(ENV, "SES_TRI_DEEP", "")
    do_trace && println("DEEP: ctx=$(ENV["SES_TRI_CTX"]) candidates=$candidates border_init=$border")
    # BALL `buildSphericTriangles` (triangulatedSES.C:941) just iterates
    # `while (!border.empty())` with no max-iter cutoff. Match BALL.
    while !isempty(border)
        iter += 1
        eid = popfirst!(border)
        edge = tr.tses.edges[eid]
        # BALL's BST main loop (triangulatedSES.C:940-944) does NOT skip
        # fully-saturated edges. Edges with both face_[0] and face_[1] set
        # should have been removed from border by the linking step
        # (border.remove). If they weren't, BALL processes them anyway and
        # silently emits an extra triangle. Match BALL.
        start_tri = tr.tses.tris[edge.f0]
        third_pid = _ball_tri_third(start_tri, edge.v0, edge.v1)
        if do_trace
            ep0 = tr.tses.points[edge.v0].point
            ep1 = tr.tses.points[edge.v1].point
            tp = tr.tses.points[third_pid].point
            println("DEEP iter $iter: eid=$eid v0=$(edge.v0)$(round.(Tuple(ep0); digits=4)) v1=$(edge.v1)$(round.(Tuple(ep1); digits=4)) third=$third_pid$(round.(Tuple(tp); digits=4)) border_size=$(length(border))")
        end
        # Orientation swap (BALL triangulatedSES.C:812-821).
        p0_idx = _ball_rel_index(start_tri, edge.v0)
        p1_idx = _ball_rel_index(start_tri, edge.v1)
        diff = p1_idx - p0_idx
        if (convex && (diff == -1 || diff == 2)) ||
           (!convex && (diff == 1 || diff == -2))
            edge.v0, edge.v1 = edge.v1, edge.v0
        end
        # BALL's BST inner loop (triangulatedSES.C:963-1006) has NO
        # saturation pre-filter. It picks the best candidate purely on the
        # `this_value > test_value` criterion and lets `buildUnambiguousTriangle`
        # / `buildAmbiguousTriangles` deal with the linking. Match BALL.
        best = Int[]
        best_normal = Vector3{T}(0, 0, 0)
        test_value = zero(T)
        have_seed = false
        for pid in candidates
            (pid == edge.v0 || pid == edge.v1 || pid == third_pid) && continue
            pp = tr.tses.points[pid].point
            ep0 = tr.tses.points[edge.v0].point
            ep1 = tr.tses.points[edge.v1].point
            n = cross(pp - ep1, pp - ep0)
            if !have_seed
                push!(best, pid)
                best_normal = n
                test_value = dot(n, ep0)
                have_seed = true
                continue
            end
            this_value = dot(best_normal, pp)
            # BALL `Maths::isGreater(a, b)` = `(a - b >= EPSILON)` (note `>=`,
            # not `>`); `Maths::isEqual(a, b)` = `(abs(a - b) < EPSILON)`
            # (strict `<`). Match exactly.
            if this_value - test_value >= EPS
                empty!(best)
                push!(best, pid)
                best_normal = n
                test_value = dot(best_normal, ep0)
            elseif abs(this_value - test_value) < EPS
                push!(best, pid)
            end
        end
        if do_trace
            tids = [(pid, round.(Tuple(tr.tses.points[pid].point); digits=4)) for pid in best]
            println("DEEP iter $iter: best (count=$(length(best))) = $tids")
        end
        if length(best) == 1
            _ball_build_unambiguous!(tr, border, eid, best[1], sphere_p, convex)
        elseif length(best) > 1
            _ball_build_ambiguous!(tr, border, eid, best, sphere_p, convex)
        end
    end
    return true
end

# Port of BALL's `buildFirstTriangle` (triangulatedSES.C:882). Picks
# the seed triangle for the advancing-front using the strict full-hull
# check: a candidate `p` is valid iff NO other candidate is strictly on
# the wrong side of the (edge, p) plane. Ties tie-break by minimum
# oriented angle to the existing toric triangle's normal.
#
# After this seeds 1-2 triangles, the main loop's iterative criterion
# takes over.
function ball_build_first_triangle!(tr::SESTriangulator{T},
                                     face::SESFace{T},
                                     border::Vector{Int},
                                     candidates::Vector{Int},
                                     sphere_p::Vector3{T},
                                     convex::Bool) where T
    # Find the first valid seed TriangleEdge — mirrors BALL's
    # `firstSESEdge` (triangulatedSES.C:996). Skip singular SES edges;
    # if a non-singular SES edge has exactly 1 TriangleEdge, skip it
    # when that TriangleEdge is degenerate (BALL: `diff*diff < 0.01`).
    # We also require f0 != 0 (BALL implicitly assumes face_[0] is set
    # because the seed comes from a toric-sampled SES edge).
    seed_eid = 0
    for se_idx in face.edges
        se_idx == 0 && continue
        te_list = tr.edge[se_idx]
        isempty(te_list) && continue
        se = tr.ses.edges[se_idx]
        se.type == SESEdgeType.Singular && continue
        # BALL's degenerate check is on the FIRST TriangleEdge when
        # only 1 exists.
        if length(te_list) == 1
            te = tr.tses.edges[te_list[1]]
            d = tr.tses.points[te.v0].point - tr.tses.points[te.v1].point
            dot(d, d) < T(0.01) && continue
        end
        # BALL's `firstSESEdge` (triangulatedSES.C:1162) returns the SES
        # edge here unconditionally once the type+degeneracy filters pass.
        # No additional f0-set check — BALL assumes the SES edge already
        # has a toric-side triangle attached because firstSESEdge is only
        # called from buildFirstTriangle (after toric phase).
        seed_eid = te_list[1]
        break
    end
    seed_eid == 0 && return false
    edge = tr.tses.edges[seed_eid]
    # Orient the edge so its traversal in the existing triangle is
    # FORWARD (= the new triangle traverses it BACKWARD). BALL does
    # this unconditionally (triangulatedSES.C:906-917).
    start_tri = tr.tses.tris[edge.f0]
    p0_idx = _ball_rel_index(start_tri, edge.v0)
    p1_idx = _ball_rel_index(start_tri, edge.v1)
    diff = p1_idx - p0_idx
    if (diff == 1) || (diff == -2)
        edge.v0, edge.v1 = edge.v1, edge.v0
    end
    p0 = tr.tses.points[edge.v0].point
    p1 = tr.tses.points[edge.v1].point
    edge_vec = p0 - p1
    # `same_edge` = ALL TPoints sampling the seed SES edge (BALL
    # triangulatedSES.C:896-903). These cannot be valid third vertices
    # because they're collinear with the seed edge.
    seed_se_idx = 0
    for se_idx in face.edges
        se_idx == 0 && continue
        if seed_eid in tr.edge[se_idx]
            seed_se_idx = se_idx
            break
        end
    end
    same_edge = Set{Int}()
    if seed_se_idx != 0
        for te_idx in tr.edge[seed_se_idx]
            te = tr.tses.edges[te_idx]
            push!(same_edge, te.v0)
            push!(same_edge, te.v1)
        end
    else
        push!(same_edge, edge.v0)
        push!(same_edge, edge.v1)
    end
    # Strict hull check (BALL triangulatedSES.C:918-951). The OUTER loop
    # excludes candidates in `same_edge`. The INNER loop excludes only
    # the candidate itself and the two edge endpoints — NOT all
    # `same_edge` members.
    # The strict comparison uses BALL's `Maths::isGreater` with default
    # epsilon (1e-6 here): a candidate fails the hull test only when
    # ANOTHER candidate is more than `eps` above its plane.
    third = Int[]
    EPS = T(1e-6)
    for pid in candidates
        pid in same_edge && continue
        pp = tr.tses.points[pid].point
        n = cross(edge_vec, p0 - pp)
        tv = dot(n, p0)
        is_hull = true
        for tid in candidates
            (tid == pid || tid == edge.v0 || tid == edge.v1) && continue
            tp = tr.tses.points[tid].point
            this_value = dot(n, tp)
            if (convex && this_value - tv >= EPS) ||
               (!convex && tv - this_value >= EPS)
                is_hull = false
                break
            end
        end
        is_hull && push!(third, pid)
    end
    # Tie-break by minimum oriented angle to the existing triangle's
    # outward normal (BALL triangulatedSES.C:957-979).
    real_third = if length(third) <= 1
        third
    else
        start_tri = tr.tses.tris[edge.f0]
        sp0 = tr.tses.points[start_tri.v0].point
        sp1 = tr.tses.points[start_tri.v1].point
        sp2 = tr.tses.points[start_tri.v2].point
        ref_normal = cross(sp0 - sp1, sp0 - sp2)
        best_ang = T(3 * π)
        out = Int[]
        for pid in third
            pp = tr.tses.points[pid].point
            new_normal = cross(edge_vec, p0 - pp)
            ang = oriented_angle(ref_normal, new_normal, edge_vec)
            if ang <= best_ang
                if ang < best_ang
                    empty!(out)
                    best_ang = ang
                end
                push!(out, pid)
            end
        end
        out
    end
    # BALL `buildFirstTriangle` (triangulatedSES.C:1148-1157): no
    # border-removal of the seed — the seed comes from the SES edge's
    # TEdge list, not from `border` (which starts empty in
    # buildSphericTriangles).
    if length(real_third) == 1
        _ball_build_unambiguous!(tr, border, seed_eid, real_third[1],
                                  sphere_p, convex)
    elseif length(real_third) > 1
        _ball_build_ambiguous!(tr, border, seed_eid, real_third,
                                sphere_p, convex)
    end
    # BALL leaves real_third empty silently when the strict-hull check
    # fails — `case 0: break` in BALL's switch (triangulatedSES.C:984).
    # No fallback; the face just produces fewer triangles in that case.
    true
end

# Mirrors BALL's `buildUnambiguousTriangle` (triangulatedSES.C:1041) +
# the linking step.
function _ball_build_unambiguous!(tr::SESTriangulator{T},
                                   border::Vector{Int},
                                   eid::Int, pid::Int,
                                   sphere_p::Vector3{T}, convex::Bool) where T
    tid, eid1, old1, eid2, old2 = _ball_create_triangle_and_edges!(
        tr, eid, pid, sphere_p, convex)
    tid == 0 && return 0
    _ball_link_new_tri!(tr, border, eid, tid, eid1, old1, eid2, old2)
    tid
end

# Mirrors BALL's linking step (the post-createTriangleAndEdges block
# in `buildUnambiguousTriangle`, triangulatedSES.C:1057-1098). For each
# of the 2 new edges:
#   * If old: write to the edge's `f0`/`f1` slot (whichever is open),
#     and drop the edge from `border` (it's now fully attached).
#   * If new: link to its endpoints, set `f0 = tid`, push on `border`.
# Also fills the popped edge's `f1 = tid`. Mirrors BALL's
# `edge->face_[1] = triangle;` at line 1099.
function _ball_link_new_tri!(tr::SESTriangulator{T},
                              border::Vector{Int}, eid::Int, tid::Int,
                              eid1::Int, old1::Bool,
                              eid2::Int, old2::Bool) where T
    tri = tr.tses.tris[tid]
    # BALL `buildUnambiguousTriangle` (triangulatedSES.C:1265-1268) and
    # `buildAmbiguousTriangles` (triangulatedSES.C:1343-1346) both set
    # `edge0->face_[1] = triangle` FIRST, then insert the triangle into
    # vertex face sets. Match the order even though the final state is
    # the same.
    tr.tses.edges[eid].f1 = tid
    push!(tr.tses.points[tri.v0].faces, tid)
    push!(tr.tses.points[tri.v1].faces, tid)
    push!(tr.tses.points[tri.v2].faces, tid)
    for (eid_new, old) in ((eid1, old1), (eid2, old2))
        e = tr.tses.edges[eid_new]
        if old
            # BALL `buildUnambiguousTriangle` (triangulatedSES.C:1058-1098):
            # `if face_[0] == NULL → face_[0] = triangle; else → face_[1] = triangle`.
            # NO check whether face_[1] was already set — BALL overwrites. Match
            # this exactly even though it may corrupt manifold topology in
            # pathological iteration orders.
            if e.f0 == 0
                e.f0 = tid
            else
                e.f1 = tid
            end
            # BALL `border.remove(edge1)` (triangulatedSES.C:1242, 1264):
            # std::list::remove removes ALL occurrences. If an edge was
            # somehow pushed to `border` more than once, BALL still purges
            # all copies. Match by filtering rather than findfirst+deleteat.
            filter!(!=(eid_new), border)
        else
            push!(tr.tses.points[e.v0].edges, eid_new)
            push!(tr.tses.points[e.v1].edges, eid_new)
            e.f0 = tid
            push!(border, eid_new)
        end
    end
end

# Mirrors BALL's `buildAmbiguousTriangles` (triangulatedSES.C:1108).
function _ball_build_ambiguous!(tr::SESTriangulator{T},
                                 border::Vector{Int}, seed_eid::Int,
                                 pids::Vector{Int}, sphere_p::Vector3{T},
                                 convex::Bool) where T
    seed_edge = tr.tses.edges[seed_eid]
    cands = copy(pids)
    push!(cands, seed_edge.v0)
    push!(cands, seed_edge.v1)
    planar = Int[seed_eid]
    do_trace = haskey(ENV, "SES_TRI_DEEP") && get(ENV, "SES_TRI_CTX", "") == get(ENV, "SES_TRI_DEEP", "")
    do_trace && println(stderr, "    AMBIG: seed_eid=$seed_eid pids=$pids cands=$cands")
    amb_iter = 0
    while !isempty(planar)
        amb_iter += 1
        e0_id = popfirst!(planar)
        e0 = tr.tses.edges[e0_id]
        built = false
        for pid in cands
            (pid == e0.v0 || pid == e0.v1) && continue
            tid, eid1, old1, eid2, old2 = _ball_create_triangle_and_edges!(
                tr, e0_id, pid, sphere_p, convex)
            keep = (e0_id == seed_eid)
            if !keep
                # BALL `buildAmbiguousTriangles` (triangulatedSES.C:1316):
                # accesses `edge0->face_[0]` unconditionally — no
                # `e0.f0 > 0` defensive check. BALL assumes the popped
                # planar edge always has face_[0] set (it must, since it
                # got onto planar by being either the seed or a newly-
                # linked edge of a kept triangle).
                old_tri = tr.tses.tris[e0.f0]
                i1 = _ball_rel_index(old_tri, e0.v0)
                i2 = _ball_rel_index(old_tri, e0.v1)
                back = (i1 - i2 == 1) || (i1 - i2 == -2)
                new_tri = tr.tses.tris[tid]
                i1 = _ball_rel_index(new_tri, e0.v0)
                i2 = _ball_rel_index(new_tri, e0.v1)
                if back
                    keep = (i1 - i2 == -1) || (i1 - i2 == 2)
                else
                    keep = (i1 - i2 == 1) || (i1 - i2 == -2)
                end
            end
            if keep
                _ball_link_new_tri!(tr, border, e0_id, tid,
                                     eid1, old1, eid2, old2)
                # BALL `planar_edges.remove(edge)` (triangulatedSES.C:1192,
                # 1215): std::list::remove removes ALL occurrences. Match
                # via filter!.
                if old1
                    filter!(!=(eid1), planar)
                else
                    push!(planar, eid1)
                end
                if old2
                    filter!(!=(eid2), planar)
                else
                    push!(planar, eid2)
                end
                # BALL `border.remove(edge0)` (triangulatedSES.C:1227): same
                # semantic — remove all copies. Without this, the main BST
                # loop will re-pop the same edge later.
                filter!(!=(e0_id), border)
                if do_trace
                    println(stderr, "    AMBIG iter $amb_iter: e0_id=$e0_id pid=$pid tid=$tid eid1=$eid1(old=$old1) eid2=$eid2(old=$old2) planar_size=$(length(planar))")
                end
                built = true
                break
            else
                # Roll back: pop the freshly created triangle and edges.
                pop!(tr.tses.tris)
                !old2 && pop!(tr.tses.edges)
                !old1 && pop!(tr.tses.edges)
            end
        end
        built || nothing
    end
end

# Port of BALL's `triangulateContactFace` (triangulatedSES.C:676).
# For a contact face on atom `face.sphere`:
#   1. Take the icosphere template at the appropriate refinement level.
#   2. Scale it by the atom radius.
#   3. Cut it against each SES edge's carrier-circle plane (keeping
#      only points on the contact-patch side, with fuzzy = 0.05).
#   4. Shift to atom centre and register surviving icosphere points as
#      new TPoints.
#   5. Run `ball_build_spheric_triangles!` with `convex=true`, passing
#      the icosphere TPoints as `extra_points` candidates.
function ball_triangulate_contact_face!(tr::SESTriangulator{T},
                                         face::SESFace{T}) where T
    sphere_p = Vector3{T}(face.sphere.center...)
    sphere_r = face.sphere.r
    # Diagnostic: dump boundary TPoint inventory (mirrors BALL's dump in
    # triangulatedSES.C:280-291). Gated on env SES_CF_INV matching the
    # 0-based atom index (BALL convention).
    rs_atom_1based = face.rs_index == 0 ? 0 :
                     tr.ses.reduced_surface.vertices[face.rs_index].atom
    if haskey(ENV, "SES_CF_INV") &&
       tryparse(Int, ENV["SES_CF_INV"]) == rs_atom_1based - 1
        println(stderr, "[CF_INV] atom_idx=$(rs_atom_1based-1) radius=$(face.sphere.r) center=$(round.(Tuple(sphere_p); digits=6))")
        println(stderr, "  SES edges ($(length(face.edges))):")
        for (i, se_idx) in enumerate(face.edges)
            if se_idx == 0
                println(stderr, "    [$(i-1)] SENTINEL")
                continue
            end
            se = tr.ses.edges[se_idx]
            n_te = length(tr.edge[se_idx])
            println(stderr, "    [$(i-1)] se_idx=$(se_idx-1) type=$(Int(se.type)) tedges=$n_te")
        end
        bset = Set{Int}()
        for se_idx in face.edges
            se_idx == 0 && continue
            for te_idx in tr.edge[se_idx]
                te = tr.tses.edges[te_idx]
                push!(bset, te.v0)
                push!(bset, te.v1)
            end
        end
        println(stderr, "  Boundary TPoint inventory ($(length(bset)) unique):")
        for pid in bset
            tp = tr.tses.points[pid]
            println(stderr, "    tp_idx=$(pid-1) pos=$(round.(Tuple(tp.point); digits=4)) n_edges=$(length(tp.edges))")
        end
    end
    # Free contact face (no boundary edges): emit the full icosphere
    # scaled and shifted. BALL emits a `TriangulatedSphere` directly.
    if isempty(face.edges)
        n_ref = ball_num_refinements(tr.density, sphere_r)
        tmpl = _ball_template_at(tr, n_ref)
        # Add each template vertex as a new TPoint, scale and shift.
        first_pid = length(tr.tses.points) + 1
        for v in tmpl.position
            outward = Vector3{T}(v...)
            p = sphere_p + sphere_r * outward
            push!(tr.tses.points, TPoint{T}(p, outward, 0))
        end
        # Emit the icosphere triangles directly (no advancing-front needed).
        # BALL `triangulateContactFace` (triangulatedSES.C:769-781) emits
        # template triangles without a degenerate guard — the icosphere
        # template never has coincident vertices.
        for f in tmpl.faces
            a = first_pid + f[1] - 1
            b = first_pid + f[2] - 1
            c = first_pid + f[3] - 1
            push!(tr.tses.tris, TTri(a, b, c, 0, 0, 0))
            tid = length(tr.tses.tris)
            push!(tr.tses.points[a].faces, tid)
            push!(tr.tses.points[b].faces, tid)
            push!(tr.tses.points[c].faces, tid)
        end
        return
    end
    # BALL special case (triangulatedSES.C:682-687): if the contact face
    # has exactly 2 SES edges and each has only 1 TEdge, BALL early-returns
    # without emitting any triangles. Mirror that here.
    if length(face.edges) == 2 &&
       length(tr.edge[face.edges[1]]) == 1 &&
       length(tr.edge[face.edges[2]]) == 1
        return
    end
    # BALL special case (triangulatedSES.C:688-716): if the contact face
    # has 3 SES edges each with 1 TEdge, emit ONE triangle from the 4
    # endpoint TPoints (= 3 unique corners). Orient via the atom-center
    # outward test.
    if length(face.edges) == 3 &&
       length(tr.edge[face.edges[1]]) == 1 &&
       length(tr.edge[face.edges[2]]) == 1 &&
       length(tr.edge[face.edges[3]]) == 1
        # BALL `triangulateContactFace` (triangulatedSES.C:814-842): insert
        # endpoints of front and back SES edges' TEdges (4 inserts → 3 unique
        # in a triangular face), take the first 3 from the set, build a
        # single triangle. The middle SES edge's endpoints aren't explicitly
        # inserted — they coincide with one of the front/back endpoints in
        # a 3-edge closed cycle.
        ept = Set{Int}()
        te_front = tr.tses.edges[tr.edge[face.edges[1]][1]]
        te_back  = tr.tses.edges[tr.edge[face.edges[3]][1]]
        push!(ept, te_front.v0); push!(ept, te_front.v1)
        push!(ept, te_back.v0);  push!(ept, te_back.v1)
        # BALL takes the first 3 from the HashSet iteration regardless of
        # count (`for (i = 0; i < 3; i++) p++;` — undefined behaviour if
        # fewer). For safety, take min(3, length).
        vs = collect(ept)
        if length(vs) >= 3
            vs = vs[1:3]
            v0p = tr.tses.points[vs[1]].point
            v1p = tr.tses.points[vs[2]].point
            v2p = tr.tses.points[vs[3]].point
            normal = cross(v0p - v1p, v0p - v2p)
            # BALL `triangulateContactFace` (triangulatedSES.C:833):
            # `if (Maths::isGreater(normal * (sphere.p - v1), 0.0))`. The
            # `isGreater` boundary is `>= EPSILON` (1e-6), not strictly > 0.
            if dot(normal, sphere_p - v1p) >= T(1e-6)
                vs[1], vs[2] = vs[2], vs[1]
            end
            push!(tr.tses.tris, TTri(vs[1], vs[2], vs[3], 0, 0, 0))
            tid = length(tr.tses.tris)
            push!(tr.tses.points[vs[1]].faces, tid)
            push!(tr.tses.points[vs[2]].faces, tid)
            push!(tr.tses.points[vs[3]].faces, tid)
            return
        end
    end
    # Normal case (BALL triangulatedSES.C:718-745): pick refinement
    # level from `numberOfRefinements(density, radius)`, cut the icosphere
    # against each SES edge's circle plane (fuzzy 0.05), shift to atom
    # center, then run `buildSphericTriangles` with the surviving interior
    # points as candidates.
    n_ref = ball_num_refinements(tr.density, sphere_r)
    tmpl = _ball_template_at(tr, n_ref)
    fuzzy = T(0.05)
    rs = tr.ses.reduced_surface
    extra_points = Int[]
    for u_t in tmpl.position
        outward = Vector3{T}(u_t...)
        point = sphere_p + sphere_r * outward
        keep = true
        for se_idx in face.edges
            se_idx == 0 && continue
            se = tr.ses.edges[se_idx]
            # BALL `TriangulatedSurface::cut(plane, 0.05)` uses
            # `plane.n * (point - plane.p) <= fuzzy` with an
            # UNNORMALIZED `plane.n` set to `circle0.p - circle1.p`
            # (createConvexEdge, solventExcludedSurface.C:1346).
            # My SES convex edges store the NORMALISED normal, so
            # I must reconstruct the unnormalised magnitude (= the
            # distance between the two contact circles) to keep the
            # clipping equally aggressive.
            ℓn = sqrt(dot(se.circle.n, se.circle.n))
            ℓn < eps(T) && continue
            n_hat = se.circle.n / ℓn
            # Recover BALL's unnormalised normal magnitude from the
            # underlying RS edge's contact circles (se.f2 == RS edge
            # index for convex edges).
            mag = one(T)
            if se.type == SESEdgeType.Convex && se.f2 != 0 && se.f2 <= length(rs.edges)
                re = rs.edges[se.f2]
                if re.v1 != 0
                    dv = re.contact_circle1.p - re.contact_circle2.p
                    mag = sqrt(dot(dv, dv))
                    mag = mag > eps(T) ? mag : one(T)
                end
            end
            # `dot(point - p, n_hat) * mag <= fuzzy` ⟺
            # `dot(point - p, n_hat) <= fuzzy / mag`
            if dot(point - se.circle.p, n_hat) * mag <= fuzzy
                keep = false; break
            end
        end
        keep || continue
        push!(tr.tses.points, TPoint{T}(point, outward, 0))
        push!(extra_points, length(tr.tses.points))
    end
    ball_build_spheric_triangles!(tr, face, sphere_p, true;
                                   extra_points = extra_points)
end

# Port of BALL's `triangulateSphericFace` (triangulatedSES.C:761) —
# thin wrapper around `ball_build_spheric_triangles!` with `convex=false`
# and no icosphere (spheric faces are bounded only by their concave
# arcs; no interior candidates are needed beyond the boundary).
function ball_triangulate_spheric_face!(tr::SESTriangulator{T},
                                         face::SESFace{T}) where T
    sphere_p = Vector3{T}(face.sphere.center...)
    # BALL `triangulateSphericFace` (triangulatedSES.C:887-900) builds into
    # a LOCAL `TriangulatedSES part`. On failure, `part` is discarded;
    # only on success does it call `tses_->join(part)`. My port builds
    # directly into `tr.tses`, so a failed attempt leaves partial
    # triangles/edges/TPoints behind that the spheric_retry then DUPLICATES.
    # Snapshot the sizes; on failure, roll back to the snapshot — matching
    # BALL's "discard part on failure" semantics.
    n_tris_before  = length(tr.tses.tris)
    n_edges_before = length(tr.tses.edges)
    n_pts_before   = length(tr.tses.points)
    ok = ball_build_spheric_triangles!(tr, face, sphere_p, false)
    if !ok
        # Roll back any partial triangles/edges/TPoints created by the
        # failed buildSphericTriangles attempt. Vertex.faces and edge
        # face_[0]/face_[1] mutations on PRE-EXISTING entities (e.g.
        # toric-side TEdges) also need to be undone, but those are
        # harder to track — for now, just trim the appended tails.
        # In practice the failed attempts don't update existing TEdges
        # (because the first triangle's seed wasn't built), so trimming
        # is sufficient.
        if length(tr.tses.tris) > n_tris_before
            resize!(tr.tses.tris, n_tris_before)
        end
        if length(tr.tses.edges) > n_edges_before
            resize!(tr.tses.edges, n_edges_before)
        end
        if length(tr.tses.points) > n_pts_before
            resize!(tr.tses.points, n_pts_before)
        end
    end
    return ok
end

# Construct a fresh SESTriangulator and run preProcessing.
function SESTriangulator{T}(ses::SolventExcludedSurface{T}, density::T) where T
    # BALL `buildTemplateSpheres` (triangulatedSES.C) pre-builds icosphere
    # templates up to a fixed level; my port extends on demand so the cap
    # is whatever `numberOfRefinements` returns for the input. Pre-seed
    # the 4 smallest since they cover typical density+radius inputs.
    templates = Any[icosphere(T, 0), icosphere(T, 1), icosphere(T, 2), icosphere(T, 3)]
    tr = SESTriangulator{T}(ses, density, sqrt(density),
                             TSES{T}(), Int[], Vector{Vector{Int}}(),
                             templates)
    ball_pre_processing!(tr)
    tr
end

# Extract the final mesh from the SESTriangulator. Mirrors BALL's
# join-into-tses_ then setIndices flow. Each TPoint becomes a
# GeometryBasics vertex (with normal); each TTri becomes a face. The
# triangle winding is oriented per `_push_triangle!` outward-hint
# logic so the final mesh has CCW-from-outside winding regardless of
# the internal TTri vertex order.
function ball_extract_mesh(tr::SESTriangulator{T}) where T
    n_pts = length(tr.tses.points)
    n_tris = length(tr.tses.tris)
    # BALL `TriangulatedSurface::deleteIsolatedPoints` (triangulatedSurface.C:736)
    # is called inside each face's triangulation (per `part` in
    # triangulatedSES.C:872) before joining into the global mesh. The end
    # result is the global mesh has no TPoint with zero incident triangles.
    # My port works directly on the global `tr.tses` and accumulates
    # all icosphere-template TPoints (including clipped-but-unused ones).
    # Apply the same removal globally here.
    used = falses(n_pts)
    for tt in tr.tses.tris
        used[tt.v0] = true
        used[tt.v1] = true
        used[tt.v2] = true
    end
    # Build old → new index remap; emit only used TPoints.
    remap = zeros(Int, n_pts)
    pos = Point3{T}[]
    norm = Vec3{T}[]
    sizehint!(pos,  count(used))
    sizehint!(norm, count(used))
    for i in 1:n_pts
        used[i] || continue
        tp = tr.tses.points[i]
        push!(pos,  Point3{T}(tp.point...))
        push!(norm, Vec3{T}(tp.normal...))
        remap[i] = length(pos)
    end
    faces = Vector{TriangleFace{Int}}(undef, n_tris)
    for i in 1:n_tris
        tt = tr.tses.tris[i]
        # Vertex order is already correct from `createTriangleAndEdges`'
        # geometric swap and the post-buildTriangles orientation checks
        # for toric faces. BALL emits triangles in their construction
        # order — we do the same.
        faces[i] = TriangleFace{Int}(remap[tt.v0], remap[tt.v1], remap[tt.v2])
    end
    GeometryBasics.Mesh(pos, faces; normal = norm)
end

"""
    triangulate_ses_ball(ses::SolventExcludedSurface; density)

BALL-faithful SES triangulator (port of `SESTriangulator::run`,
triangulatedSES.C:130). Runs:

  1. preProcessing — creates one TPoint per SES vertex.
  2. toric faces — non-singular grid + singular fan; free variant
     still pending.
  3. contact faces — cut icosphere + `buildSphericTriangles`.
  4. spheric faces — `buildSphericTriangles` only.

Returns a `GeometryBasics.Mesh`. Free toric faces are currently
skipped; for inputs without isolated atoms (most realistic cases)
the output is the full BALL-faithful triangulation.
"""

function triangulate_ses_ball(ses::SolventExcludedSurface{T};
                              density::T = T(1.0)) where T
    # BALL `SESTriangulator::preProcessing` (triangulatedSES.C:188-189):
    #   tses_->ses_->clean(tses_->density_);
    #   tses_->ses_->splitSphericFaces();
    # Mirror this order: clean small toric faces first using the user
    # density, then split disconnected spheric face boundaries.
    _ses_clean_small_toric_faces!(ses; density = density)
    split_spheric_faces!(ses)
    tr = SESTriangulator{T}(ses, density)
    rs = ses.reduced_surface
    # Diagnostic: record source of each triangle. Enabled via env var.
    debug_sources = haskey(ENV, "SES_TRI_TRACE")
    tri_source = debug_sources ? String[] : nothing
    function mark(prefix)
        debug_sources && while length(tri_source) < length(tr.tses.tris)
            push!(tri_source, prefix)
        end
        nothing
    end
    # Set the "current face context" tag so saturation traces can attribute.
    function set_ctx(s)
        debug_sources && (ENV["SES_TRI_CTX"] = s)
        nothing
    end
    face_tris_pre = length(tr.tses.tris)
    face_tris_phase_start = length(tr.tses.tris)
    for (tfi, tf) in enumerate(ses.toric_faces)
        tf.rs_index == 0 && continue   # sentinel = deleted toric face
        set_ctx("toric_$tfi")
        if tf.type == SESFaceType.Toric
            re = rs.edges[tf.rs_index]
            if re.f1 == 0 || re.f2 == 0
                ball_triangulate_free_toric!(tr, tf)
                mark("toric_free_$tfi")
            else
                # BALL `triangulateToricFace` (triangulatedSES.C:376-394)
                # ALWAYS dispatches to non-singular (or singular for
                # `rsedge_->isSingular()`). NO 4-edge guard.
                ball_triangulate_non_singular_toric!(tr, tf)
                mark("toric_ns_$tfi")
            end
        elseif tf.type == SESFaceType.ToricSingular
            ball_triangulate_singular_toric!(tr, tf)
            mark("toric_sing_$tfi")
        end
        if haskey(ENV, "RS_FACE_TRIS")
            println(stderr, "FACE_TRIS type=toric idx=$(tfi-1) added=$(length(tr.tses.tris) - face_tris_pre)")
            face_tris_pre = length(tr.tses.tris)
        end
    end
    haskey(ENV, "RS_FACE_TRIS") &&
        println(stderr, "FACE_TRIS_SUM type=toric total_added=$(length(tr.tses.tris) - face_tris_phase_start)")
    face_tris_phase_start = length(tr.tses.tris)
    ball_partition_singular_edges!(tr)
    face_tris_pre = length(tr.tses.tris)
    for (cfi, cf) in enumerate(ses.contact_faces)
        cf.rs_index == 0 && continue   # sentinel = deleted contact face
        set_ctx("contact_$cfi")
        ball_triangulate_contact_face!(tr, cf)
        mark("contact_$cfi")
        if haskey(ENV, "RS_FACE_TRIS")
            atom_idx = (cf.rs_index != 0 && cf.rs_index <= length(ses.reduced_surface.vertices)) ?
                       ses.reduced_surface.vertices[cf.rs_index].atom : -1
            println(stderr, "FACE_TRIS type=contact idx=$(cfi-1) atom=$(atom_idx-1) nE=$(length(cf.edges)) added=$(length(tr.tses.tris) - face_tris_pre)")
            face_tris_pre = length(tr.tses.tris)
        end
    end
    haskey(ENV, "RS_FACE_TRIS") &&
        println(stderr, "FACE_TRIS_SUM type=contact total_added=$(length(tr.tses.tris) - face_tris_phase_start)")
    face_tris_phase_start = length(tr.tses.tris)
    # BALL's `triangulateSphericFaces` (triangulatedSES.C:200-247): first
    # pass tries every face; failures go on a retry queue. The retry loop
    # iterates over the failed face's SES edges; if one's first TEdge has
    # f0 set (= the toric side already triangulated it, so there's an
    # anchor to seed from) AND the edge is currently SINGULAR (so
    # `firstSESEdge` would skip it), temporarily promote it to CONCAVE
    # and re-call `triangulateSphericFace`. This lets BFT seed where it
    # would have failed.
    not_triangulated = Tuple{Int, SESFace{T}}[]
    face_tris_pre = length(tr.tses.tris)
    for (sfi, sf) in enumerate(ses.spheric_faces)
        sf.rs_index == 0 && continue   # sentinel = deleted RSFace
        set_ctx("spheric_$sfi")
        ok = ball_triangulate_spheric_face!(tr, sf)
        mark("spheric_$(sfi)$(ok ? "" : "_failed")")
        ok || push!(not_triangulated, (sfi, sf))
        if haskey(ENV, "RS_FACE_TRIS")
            println(stderr, "FACE_TRIS type=spheric idx=$(sfi-1) added=$(length(tr.tses.tris) - face_tris_pre)")
            face_tris_pre = length(tr.tses.tris)
        end
    end
    haskey(ENV, "RS_FACE_TRIS") &&
        println(stderr, "FACE_TRIS_SUM type=spheric total_added=$(length(tr.tses.tris) - face_tris_phase_start)")
    finished = isempty(not_triangulated)
    while !finished
        finished = true
        new_failed = Tuple{Int, SESFace{T}}[]
        for (sfi, sf) in not_triangulated
            ok = false
            for eid in sf.edges
                eid == 0 && continue
                te_list = tr.edge[eid]
                isempty(te_list) && continue
                first_te = tr.tses.edges[te_list[1]]
                first_te.f0 == 0 && continue
                e = tr.ses.edges[eid]
                saved_type = e.type
                e.type = SESEdgeType.Concave
                ok = ball_triangulate_spheric_face!(tr, sf)
                mark("spheric_retry_$(sfi)$(ok ? "" : "_failed")")
                e.type = saved_type
                ok && break
            end
            if ok
                finished = false
            else
                push!(new_failed, (sfi, sf))
            end
        end
        not_triangulated = new_failed
    end
    if debug_sources
        global _SES_TRI_SOURCE = tri_source
    end
    ball_extract_mesh(tr)
end

# Global slot to expose per-triangle source labels when SES_TRI_TRACE is
# set; consumed by diagnostic scripts.
_SES_TRI_SOURCE = String[]
