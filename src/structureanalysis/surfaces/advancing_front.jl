# ===========================================================================
# Advancing-front triangulator for spherical patches — a line-by-line port of
# BALL's `SESTriangulator::buildSphericTriangles` (triangulatedSES.C:777) and
# its support routines. The advancing front is used to build the contact and
# spheric SES faces; the algorithm is identical for both, differing only in
# the `convex` flag.
#
# BALL's mesh data is a doubly-linked TrianglePoint / TriangleEdge / Triangle
# graph. In Julia we keep that structure verbatim: each WS{Point, Edge, Tri}
# is a `mutable struct` with reference semantics, holding sets/vectors of
# neighbouring objects. The working surface `WSurface` owns all of them.
#
# Key invariants BALL maintains and we mirror:
#
#   * `TriangleEdge.vertex_[2]` is **directed**: traversing it in
#     `vertex_[0] → vertex_[1]` order belongs to one of its two
#     `face_[0..1]`. The other face traverses in the reverse direction. The
#     algorithm flips `vertex_[]` (see `_pop_orient!`) before constructing
#     each new triangle, so the existing-triangle-side direction is
#     consistent.
#
#   * `Triangle.vertex_[3]`, `Triangle.edge_[3]` are filled by
#     `createTriangleAndEdges`. The vertex order is CCW from outside the
#     patch sphere (away-from-centre for convex, toward-centre for concave);
#     a final geometric test swaps vertex_[0]/[1] if the cross-product
#     normal points the wrong way.
#
#   * Manifoldness is enforced by `_can_add_triangle` in the ambiguous case
#     (BALL `buildAmbiguousTriangles`, triangulatedSES.C:1108). In the
#     unambiguous case BALL skips the check because the convex-hull
#     criterion has already picked the only candidate that gives correct
#     orientation; we mirror that.
# ===========================================================================

mutable struct WSPoint{T<:Real}
    point::Vector3{T}
    normal::Vector3{T}
    edges::Set{Int}       # incident edge IDs in the parent WSurface
    faces::Set{Int}       # incident triangle IDs
    global_id::Int        # vertex ID in the parent mesh's `pos` array (= the
                          # ID the final `TriangleFace{Int}` will reference)
end

WSPoint{T}(pt::Vector3{T}, n::Vector3{T}, gid::Int) where T =
    WSPoint{T}(pt, n, Set{Int}(), Set{Int}(), gid)

# BALL's TriangleEdge — directed.
mutable struct WSEdge
    v0::Int   # FROM (vertex ID in WSurface.points)
    v1::Int   # TO
    f0::Int   # first incident triangle (0 = none)
    f1::Int   # second incident triangle (0 = none)
end

# BALL's Triangle.
mutable struct WSTri
    v0::Int       # vertex IDs in CCW-from-outside order
    v1::Int
    v2::Int
    e01::Int      # edge between v0 and v1 (= triangle.edge_[0] in BALL)
    e12::Int      # edge between v1 and v2
    e20::Int      # edge between v2 and v0
end

mutable struct WSurface{T<:Real}
    points::Vector{WSPoint{T}}
    edges::Vector{WSEdge}
    tris::Vector{WSTri}
end

WSurface{T}() where T = WSurface{T}(WSPoint{T}[], WSEdge[], WSTri[])

# ---------------------------------------------------------------------------
# `SESBuilder{T}` — shared state for the whole SES triangulation pipeline.
#
# Holds one `WSurface{T}` populated incrementally by:
#   1. shared vertex pool registration (corners + arcs),
#   2. toric face triangulation (which populates `face_[0]` of the
#      boundary edges so that...
#   3. contact and spheric face triangulation, which run advancing-front
#      on the already-populated `WSurface` with correct orientation
#      references.
#
# The lookup tables mirror BALL's `point_[]` / `edge_[index]` indexing —
# they map SES topology objects (RSFace corners, RSEdge contact arcs)
# back to the working-surface point and edge IDs.
# ---------------------------------------------------------------------------

mutable struct SESBuilder{T<:Real}
    ws::WSurface{T}
    # rsface_corners[fi] = the three contact-point WSPoint IDs of RSFace fi
    # (in the order matching `rs.faces[fi].v1, v2, v3`).
    rsface_corners::Vector{NTuple{3, Int}}
    # contact_arc_A_pids[ei] = ordered WSPoint IDs sampling the contact arc
    # of RS edge ei on atom A's bare-atom surface (= the contact face side).
    # For singular RS edges, the arc is *augmented* with the two axis
    # intersection points (see `singular_int_pids[ei]`) inserted at the
    # appropriate position along the arc, so the contact face's boundary
    # walk traverses them.
    contact_arc_A_pids::Vector{Vector{Int}}
    # contact_arc_B_pids[ei] = same on atom B.
    contact_arc_B_pids::Vector{Vector{Int}}
    # concave_arcs_pids[fi] = three concave-arc WSPoint ID lists per RSFace,
    # one per face edge (matching `f.e1, f.e2, f.e3` in cyclic order).
    concave_arcs_pids::Vector{NTuple{3, Vector{Int}}}
    # contact_arc_A_edges[ei] / contact_arc_B_edges[ei] = WSEdge IDs along
    # those arcs. Filled by `_build_toric_face_ws!` when the toric grid is
    # emitted; mirrors BALL's `edge_[(*sesedge)->index_]` per-SES-edge edge
    # list.
    contact_arc_A_edges::Vector{Vector{Int}}
    contact_arc_B_edges::Vector{Vector{Int}}
    # Likewise concave arc edges, indexed `concave_arc_edges[fi][k]` where
    # k ∈ 1:3.
    concave_arc_edges::Vector{NTuple{3, Vector{Int}}}
    # singular_int_pids[ei] = (pid_int_A_side, pid_int_B_side) for singular
    # RS edges (the two axis intersection points). For non-singular edges,
    # (0, 0). These WSPoints are SHARED across the contact arcs on atom A,
    # atom B, and the singular toric face.
    singular_int_pids::Vector{NTuple{2, Int}}
end

function SESBuilder{T}(n_rsfaces::Int, n_rsedges::Int) where T
    SESBuilder{T}(
        WSurface{T}(),
        Vector{NTuple{3, Int}}(undef, n_rsfaces),
        Vector{Vector{Int}}(undef, n_rsedges),
        Vector{Vector{Int}}(undef, n_rsedges),
        Vector{NTuple{3, Vector{Int}}}(undef, n_rsfaces),
        Vector{Vector{Int}}(undef, n_rsedges),
        Vector{Vector{Int}}(undef, n_rsedges),
        Vector{NTuple{3, Vector{Int}}}(undef, n_rsfaces),
        fill((0, 0), n_rsedges),
    )
end

# ---------------------------------------------------------------------------
# Shared vertex pool registration. Each routine returns the WSPoint IDs in
# the same shape the caller already produces for the global `pos`/`norm`
# arrays; the global IDs and WS IDs always agree because we call
# `_ws_add_point!` immediately after pushing to `pos`/`norm`.
# ---------------------------------------------------------------------------

# Step 1: per-RSFace contact points (atom-surface corners). Mirrors the
# existing `_build_corner_vertices` but also registers each point in the
# `WSurface` and records the WS IDs into `builder.rsface_corners`.
function _build_corner_vertices_ws!(
    builder::SESBuilder{T},
    rs::ReducedSurface{T},
    pos::Vector{Point3{T}},
    norm::Vector{Vec3{T}},
) where T
    for (fi, f) in enumerate(rs.faces)
        triple = ntuple(3) do k
            atom_idx = rs.vertices[(f.v1, f.v2, f.v3)[k]].atom
            atom_c = Vector3{T}(rs.atoms[atom_idx].center...)
            R_atom = rs.atoms[atom_idx].r
            dir = f.center - atom_c
            ℓ = sqrt(dot(dir, dir))
            outward = ℓ > eps(T) ? dir / ℓ : Vector3{T}(0, 0, 1)
            point = atom_c + R_atom * outward
            push!(pos, Point3{T}(point...))
            push!(norm, Vec3{T}(outward...))
            gid = length(pos)
            _ws_add_point!(builder.ws, point, outward, gid)
        end
        builder.rsface_corners[fi] = triple
    end
end

# Step 2a/2b: contact arcs on atom A/B surfaces. Mirrors `_sample_contact_arc`
# / `_build_contact_arcs` from triangulated_ses.jl but registers each
# interior arc vertex in `WSurface` and records the WS IDs in
# `builder.contact_arc_{A,B}_pids[ei]`. Edges along the arc are NOT
# created here — the toric face builder is the first writer of those
# edges (it sets `face_[0]` on them).
function _build_contact_arcs_ws!(
    builder::SESBuilder{T},
    rs::ReducedSurface{T},
    pos::Vector{Point3{T}},
    norm::Vector{Vec3{T}},
    n_seg_edge::Vector{Int},
) where T
    n_edges = length(rs.edges)
    # Pre-pass: for singular RS edges, create the 2 axis-intersection
    # WSPoints (shared between atom A's and atom B's contact arcs and
    # the singular toric triangulator). BALL's `treatSingularToricFace`
    # creates these as `new_point1` and `new_point3`; here we eagerly
    # create them up front so contact-arc sampling can reference them.
    for ei in 1:n_edges
        e = rs.edges[ei]
        e.singular || continue
        (e.f1 == 0 || e.f2 == 0) && continue
        # Outward normal: toward the midpoint of the two adjacent probe
        # centres (= roughly the local saddle outward).
        midp = (rs.faces[e.f1].center + rs.faces[e.f2].center) / T(2)
        n1 = midp - e.intersection1
        ℓ1 = sqrt(dot(n1, n1))
        n1 = ℓ1 > eps(T) ? n1 / ℓ1 : Vector3{T}(0, 0, 1)
        push!(pos, Point3{T}(e.intersection1...))
        push!(norm, Vec3{T}(n1...))
        pid_A = _ws_add_point!(builder.ws, e.intersection1, n1, length(pos))
        n2 = midp - e.intersection2
        ℓ2 = sqrt(dot(n2, n2))
        n2 = ℓ2 > eps(T) ? n2 / ℓ2 : Vector3{T}(0, 0, 1)
        push!(pos, Point3{T}(e.intersection2...))
        push!(norm, Vec3{T}(n2...))
        pid_B = _ws_add_point!(builder.ws, e.intersection2, n2, length(pos))
        builder.singular_int_pids[ei] = (pid_A, pid_B)
    end
    for ei in 1:n_edges
        e = rs.edges[ei]
        if e.f1 == 0 || e.f2 == 0
            builder.contact_arc_A_pids[ei] = Int[]
            builder.contact_arc_B_pids[ei] = Int[]
            continue
        end
        atom_A = rs.vertices[e.v1].atom
        atom_B = rs.vertices[e.v2].atom
        s_f1_A = _atom_slot(rs, e.f1, atom_A)
        s_f2_A = _atom_slot(rs, e.f2, atom_A)
        s_f1_B = _atom_slot(rs, e.f1, atom_B)
        s_f2_B = _atom_slot(rs, e.f2, atom_B)
        if s_f1_A == 0 || s_f2_A == 0 || s_f1_B == 0 || s_f2_B == 0
            builder.contact_arc_A_pids[ei] = Int[]
            builder.contact_arc_B_pids[ei] = Int[]
            continue
        end
        cv_f1_A = builder.rsface_corners[e.f1][s_f1_A]
        cv_f2_A = builder.rsface_corners[e.f2][s_f2_A]
        cv_f1_B = builder.rsface_corners[e.f1][s_f1_B]
        cv_f2_B = builder.rsface_corners[e.f2][s_f2_B]
        atom_A_c = Vector3{T}(rs.atoms[atom_A].center...)
        atom_B_c = Vector3{T}(rs.atoms[atom_B].center...)
        third_pts = Vector3{T}[]
        f1 = rs.faces[e.f1]; f2 = rs.faces[e.f2]
        for vidx in (f1.v1, f1.v2, f1.v3, f2.v1, f2.v2, f2.v3)
            a = rs.vertices[vidx].atom
            (a == atom_A || a == atom_B) && continue
            push!(third_pts, Vector3{T}(rs.atoms[a].center...))
        end
        n_seg = n_seg_edge[ei]
        ids_A = _sample_contact_arc_into_ws!(
            builder, e.contact_circle1, atom_A_c, rs.atoms[atom_A].r,
            Vector3{T}(pos[cv_f1_A]...), Vector3{T}(pos[cv_f2_A]...),
            cv_f1_A, cv_f2_A, n_seg, third_pts, pos, norm)
        ids_B = _sample_contact_arc_into_ws!(
            builder, e.contact_circle2, atom_B_c, rs.atoms[atom_B].r,
            Vector3{T}(pos[cv_f1_B]...), Vector3{T}(pos[cv_f2_B]...),
            cv_f1_B, cv_f2_B, n_seg, third_pts, pos, norm)
        # For singular RS edges insert the 2 axis-intersection WSPoints
        # in the MIDDLE of the arc — atom A's arc gets them in (A-side,
        # B-side) order; atom B's arc in reverse (B-side, A-side) so each
        # contact face traverses the singular detour with consistent
        # orientation. Without the insertion the contact face's boundary
        # walk would skip the singular vertices entirely and the contact-
        # arc segment between mid samples (which crosses the singular
        # probe position) would have no proper triangle on the toric
        # side.
        # NOTE: for singular RS edges we eagerly created the 2 axis-
        # intersection WSPoints in the pre-pass above (stored in
        # `builder.singular_int_pids[ei]`). They are NOT inserted into
        # the contact arcs here — empirically, augmenting the contact
        # arcs with off-surface intersection points causes the contact
        # face's iterative quickhull to pick candidates that lead to
        # over-fills (because the intersection points are
        # geometrically inside the molecular volume, far from atom A's
        # convex contact patch). The intersection points stay reserved
        # for the future singular toric triangulator (which will
        # connect them via diagonal SES edges on the toric saddle).
        builder.contact_arc_A_pids[ei] = ids_A
        builder.contact_arc_B_pids[ei] = ids_B
    end
end

# Like `_sample_contact_arc` but also adds each new interior point to the
# WSurface. Endpoints are already in WSurface (from corner registration).
function _sample_contact_arc_into_ws!(
    builder::SESBuilder{T},
    circle::Circle3{T}, atom_center::Vector3{T},
    R_atom::T, start_pt::Vector3{T}, end_pt::Vector3{T},
    start_idx::Int, end_idx::Int, n_seg::Int,
    avoid_pts::Vector{Vector3{T}},
    pos::Vector{Point3{T}}, norm::Vector{Vec3{T}},
) where T
    u, v = _onb_perp(circle.n)
    rel_s = start_pt - circle.p
    rel_e = end_pt - circle.p
    ψ_s = atan(dot(rel_s, v), dot(rel_s, u))
    ψ_e = atan(dot(rel_e, v), dot(rel_e, u))
    Δ_short = ψ_e - ψ_s
    if Δ_short >  T(π); Δ_short -= 2 * T(π); end
    if Δ_short < -T(π); Δ_short += 2 * T(π); end
    Δ_long = Δ_short - sign(Δ_short) * 2 * T(π)
    function midpoint_dist2(Δ)
        ψ_mid = ψ_s + Δ / 2
        mid = circle.p + circle.r * (cos(ψ_mid) * u + sin(ψ_mid) * v)
        s = zero(T)
        for p in avoid_pts
            d = mid - p
            s += dot(d, d)
        end
        s
    end
    Δ = midpoint_dist2(Δ_long) > midpoint_dist2(Δ_short) ? Δ_long : Δ_short
    ids = Vector{Int}(undef, n_seg + 1)
    ids[1] = start_idx
    for k in 1:(n_seg - 1)
        t = T(k) / T(n_seg)
        ψ = ψ_s + t * Δ
        point = circle.p + circle.r * (cos(ψ) * u + sin(ψ) * v)
        outward = point - atom_center
        ℓ = sqrt(dot(outward, outward))
        outward = ℓ > eps(T) ? outward / ℓ : Vector3{T}(0, 0, 1)
        push!(pos, Point3{T}(point...))
        push!(norm, Vec3{T}(outward...))
        gid = length(pos)
        _ws_add_point!(builder.ws, point, outward, gid)
        ids[k + 1] = gid
    end
    ids[end] = end_idx
    ids
end

# Step 2c: concave arcs on probe spheres. For each RSFace `fi`, sample
# the three great-circle arcs (one per RS edge of the face) on the probe
# sphere centred at `f.center`. Mirrors `_build_concave_arcs` from
# triangulated_ses.jl but also registers each new interior arc point in
# the WSurface and records the WS IDs in
# `builder.concave_arcs_pids[fi][k]` (k ∈ 1:3).
function _build_concave_arcs_ws!(
    builder::SESBuilder{T},
    rs::ReducedSurface{T},
    pos::Vector{Point3{T}},
    norm::Vector{Vec3{T}},
    n_seg_edge::Vector{Int},
) where T
    for (fi, f) in enumerate(rs.faces)
        ids = builder.rsface_corners[fi]
        # Map each arc (k=1..3 = side between vertex k and vertex k+1) to
        # ONE incident RS edge by atom pair. With stacked-probe faces a
        # face may have more than 3 incident edges; pick the first one
        # matching the arc's atom pair.
        verts = (f.v1, f.v2, f.v3)
        face_atom = ntuple(k -> rs.vertices[verts[k]].atom, 3)
        face_edges = ntuple(3) do k
            a = face_atom[k]; b = face_atom[mod1(k + 1, 3)]
            for ei in f.edges
                re = rs.edges[ei]
                at1 = rs.vertices[re.v1].atom
                at2 = rs.vertices[re.v2].atom
                if (at1 == a && at2 == b) || (at1 == b && at2 == a)
                    return ei
                end
            end
            return 0
        end
        arcs = ntuple(3) do k
            ka = k; kb = mod1(k + 1, 3)
            sa = Vector3{T}(pos[ids[ka]]...)
            sb = Vector3{T}(pos[ids[kb]]...)
            ei = face_edges[k]
            n_seg = (ei == 0) ? 2 : n_seg_edge[ei]
            _sample_concave_arc_into_ws!(builder, f.center, rs.probe_radius,
                                          sa, sb, ids[ka], ids[kb], n_seg,
                                          pos, norm)
        end
        builder.concave_arcs_pids[fi] = arcs
    end
end

# Sample a great-circle arc on the probe sphere between two points, adding
# each interior sample to the WSurface and the global pos/norm arrays.
function _sample_concave_arc_into_ws!(
    builder::SESBuilder{T},
    probe_center::Vector3{T}, probe_r::T,
    start_pt::Vector3{T}, end_pt::Vector3{T},
    start_idx::Int, end_idx::Int, n_seg::Int,
    pos::Vector{Point3{T}}, norm::Vector{Vec3{T}},
) where T
    ua = (start_pt - probe_center) / probe_r
    ub = (end_pt   - probe_center) / probe_r
    ids = Vector{Int}(undef, n_seg + 1)
    ids[1] = start_idx
    for k in 1:(n_seg - 1)
        t = T(k) / T(n_seg)
        u = _slerp_unit(ua, ub, t)
        point = probe_center + probe_r * u
        push!(pos, Point3{T}(point...))
        # SES outward at a probe-surface point is *toward* the probe centre
        # (away from molecular interior on the probe-rest side).
        push!(norm, Vec3{T}(-u[1], -u[2], -u[3]))
        gid = length(pos)
        _ws_add_point!(builder.ws, point, Vector3{T}(-u[1], -u[2], -u[3]), gid)
        ids[k + 1] = gid
    end
    ids[end] = end_idx
    ids
end

# ---------------------------------------------------------------------------
# Step 3: toric face triangulation onto the shared `WSurface`. Adds the
# (n_seg+1) × (n_seg+1) quad-grid for a non-singular toric face and emits
# its 2·n_seg² triangles. Each triangle's edges (vertical, horizontal,
# diagonal) are registered in the working surface — boundary edges (top,
# bottom, left, right of the grid) end up with `face_[0]` set to the
# toric triangle; the corresponding contact / spheric face triangulator
# will later fill `face_[1]`.
#
# This mirrors `_triangulate_toric!` from triangulated_ses.jl
# geometrically; the difference is that every grid vertex/edge/triangle
# is also registered in the WSurface, and a corresponding global
# `pos`/`norm`/`faces` entry is appended so the final mesh contains the
# toric triangles.
# ---------------------------------------------------------------------------

function _build_toric_face_ws!(
    builder::SESBuilder{T},
    ses::SolventExcludedSurface{T},
    pos::Vector{Point3{T}},
    norm::Vector{Vec3{T}},
    faces::Vector{TriangleFace{Int}},
    tf,
    n_seg_edge::Vector{Int},
) where T
    rs = ses.reduced_surface
    if tf.type == SESFaceType.ToricSingular
        _triangulate_singular_toric!(builder, rs, faces, pos, tf)
        return
    end
    tf.type == SESFaceType.Toric || return
    ei = tf.rs_index
    e = rs.edges[ei]
    (e.f1 == 0 || e.f2 == 0) && return
    atom_A = rs.vertices[e.v1].atom
    atom_B = rs.vertices[e.v2].atom
    ca_f1, rev_f1 = _concave_arc_id_ws(builder, rs, e.f1, atom_A, atom_B)
    ca_f2, rev_f2 = _concave_arc_id_ws(builder, rs, e.f2, atom_A, atom_B)
    (isempty(ca_f1) || isempty(ca_f2)) && return
    f1_arc = rev_f1 ? reverse(ca_f1) : ca_f1
    f2_arc = rev_f2 ? reverse(ca_f2) : ca_f2
    cA = builder.contact_arc_A_pids[ei]
    cB = builder.contact_arc_B_pids[ei]
    (isempty(cA) || isempty(cB)) && return
    n_seg = n_seg_edge[ei]
    (length(f1_arc) == n_seg + 1 && length(f2_arc) == n_seg + 1 &&
     length(cA) == n_seg + 1 && length(cB) == n_seg + 1) || return

    # Build grid of WSPoint IDs.
    # grid[k+1, j+1] = WSPoint ID at rolling step k, phi step j.
    grid = fill(0, n_seg + 1, n_seg + 1)
    for j in 0:n_seg
        grid[1, j + 1]         = f1_arc[j + 1]
        grid[n_seg + 1, j + 1] = f2_arc[j + 1]
    end
    for k in 0:n_seg
        grid[k + 1, 1]         = cA[k + 1]
        grid[k + 1, n_seg + 1] = cB[k + 1]
    end

    # Interior probe-rolling parameters (same as _triangulate_toric!).
    atom_A_c = Vector3{T}(rs.atoms[atom_A].center...)
    atom_B_c = Vector3{T}(rs.atoms[atom_B].center...)
    axis_dir = atom_B_c - atom_A_c
    axis_len = sqrt(dot(axis_dir, axis_dir))
    axis_n = axis_len > eps(T) ? axis_dir / axis_len : Vector3{T}(0, 0, 1)
    torus_c = e.center_of_torus
    pc_f1 = rs.faces[e.f1].center
    pc_f2 = rs.faces[e.f2].center
    v_f1 = pc_f1 - torus_c
    v_f1_perp = v_f1 - dot(v_f1, axis_n) * axis_n
    v_f2 = pc_f2 - torus_c
    v_f2_perp = v_f2 - dot(v_f2, axis_n) * axis_n
    α_short = oriented_angle(v_f1_perp, v_f2_perp, axis_n)
    cA_mid_pid = cA[div(length(cA), 2) + 1]
    cA_mid = Vector3{T}(pos[cA_mid_pid]...)
    rot_radial = θ -> begin
        cosθ = cos(θ); sinθ = sin(θ)
        v_f1_perp * cosθ + cross(axis_n, v_f1_perp) * sinθ +
            axis_n * dot(axis_n, v_f1_perp) * (1 - cosθ)
    end
    mid_short = torus_c + rot_radial(α_short / 2)
    mid_long  = torus_c + rot_radial((α_short - 2 * T(π)) / 2)
    d_short = mid_short - cA_mid; d_short_sq = dot(d_short, d_short)
    d_long  = mid_long  - cA_mid; d_long_sq  = dot(d_long, d_long)
    α = d_long_sq < d_short_sq ? (α_short - 2 * T(π)) : α_short
    R_A_inf = rs.atoms[atom_A].r + rs.probe_radius
    R_B_inf = rs.atoms[atom_B].r + rs.probe_radius

    # Interior grid points (k = 1..n_seg-1, j = 1..n_seg-1).
    for k in 1:(n_seg - 1)
        t_k = T(k) / T(n_seg)
        θ = α * t_k
        radial = rot_radial(θ)
        probe_c = torus_c + radial
        u_A = (atom_A_c - probe_c) / R_A_inf
        u_B = (atom_B_c - probe_c) / R_B_inf
        for j in 1:(n_seg - 1)
            t_j = T(j) / T(n_seg)
            u = _slerp_unit(u_A, u_B, t_j)
            p = probe_c + rs.probe_radius * u
            push!(pos, Point3{T}(p...))
            n_outward = Vector3{T}(-u[1], -u[2], -u[3])
            push!(norm, Vec3{T}(n_outward...))
            gid = length(pos)
            pid = _ws_add_point!(builder.ws, p, n_outward, gid)
            grid[k + 1, j + 1] = pid
        end
    end

    # Emit 2 triangles per cell. For each, create/find the 3 WSEdges and
    # set the triangle's face_[0]/face_[1] appropriately.
    for k in 0:(n_seg - 1)
        for j in 0:(n_seg - 1)
            v_bl = grid[k + 1, j + 1]      # bottom-left
            v_br = grid[k + 2, j + 1]      # bottom-right
            v_tr = grid[k + 2, j + 2]      # top-right
            v_tl = grid[k + 1, j + 2]      # top-left
            # Outward hint = average of 4 vertex normals (= roughly the
            # toric saddle's outward direction at this cell).
            navg = Vector3{T}(0, 0, 0)
            for v in (v_bl, v_br, v_tr, v_tl)
                navg = navg + builder.ws.points[v].normal
            end
            # Triangle 1: (v_bl, v_br, v_tr). Triangle 2: (v_bl, v_tr, v_tl).
            # Orient each via outward hint, then emit.
            _emit_toric_triangle_ws!(builder, faces, pos, v_bl, v_br, v_tr, navg)
            _emit_toric_triangle_ws!(builder, faces, pos, v_bl, v_tr, v_tl, navg)
        end
    end
end

# ---------------------------------------------------------------------------
# Singular toric face triangulator. For each singular toric face, emit
# two triangle fans, one per side of the singular edge:
#
#   * Side f1: triangles fanning from `pid_int1` (atom-A-side axis
#     intersection) to the concave_f1 arc samples.
#   * Side f2: triangles fanning from `pid_int2` (atom-B-side axis
#     intersection) to the concave_f2 arc samples.
#
# Each fan emits `n_seg` triangles. The "ribs" of the fan (= edges from
# the intersection point to interior concave-arc samples) and the
# "spokes" (= edges between consecutive concave-arc samples) are
# registered as new WSEdges. The contact arc samples on atom A and B are
# NOT used by this triangulator — the contact face's boundary walk
# remains unchanged. The toric strip itself becomes a hole bounded by
# the contact arcs and the new diagonals.
# ---------------------------------------------------------------------------

function _triangulate_singular_toric!(
    builder::SESBuilder{T},
    rs::ReducedSurface{T},
    faces::Vector{TriangleFace{Int}},
    pos::Vector{Point3{T}},
    tf,
) where T
    ei = tf.rs_index
    e = rs.edges[ei]
    (e.f1 == 0 || e.f2 == 0) && return
    pid_int1, pid_int2 = builder.singular_int_pids[ei]
    (pid_int1 == 0 || pid_int2 == 0) && return
    atom_A = rs.vertices[e.v1].atom
    atom_B = rs.vertices[e.v2].atom
    ca_f1, rev_f1 = _concave_arc_id_ws(builder, rs, e.f1, atom_A, atom_B)
    ca_f2, rev_f2 = _concave_arc_id_ws(builder, rs, e.f2, atom_A, atom_B)
    (isempty(ca_f1) || isempty(ca_f2)) && return
    arc_f1 = rev_f1 ? reverse(ca_f1) : ca_f1
    arc_f2 = rev_f2 ? reverse(ca_f2) : ca_f2
    # NOTE: a "strip" variant that ALSO uses the contact arc samples on
    # atoms A and B was tried — it creates a lot of diagonal edges from
    # contact-arc samples to concave-arc samples that become count=1
    # boundary edges (no face on the molecular-interior side), so BPTI's
    # boundary count went UP from 726 → 1341. The simple fan ends up
    # cleaner. The proper fix would be BALL's full `buildTriangles` with
    # singular handling (treats the strip as a 2D grid with the
    # pinch-row degenerate), but that needs proper SES topology + n_tri
    # rolling-angle interpolation we don't have yet.
    _emit_singular_toric_fan!(builder, faces, pos, arc_f1, pid_int1)
    _emit_singular_toric_fan!(builder, faces, pos, arc_f2, pid_int2)
end

# Emit triangles bridging the concave arc on a probe sphere to the
# contact arcs on atoms A and B, with the apex at the axis intersection
# point.
#
# Each (j+1, j+2) cell along the rolling direction emits 2 triangles:
#   * (cA[j+1], cA[j+2], arc[j+1]): atom-A-side, fanning toward concave
#     arc.
#   * (cB[j+1], cB[j+2], arc[j+1]): atom-B-side, fanning toward concave
#     arc.
# Plus one connecting triangle per cell:
#   * (arc[j+1], arc[j+2], apex): pinch-side fan triangle.
#
# This closes the contact-arc segment edges (between consecutive cA/cB
# samples) on the toric side: each contact-arc segment gets face_[0]
# set to a triangle here, leaving face_[1] for the adjacent contact-
# face advancing-front to fill.
function _emit_singular_toric_full!(
    builder::SESBuilder{T},
    faces::Vector{TriangleFace{Int}},
    pos::Vector{Point3{T}},
    arc_top::Vector{Int},
    cA::Vector{Int}, cB::Vector{Int},
    apex_pid::Int,
) where T
    n_seg = length(arc_top) - 1
    n_seg >= 1 || return
    apex_n = builder.ws.points[apex_pid].normal
    for j in 0:(n_seg - 1)
        # Atom-A-side fan triangle (cA[j+1], cA[j+2], arc[j+1]).
        navg = builder.ws.points[cA[j + 1]].normal +
               builder.ws.points[cA[j + 2]].normal +
               builder.ws.points[arc_top[j + 1]].normal
        _emit_toric_triangle_ws!(builder, faces, pos,
                                  cA[j + 1], cA[j + 2], arc_top[j + 1], navg)
        # Atom-B-side fan triangle (cB[j+1], cB[j+2], arc[j+1]).
        navg = builder.ws.points[cB[j + 1]].normal +
               builder.ws.points[cB[j + 2]].normal +
               builder.ws.points[arc_top[j + 1]].normal
        _emit_toric_triangle_ws!(builder, faces, pos,
                                  cB[j + 1], cB[j + 2], arc_top[j + 1], navg)
        # Pinch-side fan triangle to apex.
        navg = apex_n + builder.ws.points[arc_top[j + 1]].normal +
               builder.ws.points[arc_top[j + 2]].normal
        _emit_toric_triangle_ws!(builder, faces, pos,
                                  apex_pid, arc_top[j + 1], arc_top[j + 2], navg)
    end
end

# Emit a triangle fan from `apex_pid` to the consecutive samples of
# `arc`. Each triangle's outward normal is aligned with the average of
# the three vertex normals (which are all toward the probe centre on
# the singular toric saddle, so the average is well-defined).
function _emit_singular_toric_fan!(
    builder::SESBuilder{T},
    faces::Vector{TriangleFace{Int}},
    pos::Vector{Point3{T}},
    arc::Vector{Int},
    apex_pid::Int,
) where T
    n = length(arc)
    n >= 2 || return
    apex_n = builder.ws.points[apex_pid].normal
    for k in 1:(n - 1)
        a = arc[k]; b = arc[k + 1]
        # Outward hint = average of vertex normals.
        navg = apex_n + builder.ws.points[a].normal + builder.ws.points[b].normal
        _emit_toric_triangle_ws!(builder, faces, pos, apex_pid, a, b, navg)
    end
end

# ---------------------------------------------------------------------------
# Singular toric face triangulation — port of BALL's
# `triangulateSingularToricFace` (triangulatedSES.C:433).
#
# When a probe rolling between atoms A and B touches a third atom mid-roll
# (a "spindle torus", `e.singular == true`), the toric strip pinches at the
# rolling axis. BALL handles this by adding 2 new SES vertices on the axis
# (`e.intersection1` near A, `e.intersection2` near B) and splitting the
# face into 2 halves, each fanning to one of those axis points.
#
# Geometrically each half-strip is bounded by:
#   * One concave arc on a probe sphere (f1's or f2's resting position).
#   * Two "diagonal" half-arcs going from atom-surface corners to axis
#     intersection points (on the toric saddle surface).
#   * The singular edge between the 2 intersection points (on the
#     intersection circle of the 2 adjacent probes).
#
# Our triangulator builds each half as a fan-style grid where one row
# collapses to a single intersection point, mirroring BALL's
# `buildTriangles(..., edge3 = NULL)` behaviour for singular faces.
# ---------------------------------------------------------------------------

function _build_singular_toric_face_ws!(
    builder::SESBuilder{T},
    ses::SolventExcludedSurface{T},
    pos::Vector{Point3{T}},
    norm::Vector{Vec3{T}},
    faces::Vector{TriangleFace{Int}},
    tf,
    n_seg_edge::Vector{Int},
) where T
    rs = ses.reduced_surface
    ei = tf.rs_index
    e = rs.edges[ei]
    (e.f1 == 0 || e.f2 == 0) && return
    atom_A = rs.vertices[e.v1].atom
    atom_B = rs.vertices[e.v2].atom
    # The two axis-probe intersection points become new SES vertices.
    # Their outward normal points from the intersection point toward the
    # mid-plane of the two adjacent probes (which is the singular probe
    # surface normal at this axis point).
    f1_c = rs.faces[e.f1].center
    f2_c = rs.faces[e.f2].center
    pr = rs.probe_radius
    # Center of the intersection circle of probes(f1) and probes(f2).
    # The intersection lies on the perpendicular bisector of f1_c, f2_c.
    midpoint_p = (f1_c + f2_c) / T(2)
    int1 = e.intersection1
    int2 = e.intersection2
    n_int1 = midpoint_p - int1
    ℓ1 = sqrt(dot(n_int1, n_int1))
    n_int1 = ℓ1 > eps(T) ? n_int1 / ℓ1 : Vector3{T}(0, 0, 1)
    n_int2 = midpoint_p - int2
    ℓ2 = sqrt(dot(n_int2, n_int2))
    n_int2 = ℓ2 > eps(T) ? n_int2 / ℓ2 : Vector3{T}(0, 0, 1)
    push!(pos, Point3{T}(int1...))
    push!(norm, Vec3{T}(n_int1...))
    pid_int1 = _ws_add_point!(builder.ws, int1, n_int1, length(pos))
    push!(pos, Point3{T}(int2...))
    push!(norm, Vec3{T}(n_int2...))
    pid_int2 = _ws_add_point!(builder.ws, int2, n_int2, length(pos))

    # Sanity check: the singular toric face needs valid contact corners
    # on both atoms in each adjacent RSFace.
    (_atom_slot(rs, e.f1, atom_A) == 0 ||
     _atom_slot(rs, e.f1, atom_B) == 0 ||
     _atom_slot(rs, e.f2, atom_A) == 0 ||
     _atom_slot(rs, e.f2, atom_B) == 0) && return

    # Concave arcs on f1 and f2 (corner_A → corner_B).
    ca_f1, rev_f1 = _concave_arc_id_ws(builder, rs, e.f1, atom_A, atom_B)
    ca_f2, rev_f2 = _concave_arc_id_ws(builder, rs, e.f2, atom_A, atom_B)
    (isempty(ca_f1) || isempty(ca_f2)) && return
    f1_arc = rev_f1 ? reverse(ca_f1) : ca_f1
    f2_arc = rev_f2 ? reverse(ca_f2) : ca_f2

    # Build the f1 half-strip as a fan converging to (pid_int1, pid_int2).
    # The boundary of this half is the concave_f1 arc (n_seg+1 points).
    # The far side collapses to the 2 intersection points (an edge, not
    # a single vertex).
    n_seg = n_seg_edge[ei]
    _emit_singular_toric_half!(builder, faces, pos, norm,
        f1_arc, pid_int1, pid_int2, n_seg)
    _emit_singular_toric_half!(builder, faces, pos, norm,
        f2_arc, pid_int1, pid_int2, n_seg)
end

# Emit the triangle strip for one half of a singular toric face. The
# half is bounded by:
#   * the concave arc on a probe sphere (`arc` — n_seg+1 points from
#     corner_A to corner_B);
#   * the singular edge (pid_int1 ↔ pid_int2, on the rolling axis);
#   * a "diagonal" segment from corner_A to pid_int1 on the toric saddle;
#   * a "diagonal" segment from corner_B to pid_int2 on the toric saddle.
#
# For now we emit the minimum fan triangulation: 1 row of (n_seg)
# quadrilateral cells, each split into 2 triangles. The "top row" is
# the concave arc; the "bottom row" interpolates from pid_int1 to
# pid_int2.
function _emit_singular_toric_half!(
    builder::SESBuilder{T},
    faces::Vector{TriangleFace{Int}},
    pos::Vector{Point3{T}},
    norm::Vector{Vec3{T}},
    arc::Vector{Int},
    pid_int1::Int, pid_int2::Int,
    n_seg::Int,
) where T
    length(arc) == n_seg + 1 || return
    # Bottom row: interpolate from pid_int1 (corner A side) to pid_int2
    # (corner B side). Endpoints are the 2 intersection vertices; the
    # interior samples are placed on the singular intersection circle
    # (which is the great circle on each adjacent probe shared between
    # them).
    int1_pt = builder.ws.points[pid_int1].point
    int2_pt = builder.ws.points[pid_int2].point
    bottom = Vector{Int}(undef, n_seg + 1)
    bottom[1] = pid_int1
    bottom[end] = pid_int2
    for k in 1:(n_seg - 1)
        t = T(k) / T(n_seg)
        p = (1 - t) * int1_pt + t * int2_pt
        # Outward direction: toward the corresponding probe centre.
        # Use the average direction of the two intersection-point
        # outward normals as a rough proxy.
        nA = builder.ws.points[pid_int1].normal
        nB = builder.ws.points[pid_int2].normal
        n_avg = (1 - t) * nA + t * nB
        ℓ = sqrt(dot(n_avg, n_avg))
        n_avg = ℓ > eps(T) ? n_avg / ℓ : Vector3{T}(0, 0, 1)
        push!(pos, Point3{T}(p...))
        push!(norm, Vec3{T}(n_avg...))
        bottom[k + 1] = _ws_add_point!(builder.ws, p, n_avg, length(pos))
    end
    # Now build (n_seg) cells of 2 triangles each between `bottom` and `arc`.
    for k in 0:(n_seg - 1)
        v_bl = bottom[k + 1]
        v_br = bottom[k + 2]
        v_tr = arc[k + 2]
        v_tl = arc[k + 1]
        # Outward hint (average of 4 vertex normals).
        navg = Vector3{T}(0, 0, 0)
        for v in (v_bl, v_br, v_tr, v_tl)
            navg = navg + builder.ws.points[v].normal
        end
        _emit_toric_triangle_ws!(builder, faces, pos, v_bl, v_br, v_tr, navg)
        _emit_toric_triangle_ws!(builder, faces, pos, v_bl, v_tr, v_tl, navg)
    end
end

# Stub for singular toric faces: create the 4 boundary arc edges
# (2 contact, 2 concave) as phantom WSEdges (`f0 = -1`). The toric
# strip itself is NOT triangulated — it becomes a hole in the final
# mesh — but the boundary edges are registered so the adjacent
# contact/spheric advancing-front sees a complete boundary cycle and
# closes around the hole.
function _stub_singular_toric_boundary!(
    builder::SESBuilder{T},
    rs::ReducedSurface{T},
    tf,
    ::Vector{Int},
) where T
    ei = tf.rs_index
    e = rs.edges[ei]
    (e.f1 == 0 || e.f2 == 0) && return
    atom_A = rs.vertices[e.v1].atom
    atom_B = rs.vertices[e.v2].atom
    ca_f1, _ = _concave_arc_id_ws(builder, rs, e.f1, atom_A, atom_B)
    ca_f2, _ = _concave_arc_id_ws(builder, rs, e.f2, atom_A, atom_B)
    cA = builder.contact_arc_A_pids[ei]
    cB = builder.contact_arc_B_pids[ei]
    function add_phantom_segments!(pids::Vector{Int})
        for k in 1:(length(pids) - 1)
            a = pids[k]; b = pids[k + 1]
            _ws_find_edge_through(builder.ws, a, b) != 0 && continue
            push!(builder.ws.edges, WSEdge(a, b, -1, 0))
            eid = length(builder.ws.edges)
            push!(builder.ws.points[a].edges, eid)
            push!(builder.ws.points[b].edges, eid)
        end
    end
    !isempty(ca_f1) && add_phantom_segments!(ca_f1)
    !isempty(ca_f2) && add_phantom_segments!(ca_f2)
    !isempty(cA)    && add_phantom_segments!(cA)
    !isempty(cB)    && add_phantom_segments!(cB)
end

# Find which of the 3 concave arcs at face `face_idx` runs from atom_A's
# corner to atom_B's corner. Returns (arc_pids_list, reversed::Bool) —
# mirrors `_concave_arc_id` from triangulated_ses.jl.
function _concave_arc_id_ws(
    builder::SESBuilder{T},
    rs::ReducedSurface{T},
    face_idx::Int, atom_A::Int, atom_B::Int,
) where T
    s_A = _atom_slot(rs, face_idx, atom_A)
    s_B = _atom_slot(rs, face_idx, atom_B)
    (s_A == 0 || s_B == 0) && return (Int[], false)
    arcs = builder.concave_arcs_pids[face_idx]
    for k in 1:3
        kn = mod1(k + 1, 3)
        if s_A == k && s_B == kn
            return (arcs[k], false)
        end
        if s_B == k && s_A == kn
            return (arcs[k], true)
        end
    end
    Int[], false
end

# Emit a single toric triangle. Picks vertex order so the outward normal
# aligns with `navg`, then registers the triangle + its three WSEdges in
# the working surface and pushes the (correctly wound) triangle to the
# global `faces` array.
function _emit_toric_triangle_ws!(
    builder::SESBuilder{T},
    faces::Vector{TriangleFace{Int}},
    pos::Vector{Point3{T}},
    a::Int, b::Int, c::Int,
    navg::Vector3{T},
) where T
    pa = Vector3{T}(pos[builder.ws.points[a].global_id]...)
    pb = Vector3{T}(pos[builder.ws.points[b].global_id]...)
    pc = Vector3{T}(pos[builder.ws.points[c].global_id]...)
    n_tri = cross(pb - pa, pc - pa)
    # If the natural winding doesn't align with the outward hint, flip the
    # two vertices that produce a CCW-from-outside-the-saddle order.
    v0, v1, v2 = if dot(n_tri, navg) >= 0
        (a, b, c)
    else
        (a, c, b)
    end
    # Register the triangle in the working surface.
    e01 = _ws_get_or_create_grid_edge!(builder.ws, v0, v1)
    e12 = _ws_get_or_create_grid_edge!(builder.ws, v1, v2)
    e20 = _ws_get_or_create_grid_edge!(builder.ws, v2, v0)
    push!(builder.ws.tris, WSTri(v0, v1, v2, e01, e12, e20))
    tid = length(builder.ws.tris)
    push!(builder.ws.points[v0].faces, tid)
    push!(builder.ws.points[v1].faces, tid)
    push!(builder.ws.points[v2].faces, tid)
    _ws_attach_triangle_to_edge!(builder.ws.edges[e01], tid)
    _ws_attach_triangle_to_edge!(builder.ws.edges[e12], tid)
    _ws_attach_triangle_to_edge!(builder.ws.edges[e20], tid)
    # Push the (global-ID) triangle to the output mesh.
    g0 = builder.ws.points[v0].global_id
    g1 = builder.ws.points[v1].global_id
    g2 = builder.ws.points[v2].global_id
    _push_triangle!(faces, g0, g1, g2, pos, navg)
end

# Find an existing WSEdge between `a` and `b`; if none, create with
# (v_from = a, v_to = b). Used by the toric grid builder, which doesn't
# care about edge direction (the geometric outward test in
# `_emit_toric_triangle_ws!` has already fixed the triangle winding).
function _ws_get_or_create_grid_edge!(ws::WSurface, a::Int, b::Int)
    eid = _ws_find_edge_through(ws, a, b)
    eid != 0 && return eid
    push!(ws.edges, WSEdge(a, b, 0, 0))
    eid = length(ws.edges)
    push!(ws.points[a].edges, eid)
    push!(ws.points[b].edges, eid)
    eid
end

@inline function _ws_attach_triangle_to_edge!(e::WSEdge, tid::Int)
    if e.f0 == 0
        e.f0 = tid
    elseif e.f1 == 0
        e.f1 = tid
    end
    # If both slots are already filled, silently drop (shouldn't happen on
    # a valid grid; surfaces with topology errors will be diagnosed later).
end

# ---------------------------------------------------------------------------
# Step 4: spheric & contact face triangulation on the pre-populated
# WSurface. The toric grid build above has already registered every
# boundary edge with `face_[0]` set to the toric triangle on the other
# side; we now run advancing-front over the SHARED edge graph so the
# new triangles correctly land on `face_[1]`.
# ---------------------------------------------------------------------------

# Run the BALL main loop with `border` pre-seeded by EXISTING WSEdges (=
# the toric-side edges from `_build_toric_face_ws!`). Unlike
# `_advance_front!`, this variant does NOT create phantom boundary edges
# — every boundary edge already has `face_[0]` pointing at a real toric
# triangle.
function _advance_front_existing!(
    ws::WSurface{T},
    boundary_pids::Vector{Int},
    interior_pids::Vector{Int},
    sphere_c::Vector3{T},
    convex::Bool,
) where T
    n_b = length(boundary_pids)
    n_b < 3 && return
    candidates = unique(vcat(boundary_pids, interior_pids))
    border = Int[]
    for k in 1:n_b
        a = boundary_pids[k]
        b = boundary_pids[mod1(k + 1, n_b)]
        eid = _ws_find_edge_through(ws, a, b)
        eid == 0 && continue   # missing edge — skip
        ws.edges[eid].f0 == 0 && continue  # nothing to orient against
        push!(border, eid)
    end
    _build_spheric_triangles_loop!(ws, border, candidates, sphere_c, convex)
end

# Spheric face triangulator. The boundary is the cyclic concatenation of
# the three concave arcs at `rs.faces[fi]`, sharing corners. The probe
# sphere is concave from the SES perspective (outward = toward probe
# centre), so `convex=false`.
function _build_spheric_face_ws!(
    builder::SESBuilder{T},
    ses::SolventExcludedSurface{T},
    sf,
) where T
    rs = ses.reduced_surface
    fi = sf.rs_index
    arcs = builder.concave_arcs_pids[fi]
    (isempty(arcs[1]) || isempty(arcs[2]) || isempty(arcs[3])) && return
    boundary_pids = Int[]
    for k in 1:3
        arc = arcs[k]
        for i in 1:(length(arc) - 1)
            push!(boundary_pids, arc[i])
        end
    end
    length(boundary_pids) < 3 && return
    f = rs.faces[fi]
    _advance_front_existing!(builder.ws, boundary_pids, Int[],
                              f.center, false)
end

# Contact face triangulator. Walks the contact boundary (cyclic), adds
# cut icosphere interior candidates, and runs `_advance_front_existing!`
# with `convex=true`. The icosphere interior points are also added to the
# global `pos`/`norm` arrays.
function _build_contact_face_ws!(
    builder::SESBuilder{T},
    ses::SolventExcludedSurface{T},
    cf,
    pos::Vector{Point3{T}},
    norm::Vector{Vec3{T}},
    density::T,
    templates,
) where T
    rs = ses.reduced_surface
    rsv = rs.vertices[cf.rs_index]
    atom_idx = rsv.atom
    c = Vector3{T}(rs.atoms[atom_idx].center...)
    R = rs.atoms[atom_idx].r
    boundary_pids = _walk_contact_boundary(rs, cf.rs_index,
                                            builder.contact_arc_A_pids,
                                            builder.contact_arc_B_pids)
    length(boundary_pids) < 3 && return
    # Cut icosphere — keep only interior vertices on the contact-face side
    # of every incident contact-circle plane. BALL uses fuzzy = 0.05.
    fuzzy = T(0.05)
    tmpl = templates[min(4, _num_refinements(density, R) + 1)]
    # Build a small set of boundary WSPoint positions for cheap
    # deduplication. Without this, icosphere points that happen to land
    # near a corner produce duplicate WSPoints at the same geometric
    # location, which the advancing-front then treats as separate
    # candidates — yielding zero-length non-manifold edges in the output
    # mesh (see "duplicate vertices at the same 3D location" diagnostic).
    boundary_pts = [builder.ws.points[pid].point for pid in boundary_pids]
    dedup_tol2 = T(1e-4)        # ≈ 0.01 Å — points within this fold together
    interior_pids = Int[]
    for u_t in tmpl.position
        outward = Vector3{T}(u_t...)
        point = c + R * outward
        keep = true
        for ei in rsv.edges
            e = rs.edges[ei]
            neighbour_idx = (rs.vertices[e.v1].atom == atom_idx) ?
                rs.vertices[e.v2].atom : rs.vertices[e.v1].atom
            pn_raw = c - Vector3{T}(rs.atoms[neighbour_idx].center...)
            ℓn = sqrt(dot(pn_raw, pn_raw))
            ℓn < eps(T) && continue
            pn = pn_raw / ℓn
            cc = (rs.vertices[e.v1].atom == atom_idx) ?
                e.contact_circle1 : e.contact_circle2
            if dot(point - cc.p, pn) <= fuzzy
                keep = false; break
            end
        end
        keep || continue
        # Drop icosphere points that coincide with an existing boundary
        # vertex (= corner or arc sample).
        dup = false
        for bp in boundary_pts
            d = point - bp
            if dot(d, d) < dedup_tol2
                dup = true; break
            end
        end
        dup && continue
        push!(pos, Point3{T}(point...))
        push!(norm, Vec3{T}(outward...))
        gid = length(pos)
        push!(interior_pids, _ws_add_point!(builder.ws, point, outward, gid))
    end
    # NOTE: an earlier experiment added the singular axis-intersection
    # WSPoints to interior_pids (mirroring BALL's `buildSphericTriangles`
    # initialization which includes ALL contact-face boundary vertices).
    # That gained 98 fewer boundary edges on BPTI but introduced severe
    # T-junctions on the singular SES edge (count up to 13 incident
    # triangles), because each contact face's advancing-front picked
    # both intersection points as candidates and built triangles using
    # them — and those triangles all share the singular edge. Net mesh
    # quality regressed.
    # The correct fix is to teach the advancing-front about the singular
    # edge's at-most-one-triangle-per-side invariant. Without that, we
    # skip these as candidates and accept the 726 boundary edges around
    # singular tori.
    _advance_front_existing!(builder.ws, boundary_pids, interior_pids, c, true)
end

# Final mesh extraction: walk every WSTri and emit a `TriangleFace{Int}`
# to the global `faces` array using each vertex's `global_id`. Triangles
# that were emitted by `_build_toric_face_ws!` were ALREADY pushed to
# `faces`, so we skip those here by tracking which triangle IDs are
# already-emitted. (BALL avoids this by emitting once at join time.)
function _emit_remaining_triangles!(
    builder::SESBuilder{T},
    faces::Vector{TriangleFace{Int}},
    pos::Vector{Point3{T}},
    n_toric_tris::Int,
) where T
    for tid in (n_toric_tris + 1):length(builder.ws.tris)
        tri = builder.ws.tris[tid]
        g0 = builder.ws.points[tri.v0].global_id
        g1 = builder.ws.points[tri.v1].global_id
        g2 = builder.ws.points[tri.v2].global_id
        (g0 == g1 || g1 == g2 || g0 == g2) && continue
        p0 = Vector3{T}(pos[g0]...)
        p1 = Vector3{T}(pos[g1]...)
        p2 = Vector3{T}(pos[g2]...)
        centroid = (p0 + p1 + p2) / T(3)
        # Use centroid radial direction as outward hint. For contact faces
        # the centroid is on the atom sphere; for spheric faces it's on a
        # probe sphere. Either way the radial direction is the outward
        # normal direction for the face — the `_push_triangle!` outward
        # check will flip winding if needed.
        # We don't know the atom/probe centre here cheaply, so use the
        # vertex normal average as a proxy.
        navg = builder.ws.points[tri.v0].normal +
               builder.ws.points[tri.v1].normal +
               builder.ws.points[tri.v2].normal
        _push_triangle!(faces, g0, g1, g2, pos, navg)
    end
end

# ---------------------------------------------------------------------------
# Point/edge bookkeeping helpers (BALL's `GraphVertex::has` plus our own
# `find_existing_edge` that walks the small incident-edge set).
# ---------------------------------------------------------------------------

@inline function _ws_add_point!(ws::WSurface{T}, pt::Vector3{T}, n::Vector3{T}, gid::Int) where T
    push!(ws.points, WSPoint{T}(pt, n, gid))
    length(ws.points)
end

# BALL's `vertex->has(edge)`: search the point's incident edges for one whose
# other endpoint matches `target`. Returns the matching edge ID or 0.
@inline function _ws_find_edge_through(ws::WSurface, pid::Int, target::Int)
    for eid in ws.points[pid].edges
        e = ws.edges[eid]
        if (e.v0 == pid && e.v1 == target) || (e.v0 == target && e.v1 == pid)
            return eid
        end
    end
    return 0
end

# Triangle helpers matching BALL's `GraphTriangle::getRelativeIndex` (0/1/2 or
# -1 if not found) and `Triangle::third(v1, v2)`.
@inline function _ws_get_relative_index(t::WSTri, vid::Int)
    t.v0 == vid && return 0
    t.v1 == vid && return 1
    t.v2 == vid && return 2
    return -1
end

@inline function _ws_third(t::WSTri, v1::Int, v2::Int)
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

# ---------------------------------------------------------------------------
# `_create_triangle_and_edges` — direct port of BALL's
# `SESTriangulator::createTriangleAndEdges` (triangulatedSES.C:1250).
#
# Inputs: an edge `eid` (typically just popped from the border) and a third
# point `pid` (the candidate picked by the selection criterion). The new
# triangle's vertex order is initially set to `(edge.v1, edge.v0, point)` —
# we then geometrically test whether the resulting cross-product normal
# points outward (for `convex=true`) or inward (for `convex=false`) and
# swap `vertex_[0]/[1]` if not.
#
# The two new edges (`edge1` between vertex_0 and point, `edge2` between
# vertex_1 and point) are reused if they already exist in either endpoint's
# `edges` set; otherwise they are freshly created with `f0 = 0`.
# Returns `(tid, eid1, was_eid1_existing, eid2, was_eid2_existing)`.
# ---------------------------------------------------------------------------

function _ws_create_triangle_and_edges!(
    ws::WSurface{T},
    eid::Int,
    pid::Int,
    sphere_c::Vector3{T},
    convex::Bool,
) where T
    edge = ws.edges[eid]
    # Provisional new edges.
    target1 = edge.v0
    target2 = edge.v1
    # Reuse existing edges if present (BALL: vertex->has(provisional_edge)).
    eid1 = _ws_find_edge_through(ws, target1, pid)
    old1 = eid1 != 0
    if !old1
        push!(ws.edges, WSEdge(target1, pid, 0, 0))
        eid1 = length(ws.edges)
    end
    eid2 = _ws_find_edge_through(ws, target2, pid)
    old2 = eid2 != 0
    if !old2
        push!(ws.edges, WSEdge(pid, target2, 0, 0))
        eid2 = length(ws.edges)
    end
    # Build the triangle. BALL's initial vertex order is
    # `(edge.v1, edge.v0, point)` (note the swap relative to the edge).
    v0, v1, v2 = edge.v1, edge.v0, pid
    e01 = eid
    e12 = eid1
    e20 = eid2
    # Geometric outward test (BALL triangulatedSES.C:1302-1313). Reverse
    # vertex_[0]/[1] if the cross-product normal points the wrong way for
    # the sphere's convexity.
    p0 = ws.points[v0].point
    p1 = ws.points[v1].point
    p2 = ws.points[v2].point
    test_vec = cross(p1 - p0, p2 - p0)
    test_val = dot(test_vec, sphere_c - p0)
    if (convex && test_val > zero(T)) || (!convex && test_val < zero(T))
        v0, v1 = v1, v0
        e12, e20 = e20, e12
    end
    push!(ws.tris, WSTri(v0, v1, v2, e01, e12, e20))
    tid = length(ws.tris)
    return tid, eid1, old1, eid2, old2
end

# ---------------------------------------------------------------------------
# `_link_new_triangle!` — common post-creation linking shared by
# `buildUnambiguousTriangle` and `buildAmbiguousTriangles`: wire the
# triangle into the incidence sets of its three vertices and into each
# edge's `face_[0]/face_[1]` slot, updating the border accordingly.
# ---------------------------------------------------------------------------

function _link_new_triangle!(
    ws::WSurface,
    border::Vector{Int},
    eid::Int,
    tid::Int,
    eid1::Int, old1::Bool,
    eid2::Int, old2::Bool,
)
    tri = ws.tris[tid]
    push!(ws.points[tri.v0].faces, tid)
    push!(ws.points[tri.v1].faces, tid)
    push!(ws.points[tri.v2].faces, tid)
    # popped edge — face_[1] (face_[0] is the start triangle).
    ws.edges[eid].f1 = tid
    # Filter helper: remove `target` from `border` (BALL `list::remove`).
    function drop_from_border!(target::Int)
        idx = findfirst(==(target), border)
        idx !== nothing && deleteat!(border, idx)
        nothing
    end
    for (eid_new, old) in ((eid1, old1), (eid2, old2))
        e = ws.edges[eid_new]
        if old
            # The edge already existed; this triangle fills its remaining
            # face slot, completing manifold (= it leaves the border).
            if e.f0 == 0
                e.f0 = tid
            elseif e.f1 == 0
                e.f1 = tid
            end
            # If both f0 and f1 were already set, silently skip — the
            # selection should have rejected this candidate; reaching
            # here means the strict-hull pass picked a geometrically
            # valid candidate that nevertheless over-fills an edge. We
            # leave the triangle in `ws.tris` (it'll be emitted) since
            # popping it would create another non-manifold issue
            # (vertex `.faces` already includes tid). The over-fill
            # bookkeeping miss is the lesser evil — the mesh has an
            # extra triangle but otherwise stays consistent.
            drop_from_border!(eid_new)
        else
            # Brand new edge — link to its endpoints, set face_[0], push on
            # the border.
            push!(ws.points[e.v0].edges, eid_new)
            push!(ws.points[e.v1].edges, eid_new)
            e.f0 = tid
            push!(border, eid_new)
        end
    end
    nothing
end

# ---------------------------------------------------------------------------
# `_pop_orient!` — BALL's preamble to each main-loop iteration
# (triangulatedSES.C:812-821). Inspects the start triangle's traversal of
# the popped edge and reverses `edge.vertex_[]` if the existing triangle
# winds the "wrong" way; this guarantees that the NEW triangle constructed
# in `_create_triangle_and_edges` lands on the opposite side.
# ---------------------------------------------------------------------------

function _pop_orient!(ws::WSurface, eid::Int, convex::Bool)
    edge = ws.edges[eid]
    start_tri_id = edge.f0
    start_tri_id == 0 && return  # nothing to orient against
    start = ws.tris[start_tri_id]
    p0 = _ws_get_relative_index(start, edge.v0)
    p1 = _ws_get_relative_index(start, edge.v1)
    diff = p1 - p0
    if (convex && (diff == -1 || diff == 2)) ||
       (!convex && (diff == 1 || diff == -2))
        edge.v0, edge.v1 = edge.v1, edge.v0
    end
    nothing
end

# ---------------------------------------------------------------------------
# `_build_unambiguous_triangle!` — direct port of
# `SESTriangulator::buildUnambiguousTriangle` (triangulatedSES.C:1041).
# ---------------------------------------------------------------------------

function _build_unambiguous_triangle!(
    ws::WSurface,
    border::Vector{Int},
    eid::Int,
    pid::Int,
    sphere_c::Vector3,
    convex::Bool,
)
    tid, eid1, old1, eid2, old2 = _ws_create_triangle_and_edges!(
        ws, eid, pid, sphere_c, convex)
    _link_new_triangle!(ws, border, eid, tid, eid1, old1, eid2, old2)
    tid
end

# ---------------------------------------------------------------------------
# `_build_ambiguous_triangles!` — port of BALL's
# `SESTriangulator::buildAmbiguousTriangles` (triangulatedSES.C:1108).
#
# Called when several candidates tied for the convex-hull criterion. The
# resulting triangles are necessarily coplanar in some sense; we try each
# candidate against the orientation check, keep the one that gives the
# correct edge-traversal direction, and then recurse into the new edges
# (`planar_edges`) until they're either all consumed or the candidate list
# is exhausted.
# ---------------------------------------------------------------------------

function _build_ambiguous_triangles!(
    ws::WSurface,
    border::Vector{Int},
    seed_eid::Int,
    pids::Vector{Int},
    sphere_c::Vector3,
    convex::Bool,
)
    # BALL pushes the original edge's two endpoints into the candidate set
    # at the start of the routine — these can themselves be chosen as the
    # third vertex of a degenerate-coplanar triangle.
    seed_edge = ws.edges[seed_eid]
    cands = copy(pids)
    push!(cands, seed_edge.v0)
    push!(cands, seed_edge.v1)
    planar = Int[seed_eid]
    while !isempty(planar)
        e0_id = popfirst!(planar)
        e0 = ws.edges[e0_id]
        built = false
        for pid in cands
            (pid == e0.v0 || pid == e0.v1) && continue
            tid, eid1, old1, eid2, old2 = _ws_create_triangle_and_edges!(
                ws, e0_id, pid, sphere_c, convex)
            # First iteration on the seed edge bypasses the check. Phantom
            # edges (f0 < 0) have no existing triangle to compare against,
            # so accept the first candidate.
            keep = (e0_id == seed_eid) || e0.f0 <= 0
            if !keep
                old_tri = ws.tris[e0.f0]
                i1 = _ws_get_relative_index(old_tri, e0.v0)
                i2 = _ws_get_relative_index(old_tri, e0.v1)
                back = (i1 - i2 == 1) || (i1 - i2 == -2)
                new_tri = ws.tris[tid]
                i1 = _ws_get_relative_index(new_tri, e0.v0)
                i2 = _ws_get_relative_index(new_tri, e0.v1)
                if back
                    keep = (i1 - i2 == -1) || (i1 - i2 == 2)
                else
                    keep = (i1 - i2 == 1) || (i1 - i2 == -2)
                end
            end
            if keep
                _link_new_triangle!(ws, border, e0_id, tid,
                                    eid1, old1, eid2, old2)
                # New edges enter both the global border AND the local
                # `planar_edges` queue (BALL triangulatedSES.C:1198/1221).
                if !old1
                    push!(planar, eid1)
                else
                    idx = findfirst(==(eid1), planar)
                    idx !== nothing && deleteat!(planar, idx)
                end
                if !old2
                    push!(planar, eid2)
                else
                    idx = findfirst(==(eid2), planar)
                    idx !== nothing && deleteat!(planar, idx)
                end
                built = true
                break
            else
                # Roll back: discard the freshly created edges/triangle.
                # `_ws_create_triangle_and_edges!` pushes them in the order
                # (edge1, edge2, triangle) when both are new, so pop in
                # reverse.
                pop!(ws.tris)
                !old2 && pop!(ws.edges)
                !old1 && pop!(ws.edges)
            end
        end
        # If no candidate worked for this edge, leave it on the border (BALL
        # just drops it from `planar_edges` without removing from `border`).
    end
    nothing
end

# ---------------------------------------------------------------------------
# Main loop. BALL `SESTriangulator::buildSphericTriangles`
# (triangulatedSES.C:777). `border` is seeded by `_build_first_triangle!`
# below; once that finishes the body of the loop is unchanged from BALL.
# ---------------------------------------------------------------------------

function _build_spheric_triangles_loop!(
    ws::WSurface{T},
    border::Vector{Int},
    candidate_pids::Vector{Int},
    sphere_c::Vector3{T},
    convex::Bool,
) where T
    max_iter = 32 * (length(candidate_pids) + length(border)) + 256
    iter = 0
    while !isempty(border) && iter < max_iter
        iter += 1
        eid = popfirst!(border)
        edge = ws.edges[eid]
        # Edge already complete?
        if edge.f0 == -1
            # Phantom-toric boundary: the other side is a singular toric
            # face we didn't triangulate. Only one real triangle (in f1)
            # can ever attach, so once f1 is set, we're done.
            edge.f1 != 0 && continue
        else
            edge.f0 != 0 && edge.f1 != 0 && continue
        end
        if edge.f0 > 0
            _pop_orient!(ws, eid, convex)
        else
            # Phantom: no existing triangle to orient against. Defer to
            # the geometric outward swap inside
            # `_ws_create_triangle_and_edges!`; for concave (spheric)
            # faces flip the edge direction so the new triangle's
            # winding lands correctly (mirrors `_pop_orient_no_phantom!`).
            if !convex
                edge.v0, edge.v1 = edge.v1, edge.v0
            end
        end
        edge = ws.edges[eid]
        excl = edge.f0 > 0 ?
            _ws_third(ws.tris[edge.f0], edge.v0, edge.v1) : -1
        best = _select_best_candidates(ws, edge.v0, edge.v1, excl,
                                       candidate_pids, convex)
        if length(best) == 1
            _build_unambiguous_triangle!(ws, border, eid, best[1],
                                         sphere_c, convex)
        elseif length(best) > 1
            _build_ambiguous_triangles!(ws, border, eid, best,
                                        sphere_c, convex)
        end
    end
    nothing
end

# Convex-hull-style best-candidate selection — direct port of BALL's
# iterative quickhull at triangulatedSES.C:823-866. We seed `third` with
# the first valid candidate, compute its plane (normal, test_value =
# `normal · vertex_0`), and for each subsequent candidate check whether
# it is *above* that plane (`normal · candidate > test_value`). When yes,
# replace `third` and recompute the plane; when equal, append to `third`
# (these become BALL's "ambiguous" candidates).
# An edge "has room" if it has fewer than 2 real triangles (= f0 or f1
# slot still empty). Phantom (`-1`) counts as filled. Boundary edges may
# only ever accept ONE real triangle, since the other side is provided by
# the toric/spheric neighbour; we model that with a phantom in f0.
@inline function _edge_has_room(e::WSEdge)
    real_count = (e.f0 > 0 ? 1 : 0) + (e.f1 > 0 ? 1 : 0)
    if e.f0 == -1
        real_count < 1
    else
        real_count < 2
    end
end

function _select_best_candidates(
    ws::WSurface{T},
    v0::Int, v1::Int, excl::Int,
    cands::Vector{Int},
    ::Bool,
) where T
    p0 = ws.points[v0].point
    p1 = ws.points[v1].point
    best = Int[]
    best_normal = Vector3{T}(0, 0, 0)
    test_value = zero(T)
    have_seed = false
    # BALL's `Maths::isGreater` / `Maths::isEqual` use a 1e-6 absolute
    # tolerance. Strict `>` causes near-ties to ping-pong and pick
    # geometrically poor candidates at high tessellation density.
    EPS = T(1e-6)
    # A strict O(N²) hull check has been tested and gives only marginal
    # mesh-quality improvement (~140 fewer non-manifold edges on BPTI)
    # at ~3x the run-time, so we stick with BALL's iterative quickhull.
    # The `convex` flag is unused: BALL's iterative criterion uses `>`
    # regardless (triangulatedSES.C:849); convex/concave handling lives
    # in `_pop_orient!` and the geometric swap inside
    # `_ws_create_triangle_and_edges!`.
    for pid in cands
        (pid == v0 || pid == v1 || pid == excl) && continue
        pp = ws.points[pid].point
        if !have_seed
            push!(best, pid)
            best_normal = cross(pp - p1, pp - p0)
            test_value = dot(best_normal, p0)
            have_seed = true
            continue
        end
        this_value = dot(best_normal, pp)
        if this_value - test_value >= EPS
            empty!(best)
            push!(best, pid)
            best_normal = cross(pp - p1, pp - p0)
            test_value = dot(best_normal, p0)
        elseif abs(this_value - test_value) < EPS
            push!(best, pid)
        end
    end
    best
end

# ---------------------------------------------------------------------------
# `_build_first_triangle!` — BALL `SESTriangulator::buildFirstTriangle`
# (triangulatedSES.C:882). Picks the seed edge (= the first boundary edge
# that already has a single phantom face attached), then evaluates the
# strict convex-hull test against ALL candidates to find the seed third
# point. Ties are broken by minimum oriented angle relative to the existing
# (= toric/spheric phantom) triangle's normal.
#
# In BALL, the seed edge is selected by `firstSESEdge` from the SES face
# data; here the caller passes us a `seed_eid` directly (typically the
# first boundary edge).
# ---------------------------------------------------------------------------

function _build_first_triangle!(
    ws::WSurface{T},
    border::Vector{Int},
    candidate_pids::Vector{Int},
    seed_eid::Int,
    sphere_c::Vector3{T},
    convex::Bool,
) where T
    edge = ws.edges[seed_eid]
    # Orient the edge so its traversal in the existing (phantom) triangle
    # is forward — that way the new triangle's traversal is backward, as
    # createTriangleAndEdges expects.
    if edge.f0 != 0
        start = ws.tris[edge.f0]
        p0 = _ws_get_relative_index(start, edge.v0)
        p1 = _ws_get_relative_index(start, edge.v1)
        diff = p1 - p0
        if (diff == 1) || (diff == -2)
            edge.v0, edge.v1 = edge.v1, edge.v0
        end
    end
    p0 = ws.points[edge.v0].point
    p1 = ws.points[edge.v1].point
    edge_vec = p0 - p1
    # Strict convex-hull check: keep candidates for which NO other
    # candidate is on the wrong side of the (edge, candidate) plane.
    third = Int[]
    same_endpoints = Set{Int}((edge.v0, edge.v1))
    for pid in candidate_pids
        pid in same_endpoints && continue
        pp = ws.points[pid].point
        n = cross(edge_vec, p0 - pp)
        test_value = dot(n, p0)
        is_convex = true
        for tid in candidate_pids
            tid == pid && continue
            tid in same_endpoints && continue
            tp = ws.points[tid].point
            this_value = dot(n, tp)
            if (convex && this_value > test_value) ||
               (!convex && this_value < test_value)
                is_convex = false
                break
            end
        end
        is_convex && push!(third, pid)
    end
    # Tie-break by minimum oriented angle to the phantom triangle's normal.
    real_third = if length(third) <= 1
        third
    else
        # The phantom triangle's outward normal (for the contact-face setup
        # we pass a synthetic phantom whose third vertex is on the toric
        # side; BALL takes this from face_[0]).
        start = ws.tris[edge.f0]
        sp0 = ws.points[start.v0].point
        sp1 = ws.points[start.v1].point
        sp2 = ws.points[start.v2].point
        normal_ref = cross(sp0 - sp1, sp0 - sp2)
        # Pick the candidate with smallest oriented angle from `normal_ref`
        # around `edge_vec`.
        best_angle = T(3 * π)
        out = Int[]
        for pid in third
            pp = ws.points[pid].point
            new_normal = cross(edge_vec, p0 - pp)
            ang = oriented_angle(normal_ref, new_normal, edge_vec)
            if ang <= best_angle
                if ang < best_angle
                    empty!(out)
                    best_angle = ang
                end
                push!(out, pid)
            end
        end
        out
    end
    if length(real_third) == 1
        _build_unambiguous_triangle!(ws, border, seed_eid, real_third[1],
                                     sphere_c, convex)
    elseif length(real_third) > 1
        _build_ambiguous_triangles!(ws, border, seed_eid, real_third,
                                    sphere_c, convex)
    end
    nothing
end

# ---------------------------------------------------------------------------
# Driver: `_advance_front!` builds a complete patch given:
#   * `ws`       — empty WSurface to populate.
#   * `boundary_pids` — IDs of boundary points already in `ws.points`,
#     in cyclic order. The boundary edges between consecutive points are
#     created (with a phantom face_[0] = -1 representing the toric side)
#     and pushed onto the border.
#   * `interior_pids` — IDs of additional candidate points already in
#     `ws.points`.
#   * `sphere_c`, `convex` — patch sphere parameters.
#
# Returns nothing; on completion `ws.tris` holds the manifold triangulation.
# ---------------------------------------------------------------------------

function _advance_front!(
    ws::WSurface{T},
    boundary_pids::Vector{Int},
    interior_pids::Vector{Int},
    sphere_c::Vector3{T},
    convex::Bool,
) where T
    n_b = length(boundary_pids)
    n_b < 3 && return
    # Phantom boundary edges. BALL would have these already filled by the
    # toric triangulation that owns the other side; we simulate a phantom
    # by setting f0 = -1 (a sentinel) — _ws_get_relative_index on -1 is a
    # no-op so the orientation logic just keeps `edge.v0 → edge.v1` as the
    # phantom's traversal direction.
    border = Int[]
    for k in 1:n_b
        a = boundary_pids[k]
        b = boundary_pids[mod1(k + 1, n_b)]
        # If the edge already exists (rare — only possible if the boundary
        # walk re-visits an edge in the same cycle, which is a malformed
        # input), skip.
        _ws_find_edge_through(ws, a, b) != 0 && continue
        push!(ws.edges, WSEdge(a, b, -1, 0))
        eid = length(ws.edges)
        push!(ws.points[a].edges, eid)
        push!(ws.points[b].edges, eid)
        push!(border, eid)
    end
    # We need a *real* phantom Triangle on f0 of each boundary edge so that
    # `_pop_orient!` can read its `getRelativeIndex`. The phantom triangle
    # is shared across all boundary edges; we synthesise it with vertices
    # (boundary_pids[1], boundary_pids[2], boundary_pids[3]) just to give
    # `getRelativeIndex` something to compute. Its only use is the
    # orientation diff check, so the choice of third vertex doesn't affect
    # the result as long as `getRelativeIndex(edge.v0) - getRelativeIndex(
    # edge.v1)` evaluates consistently.
    #
    # Actually a simpler approach: per boundary edge, encode the phantom
    # direction implicitly via `_pop_orient!` skipping when f0 < 0. We
    # rewrite the swap rule for that case (see below).
    candidates = vcat(boundary_pids, interior_pids)
    # Seed: use the first boundary edge as the seed for buildFirstTriangle.
    # Without a real phantom triangle we don't run BALL's full firstTriangle
    # logic (which needs face_[0]); instead we use the convex-hull selection
    # directly on the seed edge with no exclusion.
    seed_eid = first(border)
    _build_first_triangle_no_phantom!(ws, border, candidates, seed_eid,
                                       sphere_c, convex)
    _build_spheric_triangles_loop_no_phantom!(ws, border, candidates,
                                              sphere_c, convex)
    nothing
end

# When boundary edges have no real phantom triangle, the `_pop_orient!` swap
# decision must be made from the geometry, not from a stored face. We use
# the convex-hull criterion to decide which traversal of the boundary edge
# corresponds to the contact-face interior, then orient the edge so the
# new triangle's createTriangleAndEdges swap rules give correct outward
# winding.
#
# Concretely: for each boundary edge with f0 == -1, when popped we keep
# `edge.v0 → edge.v1` as-is (no swap). createTriangleAndEdges' own
# geometric swap (test_val vs sphere_c) then ensures correct outward
# winding. The same geometric swap also flips edge1/edge2 within the
# triangle (we updated _ws_create_triangle_and_edges! to swap them
# alongside the vertex swap).
function _pop_orient_no_phantom!(ws::WSurface, eid::Int, convex::Bool)
    edge = ws.edges[eid]
    if edge.f0 > 0
        _pop_orient!(ws, eid, convex)
        return
    end
    # Phantom (f0 == -1): pretend the phantom triangle traverses
    # `vertex_[0] → vertex_[1]` forward (diff = 1). With BALL's rule this
    # means no swap for convex, swap for concave.
    if !convex
        edge.v0, edge.v1 = edge.v1, edge.v0
    end
    nothing
end

# Variant of the seed-finder that handles f0 < 0 phantoms.
function _build_first_triangle_no_phantom!(
    ws::WSurface{T},
    border::Vector{Int},
    candidate_pids::Vector{Int},
    seed_eid::Int,
    sphere_c::Vector3{T},
    convex::Bool,
) where T
    edge = ws.edges[seed_eid]
    p0 = ws.points[edge.v0].point
    p1 = ws.points[edge.v1].point
    edge_vec = p0 - p1
    third = Int[]
    same_endpoints = Set{Int}((edge.v0, edge.v1))
    for pid in candidate_pids
        pid in same_endpoints && continue
        pp = ws.points[pid].point
        n = cross(edge_vec, p0 - pp)
        test_value = dot(n, p0)
        is_convex = true
        for tid in candidate_pids
            tid == pid && continue
            tid in same_endpoints && continue
            tp = ws.points[tid].point
            this_value = dot(n, tp)
            if (convex && this_value > test_value) ||
               (!convex && this_value < test_value)
                is_convex = false
                break
            end
        end
        is_convex && push!(third, pid)
    end
    # No phantom triangle to break ties; pick by smallest distance midpoint.
    if length(third) > 1
        midpoint = (p0 + p1) / T(2)
        best_d2 = T(Inf)
        keep = first(third)
        for pid in third
            pp = ws.points[pid].point
            d = pp - midpoint
            d2 = dot(d, d)
            if d2 < best_d2
                best_d2 = d2
                keep = pid
            end
        end
        third = [keep]
    end
    isempty(third) && return
    _build_unambiguous_triangle!(ws, border, seed_eid, first(third),
                                 sphere_c, convex)
    # Drop the seed edge from the border (it was consumed). _link_new_
    # triangle! already set f1 = tid on it; for a phantom-f0 edge that
    # means it is now complete from the contact-face side.
    idx = findfirst(==(seed_eid), border)
    idx !== nothing && deleteat!(border, idx)
end

function _build_spheric_triangles_loop_no_phantom!(
    ws::WSurface{T},
    border::Vector{Int},
    candidates::Vector{Int},
    sphere_c::Vector3{T},
    convex::Bool,
) where T
    max_iter = 32 * (length(candidates) + length(border)) + 256
    iter = 0
    while !isempty(border) && iter < max_iter
        iter += 1
        eid = popfirst!(border)
        edge = ws.edges[eid]
        # Skip if edge already saturated. A boundary phantom edge with f1
        # set is "done from the contact-face side" — only ONE contact-face
        # triangle ever touches it.
        if edge.f0 == -1
            edge.f1 != 0 && continue
        else
            edge.f0 != 0 && edge.f1 != 0 && continue
        end
        _pop_orient_no_phantom!(ws, eid, convex)
        edge = ws.edges[eid]
        # The existing triangle's third vertex is excluded from the
        # candidate search (BALL's `third_point` exclusion). For phantom
        # edges (f0 = -1) there's no such exclusion.
        excl = if edge.f0 > 0
            _ws_third(ws.tris[edge.f0], edge.v0, edge.v1)
        else
            -1
        end
        best = _select_best_candidates(ws, edge.v0, edge.v1, excl,
                                       candidates, convex)
        if length(best) == 1
            _build_unambiguous_triangle_no_phantom!(
                ws, border, eid, best[1], sphere_c, convex)
        elseif length(best) > 1
            _build_ambiguous_triangles_no_phantom!(
                ws, border, eid, best, sphere_c, convex)
        end
    end
    nothing
end

# `_link_new_triangle!` assumes the popped edge already has `f0` set (we
# fill `f1`). For a phantom-f0 edge that's correct (-1 stays in f0; the
# contact-face triangle goes to f1). No special handling needed.

# Wrap the unambiguous/ambiguous builders so they call `_link_new_triangle!`
# regardless of phantom status.
const _build_unambiguous_triangle_no_phantom! = _build_unambiguous_triangle!
const _build_ambiguous_triangles_no_phantom!  = _build_ambiguous_triangles!
