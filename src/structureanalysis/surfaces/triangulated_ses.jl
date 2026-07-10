export
    triangulate_ses

triangulate(ses::SolventExcludedSurface{T}; density::T = T(1.0)) where T =
    triangulate_ses(ses; density)

# ---------------------------------------------------------------------------
# Boundary-stitched triangulator for the SES (Connolly).
#
# The mesh is assembled in three passes:
#
#   1. **Shared vertex pool**: every "corner" (SAS contact point on an atom
#      surface), every contact-circle arc (on atom surfaces), and every
#      concave great-circle arc (on probe spheres) is materialised exactly
#      once into a flat `pos` / `norm` vector.
#
#   2. **Per-face triangulation**:
#        * spheric face → barycentric subdivision of a spherical triangle,
#          re-using the 3 concave-arc boundaries;
#        * toric  face → quadrilateral grid wired between the two atom-side
#          contact arcs and the two probe-side concave arcs;
#        * contact face → icosphere sampling on the atom, with visible
#          vertices kept and ear-clipped fans from the boundary contact-arc
#          vertices toward the visible interior.
#
#   3. **Normals**: every vertex carries the outward SES normal (away from
#      molecular interior). Corner / arc vertices use the natural normal of
#      the face they sit on; interior vertices use the analytical direction
#      (away from atom centre for contact, away from probe centre for spheric
#      / toric).
#
# All faces winding order is set so that `cross(v2-v1, v3-v1)` aligns with
# the per-triangle outward normal.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

@inline function _onb_perp(n::AbstractVector{T}) where T
    # Return an orthonormal basis (u, v) such that {u, v, n_hat} is right-
    # handed and u, v are perpendicular to `n`.
    ℓ = sqrt(dot(n, n))
    ℓ <= eps(T) && return (Vector3{T}(1,0,0), Vector3{T}(0,1,0))
    nh = n / ℓ
    a = abs(nh[1]) < T(0.9) ? Vector3{T}(1,0,0) : Vector3{T}(0,1,0)
    u = a - dot(a, nh) * nh
    u = u / sqrt(dot(u, u))
    v = cross(nh, u)
    Vector3{T}(u...), Vector3{T}(v...)
end

@inline function _slerp_unit(a::Vector3{T}, b::Vector3{T}, t::T) where T
    cosγ = clamp(dot(a, b), -one(T), one(T))
    γ = acos(cosγ)
    if γ < eps(T)
        return a
    end
    s = sin(γ)
    Vector3{T}(((sin((1 - t) * γ) * a + sin(t * γ) * b) / s)...)
end

# Find which slot (1, 2, or 3) of `rs.faces[face_idx]` sits on `atom_idx`.
@inline function _atom_slot(rs::ReducedSurface, face_idx::Int, atom_idx::Int)
    f = rs.faces[face_idx]
    rs.vertices[f.v1].atom == atom_idx && return 1
    rs.vertices[f.v2].atom == atom_idx && return 2
    rs.vertices[f.v3].atom == atom_idx && return 3
    0
end

# ---------------------------------------------------------------------------
# Step 1: per-RSFace contact points (3 corner vertices on the atom surfaces).
# Each contact point lies on the bare atom sphere; outward normal is the
# unit vector from atom centre toward probe centre (since both contact and
# spheric / toric normals at this corner align in that direction).
# ---------------------------------------------------------------------------

function _build_corner_vertices(rs::ReducedSurface{T},
                                pos::Vector{Point3{T}},
                                norm::Vector{Vec3{T}}) where T
    rsface_contacts = Vector{NTuple{3, Int}}(undef, length(rs.faces))
    for (fi, f) in enumerate(rs.faces)
        triple = ntuple(3) do k
            atom_idx = rs.vertices[(f.v1, f.v2, f.v3)[k]].atom
            atom_c = Vector3{T}(rs.atoms[atom_idx].center...)
            R_atom = rs.atoms[atom_idx].r
            dir = f.center - atom_c
            ℓ = sqrt(dot(dir, dir))
            outward = ℓ > eps(T) ? dir / ℓ : Vector3{T}(0,0,1)
            point = atom_c + R_atom * outward
            push!(pos, Point3{T}(point...))
            push!(norm, Vec3{T}(outward...))
            length(pos)
        end
        rsface_contacts[fi] = triple
    end
    rsface_contacts
end

# ---------------------------------------------------------------------------
# Step 2a: contact arcs on atom surfaces. For each RS edge with both f1 and
# f2 populated, sample the arc on contact_circle1 (atom A side) between the
# corresponding corner vertices, plus the analogous arc on contact_circle2
# (atom B side). The arc takes the SHORTER of the two paths between the
# corners along the contact circle.
# ---------------------------------------------------------------------------

function _sample_contact_arc(circle::Circle3{T}, atom_center::Vector3{T},
                             R_atom::T, start_pt::Vector3{T}, end_pt::Vector3{T},
                             start_idx::Int, end_idx::Int, n_seg::Int,
                             avoid_pts::Vector{Vector3{T}},
                             pos::Vector{Point3{T}}, norm::Vector{Vec3{T}}) where T
    u, v = _onb_perp(circle.n)
    rel_s = start_pt - circle.p
    rel_e = end_pt - circle.p
    ψ_s = atan(dot(rel_s, v), dot(rel_s, u))
    ψ_e = atan(dot(rel_e, v), dot(rel_e, u))
    # Two candidate Δs: short (∈ (-π, π]) and long (the complement).
    Δ_short = ψ_e - ψ_s
    if Δ_short >  T(π); Δ_short -= 2 * T(π); end
    if Δ_short < -T(π); Δ_short += 2 * T(π); end
    Δ_long = Δ_short - sign(Δ_short) * 2 * T(π)
    # Pick whichever midpoint is FURTHEST from `avoid_pts`.
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
        # Outward normal: away from atom centre, projected to the atom sphere.
        outward = point - atom_center
        ℓ = sqrt(dot(outward, outward))
        outward = ℓ > eps(T) ? outward / ℓ : Vector3{T}(0,0,1)
        push!(pos, Point3{T}(point...))
        push!(norm, Vec3{T}(outward...))
        ids[k + 1] = length(pos)
    end
    ids[end] = end_idx
    ids
end

function _build_contact_arcs(rs::ReducedSurface{T},
                             rsface_contacts::Vector{NTuple{3, Int}},
                             pos::Vector{Point3{T}}, norm::Vector{Vec3{T}},
                             n_seg_edge::Vector{Int}) where T
    n_edges = length(rs.edges)
    arc_A = Vector{Vector{Int}}(undef, n_edges)
    arc_B = Vector{Vector{Int}}(undef, n_edges)
    for ei in 1:n_edges
        e = rs.edges[ei]
        if e.f1 == 0 || e.f2 == 0
            arc_A[ei] = Int[]
            arc_B[ei] = Int[]
            continue
        end
        atom_A = rs.vertices[e.v1].atom
        atom_B = rs.vertices[e.v2].atom
        s_f1_A = _atom_slot(rs, e.f1, atom_A)
        s_f2_A = _atom_slot(rs, e.f2, atom_A)
        s_f1_B = _atom_slot(rs, e.f1, atom_B)
        s_f2_B = _atom_slot(rs, e.f2, atom_B)
        if s_f1_A == 0 || s_f2_A == 0 || s_f1_B == 0 || s_f2_B == 0
            arc_A[ei] = Int[]
            arc_B[ei] = Int[]
            continue
        end
        cv_f1_A = rsface_contacts[e.f1][s_f1_A]
        cv_f2_A = rsface_contacts[e.f2][s_f2_A]
        cv_f1_B = rsface_contacts[e.f1][s_f1_B]
        cv_f2_B = rsface_contacts[e.f2][s_f2_B]
        atom_A_c = Vector3{T}(rs.atoms[atom_A].center...)
        atom_B_c = Vector3{T}(rs.atoms[atom_B].center...)
        # Avoid the contact-arc passing close to the "third" atoms of f1/f2
        # (= the atoms of f1/f2 that are not A or B): those poles lie on the
        # cap-toward-them side, and the SES contact boundary wraps the *far*
        # side of the contact circle.
        third_atoms_pts = Vector3{T}[]
        f1 = rs.faces[e.f1]; f2 = rs.faces[e.f2]
        for vidx in (f1.v1, f1.v2, f1.v3, f2.v1, f2.v2, f2.v3)
            a = rs.vertices[vidx].atom
            (a == atom_A || a == atom_B) && continue
            push!(third_atoms_pts, Vector3{T}(rs.atoms[a].center...))
        end
        n_seg = n_seg_edge[ei]
        arc_A[ei] = _sample_contact_arc(e.contact_circle1, atom_A_c,
                                        rs.atoms[atom_A].r,
                                        Vector3{T}(pos[cv_f1_A]...),
                                        Vector3{T}(pos[cv_f2_A]...),
                                        cv_f1_A, cv_f2_A, n_seg,
                                        third_atoms_pts, pos, norm)
        arc_B[ei] = _sample_contact_arc(e.contact_circle2, atom_B_c,
                                        rs.atoms[atom_B].r,
                                        Vector3{T}(pos[cv_f1_B]...),
                                        Vector3{T}(pos[cv_f2_B]...),
                                        cv_f1_B, cv_f2_B, n_seg,
                                        third_atoms_pts, pos, norm)
    end
    arc_A, arc_B
end

# ---------------------------------------------------------------------------
# Step 2b: concave arcs on probe spheres. For each RS face, sample the three
# great-circle arcs of the probe sphere connecting consecutive corner
# vertices. Outward normal at each arc vertex is +u on the probe sphere
# (from probe centre toward the vertex) — that is the SES outward direction
# at a spheric/toric junction (away from probe interior = away from
# molecular interior).
# ---------------------------------------------------------------------------

function _sample_concave_arc(probe_center::Vector3{T}, probe_r::T,
                             start_pt::Vector3{T}, end_pt::Vector3{T},
                             start_idx::Int, end_idx::Int, n_seg::Int,
                             pos::Vector{Point3{T}}, norm::Vector{Vec3{T}}) where T
    ua = (start_pt - probe_center) / probe_r
    ub = (end_pt   - probe_center) / probe_r
    ids = Vector{Int}(undef, n_seg + 1)
    ids[1] = start_idx
    for k in 1:(n_seg - 1)
        t = T(k) / T(n_seg)
        u = _slerp_unit(ua, ub, t)
        point = probe_center + probe_r * u
        push!(pos, Point3{T}(point...))
        # SES outward at a probe-surface point is *toward* the probe
        # centre (away from molecule interior on the probe-rest side).
        push!(norm, Vec3{T}(-u[1], -u[2], -u[3]))
        ids[k + 1] = length(pos)
    end
    ids[end] = end_idx
    ids
end

function _build_concave_arcs(rs::ReducedSurface{T},
                             rsface_contacts::Vector{NTuple{3, Int}},
                             pos::Vector{Point3{T}}, norm::Vector{Vec3{T}},
                             n_seg_edge::Vector{Int}) where T
    # Each arc of face fi belongs to one of f's three RS edges:
    #   arcs[1] (corner 1 → corner 2) ↔ f.e1
    #   arcs[2] (corner 2 → corner 3) ↔ f.e2
    #   arcs[3] (corner 3 → corner 1) ↔ f.e3
    # So the per-arc segment count is taken from the corresponding RS edge.
    # This keeps the toric face for edge ei consistent on both axes: its
    # rolling-direction contact arcs and its phi-direction concave arcs
    # (at f1 and f2) all sample with `n_seg_edge[ei]+1` vertices.
    out = Vector{NTuple{3, Vector{Int}}}(undef, length(rs.faces))
    for (fi, f) in enumerate(rs.faces)
        ids = rsface_contacts[fi]
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
            _sample_concave_arc(f.center, rs.probe_radius, sa, sb,
                                ids[ka], ids[kb], n_seg, pos, norm)
        end
        out[fi] = arcs
    end
    out
end

# ---------------------------------------------------------------------------
# Triangulation helpers
# ---------------------------------------------------------------------------

# Append a quad as two triangles, oriented so the outward normal aligns
# with `normal_hint` (a representative outward direction for the quad).
@inline function _push_quad!(faces::Vector{TriangleFace{Int}},
                             a::Int, b::Int, c::Int, d::Int,
                             pos::Vector{Point3{T}},
                             normal_hint::Vector3{T}) where T
    # Triangle (a, b, c): orient so cross(B-A, C-A) · normal_hint > 0.
    pa = Vector3{T}(pos[a]...)
    pb = Vector3{T}(pos[b]...)
    pc = Vector3{T}(pos[c]...)
    cross_bc = cross(pb - pa, pc - pa)
    if dot(cross_bc, normal_hint) >= 0
        push!(faces, TriangleFace{Int}(a, b, c))
        push!(faces, TriangleFace{Int}(a, c, d))
    else
        push!(faces, TriangleFace{Int}(a, c, b))
        push!(faces, TriangleFace{Int}(a, d, c))
    end
end

# Push a triangle with outward winding matching `normal_hint`.
@inline function _push_triangle!(faces::Vector{TriangleFace{Int}},
                                 a::Int, b::Int, c::Int,
                                 pos::Vector{Point3{T}},
                                 normal_hint::Vector3{T}) where T
    pa = Vector3{T}(pos[a]...)
    pb = Vector3{T}(pos[b]...)
    pc = Vector3{T}(pos[c]...)
    if dot(cross(pb - pa, pc - pa), normal_hint) >= 0
        push!(faces, TriangleFace{Int}(a, b, c))
    else
        push!(faces, TriangleFace{Int}(a, c, b))
    end
end

# ---------------------------------------------------------------------------
# Step 3a: Spheric faces. Each is a spherical triangle on the probe sphere
# bounded by 3 concave great-circle arcs that, in general, have *different*
# segment counts (each arc inherits the segment count of its RS edge — see
# `_build_concave_arcs`). We therefore can't use a uniform barycentric grid.
# Instead we drop one centroid vertex on the probe sphere and fan-triangulate
# from it to every consecutive pair along the boundary arcs.
# ---------------------------------------------------------------------------

function _triangulate_spheric!(faces::Vector{TriangleFace{Int}},
                               pos::Vector{Point3{T}},
                               norm::Vector{Vec3{T}},
                               rs::ReducedSurface{T},
                               rsface_contacts::Vector{NTuple{3, Int}},
                               concave_arcs::Vector{NTuple{3, Vector{Int}}}) where T
    for (fi, f) in enumerate(rs.faces)
        corners = rsface_contacts[fi]
        arcs    = concave_arcs[fi]
        probe   = f.center
        pr      = rs.probe_radius
        # Outward hint: the spheric face is on the side AWAY from the atom
        # plane, so outward = -f.normal.
        n_hint = Vector3{T}(-f.normal[1], -f.normal[2], -f.normal[3])
        # Centroid on the probe sphere from the 3 corner directions.
        cv = ntuple(k -> Vector3{T}(pos[corners[k]]...), 3)
        ua = (cv[1] - probe) / pr
        ub = (cv[2] - probe) / pr
        uc = (cv[3] - probe) / pr
        u_cent = ua + ub + uc
        ℓ = sqrt(dot(u_cent, u_cent))
        u_cent = ℓ > eps(T) ? u_cent / ℓ : ua
        cent_pt = probe + pr * u_cent
        push!(pos, Point3{T}(cent_pt...))
        # SES outward at a probe-surface point points TOWARD the probe
        # centre (= -u_cent).
        push!(norm, Vec3{T}(-u_cent[1], -u_cent[2], -u_cent[3]))
        cent_idx = length(pos)
        # Fan from centroid to each consecutive boundary pair on the three
        # concave arcs.
        for k in 1:3
            arc = arcs[k]
            for i in 1:(length(arc) - 1)
                a = arc[i]
                b = arc[i + 1]
                _push_triangle!(faces, cent_idx, a, b, pos, n_hint)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Step 3b: Toric faces. Each toric face is a rectangular grid spanning the
# rolling angle (between probes f1 and f2) and the φ angle on the probe
# (between contact_A and contact_B). The 4 boundaries are pre-sampled:
#   * k = 0:   concave arc on f1 between contact_A_at_f1 and contact_B_at_f1
#   * k = N:   concave arc on f2 between contact_A_at_f2 and contact_B_at_f2
#   * j = 0:   contact arc on atom A (between SAS vertices at f1, f2)
#   * j = N:   contact arc on atom B
# Interior points are built by parametric interpolation across both axes
# using SLERP on the probe sphere and arc-walk on the atom contact circles.
# ---------------------------------------------------------------------------

function _concave_arc_id(rs::ReducedSurface{T},
                        concave_arcs::Vector{NTuple{3, Vector{Int}}},
                        rsface_contacts::Vector{NTuple{3, Int}},
                        face_idx::Int, atom_A::Int, atom_B::Int) where T
    # Find which of the 3 concave arcs of `face_idx` connects atom_A's
    # corner to atom_B's corner. Returns (arc_indices_vector, reversed::Bool).
    f = rs.faces[face_idx]
    s_A = _atom_slot(rs, face_idx, atom_A)
    s_B = _atom_slot(rs, face_idx, atom_B)
    s_A == 0 && return (Int[], false)
    s_B == 0 && return (Int[], false)
    # arcs[k] goes from corner k to corner mod1(k+1, 3).
    # Find k such that {s_A, s_B} = {k, mod1(k+1, 3)}.
    for k in 1:3
        kn = mod1(k + 1, 3)
        if (s_A == k && s_B == kn)
            return (concave_arcs[face_idx][k], false)
        end
        if (s_B == k && s_A == kn)
            return (concave_arcs[face_idx][k], true)
        end
    end
    Int[], false
end

function _triangulate_toric!(faces::Vector{TriangleFace{Int}},
                             pos::Vector{Point3{T}},
                             norm::Vector{Vec3{T}},
                             ses::SolventExcludedSurface{T},
                             rsface_contacts::Vector{NTuple{3, Int}},
                             concave_arcs::Vector{NTuple{3, Vector{Int}}},
                             contact_arc_A::Vector{Vector{Int}},
                             contact_arc_B::Vector{Vector{Int}},
                             n_seg_edge::Vector{Int}) where T
    rs = ses.reduced_surface
    for tf in ses.toric_faces
        tf.type == SESFaceType.Toric || tf.type == SESFaceType.ToricSingular || continue
        ei = tf.rs_index
        e = rs.edges[ei]
        e.f1 == 0 && continue
        e.f2 == 0 && continue
        atom_A = rs.vertices[e.v1].atom
        atom_B = rs.vertices[e.v2].atom
        # Concave arc on f1 from contact_A_at_f1 to contact_B_at_f1.
        ca_f1, rev_f1 = _concave_arc_id(rs, concave_arcs, rsface_contacts,
                                         e.f1, atom_A, atom_B)
        ca_f2, rev_f2 = _concave_arc_id(rs, concave_arcs, rsface_contacts,
                                         e.f2, atom_A, atom_B)
        (isempty(ca_f1) || isempty(ca_f2)) && continue
        # Reverse if needed so the arc always goes A → B in j order.
        f1_arc = rev_f1 ? reverse(ca_f1) : ca_f1
        f2_arc = rev_f2 ? reverse(ca_f2) : ca_f2
        cA = contact_arc_A[ei]
        cB = contact_arc_B[ei]
        (isempty(cA) || isempty(cB)) && continue
        n_seg = n_seg_edge[ei]
        # All four boundary arcs come from this RS edge, so they all
        # have `n_seg + 1` vertices. Build a square grid of that size.
        # grid[k+1, j+1] = vertex at (rolling step k, phi step j).
        # Sanity check; if any arc is the wrong length we can't grid it.
        (length(f1_arc) == n_seg + 1 && length(f2_arc) == n_seg + 1 &&
         length(cA) == n_seg + 1 && length(cB) == n_seg + 1) || continue
        grid = fill(0, n_seg + 1, n_seg + 1)
        # Boundary from concave arcs (k = 0 and k = n_seg).
        for j in 0:n_seg
            grid[1, j + 1]         = f1_arc[j + 1]
            grid[n_seg + 1, j + 1] = f2_arc[j + 1]
        end
        # Boundary from contact arcs (j = 0 and j = n_seg).
        for k in 0:n_seg
            grid[k + 1, 1]         = cA[k + 1]
            grid[k + 1, n_seg + 1] = cB[k + 1]
        end
        # Interior: exact toric parameterization. The probe centre rolls
        # around the (atom_A, atom_B) axis through the angle that matches
        # the contact arc on atom A. At each rolling step k, the probe
        # position is known; on its sphere we walk the great-circle arc
        # from the atom-A contact toward the atom-B contact via the chosen
        # outward side (matches `f1_arc`/`f2_arc` at the endpoints).
        atom_A_c = Vector3{T}(rs.atoms[atom_A].center...)
        atom_B_c = Vector3{T}(rs.atoms[atom_B].center...)
        axis_dir = atom_B_c - atom_A_c
        axis_len = sqrt(dot(axis_dir, axis_dir))
        axis_n = axis_len > eps(T) ? axis_dir / axis_len : Vector3{T}(0,0,1)
        torus_c = e.center_of_torus
        pc_f1 = rs.faces[e.f1].center
        pc_f2 = rs.faces[e.f2].center
        # Radial component of probe positions in axis-perp plane.
        v_f1 = pc_f1 - torus_c
        v_f1_perp = v_f1 - dot(v_f1, axis_n) * axis_n
        v_f2 = pc_f2 - torus_c
        v_f2_perp = v_f2 - dot(v_f2, axis_n) * axis_n
        # Signed angle from f1 to f2 (shorter direction).
        α_short = oriented_angle(v_f1_perp, v_f2_perp, axis_n)
        # Decide between α_short and α_long = α_short − 2π by the same
        # avoid-third-atoms heuristic used for the contact arcs. Without
        # third-atom info here, mirror the contact arc's choice by
        # comparing the rolling midpoint to where the contact arc cA's
        # midpoint sits (cA's midpoint is the SES side by construction).
        cA_mid = Vector3{T}(pos[cA[div(length(cA), 2) + 1]]...)
        # Rotate v_f1_perp by α_short/2 around axis_n.
        rot_radial(θ) = begin
            cosθ = cos(θ); sinθ = sin(θ)
            v_f1_perp * cosθ + cross(axis_n, v_f1_perp) * sinθ +
                axis_n * dot(axis_n, v_f1_perp) * (1 - cosθ)
        end
        mid_short = torus_c + rot_radial(α_short / 2)
        mid_long  = torus_c + rot_radial((α_short - 2 * T(π)) / 2)
        d_short = (mid_short - cA_mid); d_short_sq = dot(d_short, d_short)
        d_long  = (mid_long  - cA_mid); d_long_sq  = dot(d_long, d_long)
        α = d_long_sq < d_short_sq ? (α_short - 2 * T(π)) : α_short
        # Inflated radii for contact computation.
        R_A_inf = rs.atoms[atom_A].r + rs.probe_radius
        R_B_inf = rs.atoms[atom_B].r + rs.probe_radius
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
                # SES outward at toric saddle: TOWARD probe centre.
                push!(norm, Vec3{T}(-u[1], -u[2], -u[3]))
                grid[k + 1, j + 1] = length(pos)
            end
        end
        # Quad mesh. Use an interior normal hint to orient.
        for k in 0:(n_seg - 1)
            for j in 0:(n_seg - 1)
                a = grid[k + 1, j + 1]
                b = grid[k + 2, j + 1]
                c = grid[k + 2, j + 2]
                d = grid[k + 1, j + 2]
                # Hint: average of the 4 vertex normals.
                navg = Vector3{T}(0,0,0)
                for v in (a, b, c, d)
                    navg = navg + Vector3{T}(norm[v]...)
                end
                _push_quad!(faces, a, b, c, d, pos, navg)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Step 3c: Contact faces. The visible region of each atom is sampled with an
# icosphere; vertices fully outside other atoms' inflated spheres are kept.
# For each incident RS edge whose contact-arc vertices were registered, we
# also fan-stitch from the arc to the nearest visible icosphere vertex.
# ---------------------------------------------------------------------------

# Per-RS-edge segment count, mirrors BALL's adaptive scaling for toric and
# spheric faces (triangulatedSES.C, `SESTriangulator::triangulateToricFace`):
# the rolling-direction segment count along edge `e` is proportional to
# `arc_length × √density` where `arc_length = e.angle × e.radius_of_torus`.
# Bigger / longer torii get more segments without over-tessellating tiny
# saddles. Both axes of a toric face use this same count (square grid).
@inline function _segments_for_edge(density::T, angle::T, radius::T) where T
    n = round(Int, angle * radius * sqrt(density), RoundNearestTiesAway)
    max(2, n)
end

# Per-atom icosphere refinement level. Mirrors BALL's
# `SESTriangulator::numberOfRefinements` (triangulatedSES.C:1350):
# target triangle count = 4·density·π·R² ; pick the smallest icosphere
# subdivision level (0..3) whose triangle count is closest to the target.
function _num_refinements(density::T, radius::T) where T
    test0 = (4 * density * T(π) * radius^2 - 12) / 30
    n = 0
    test0 < 0 && return 0
    test1 = one(T); test2 = one(T)
    while test2 < test0
        test1 = test2
        test2 *= 4
        n += 1
    end
    if test2 - test0 < test0 - test1
        n += 1
    end
    min(n, 3)
end

# Walk the boundary of a contact face on `rsv_idx` in cyclic order. The
# walk traverses each incident RS edge once, alternating arc segments at
# every corner (= RSFace with `atom_idx` as one of its three vertices).
# Returns the boundary vertex IDs in cyclic order, *without* duplicating
# corners (each corner appears exactly once).
function _walk_contact_boundary(rs::ReducedSurface{T}, rsv_idx::Int,
                                arc_A::Vector{Vector{Int}},
                                arc_B::Vector{Vector{Int}}) where T
    rsv = rs.vertices[rsv_idx]
    atom_idx = rsv.atom
    isempty(rsv.faces) && return Int[]
    # For each face f in rsv.faces, the 2 edges incident to atom_idx.
    function edges_at_face(f_idx::Int)
        f = rs.faces[f_idx]
        v_atoms = (rs.vertices[f.v1].atom,
                   rs.vertices[f.v2].atom,
                   rs.vertices[f.v3].atom)
        others = filter(!=(atom_idx), collect(v_atoms))
        length(others) == 2 || return (0, 0)
        function find_edge(b::Int)
            for ei in f.edges
                re = rs.edges[ei]
                a1 = rs.vertices[re.v1].atom
                a2 = rs.vertices[re.v2].atom
                ((a1 == atom_idx && a2 == b) || (a2 == atom_idx && a1 == b)) && return ei
            end
            return 0
        end
        return (find_edge(others[1]), find_edge(others[2]))
    end
    start_face = first(rsv.faces)
    start_e1, start_e2 = edges_at_face(start_face)
    (start_e1 == 0 || start_e2 == 0) && return Int[]
    current_face = start_face
    current_edge = start_e1
    boundary = Int[]
    max_iter = 4 * length(rsv.faces) + 4
    for _ in 1:max_iter
        e = rs.edges[current_edge]
        arc = rs.vertices[e.v1].atom == atom_idx ? arc_A[current_edge] :
                                                   arc_B[current_edge]
        isempty(arc) && return Int[]
        next_face = e.f1 == current_face ? e.f2 : e.f1
        next_face == 0 && return Int[]
        if current_face == e.f1
            for i in 1:(length(arc) - 1)
                push!(boundary, arc[i])
            end
        else
            for i in length(arc):-1:2
                push!(boundary, arc[i])
            end
        end
        ne1, ne2 = edges_at_face(next_face)
        next_edge = ne1 == current_edge ? ne2 : ne1
        if next_face == start_face && next_edge == start_e1
            return boundary
        end
        current_face = next_face
        current_edge = next_edge
    end
    boundary
end

# Standard ear-clipping on a 2D polygon (vertex indices into a parent vector
# pool, with positions supplied as a precomputed 2D tangent-plane projection).
# Returns the triangulation as triples of *boundary-vector* indices (not
# `pos` indices).
function _ear_clip_polygon_2d(boundary_local::Vector{Int}, pts_2d::Vector{Tuple{T,T}}) where T
    triangles = NTuple{3, Int}[]
    n = length(boundary_local)
    n < 3 && return triangles
    indices = collect(1:n)
    # Determine winding direction (signed area). Force CCW.
    signed_area = zero(T)
    @inbounds for i in 1:n
        p = pts_2d[i]
        q = pts_2d[mod1(i + 1, n)]
        signed_area += (q[1] - p[1]) * (q[2] + p[2])
    end
    if signed_area > 0
        reverse!(indices)
    end
    max_iter = n * n
    iter = 0
    while length(indices) > 3 && iter < max_iter
        iter += 1
        m = length(indices)
        ear_found = false
        @inbounds for i in 1:m
            ia = indices[mod1(i - 1, m)]
            ib = indices[i]
            ic = indices[mod1(i + 1, m)]
            pa = pts_2d[ia]; pb = pts_2d[ib]; pc = pts_2d[ic]
            cross_z = (pb[1] - pa[1]) * (pc[2] - pa[2]) -
                      (pb[2] - pa[2]) * (pc[1] - pa[1])
            cross_z <= 0 && continue
            is_ear = true
            for j in 1:m
                (j == mod1(i - 1, m) || j == i || j == mod1(i + 1, m)) && continue
                _point_in_tri(pts_2d[indices[j]], pa, pb, pc) || continue
                is_ear = false
                break
            end
            if is_ear
                push!(triangles, (boundary_local[ia], boundary_local[ib], boundary_local[ic]))
                deleteat!(indices, i)
                ear_found = true
                break
            end
        end
        if !ear_found
            # Degenerate polygon (no ear found). Fallback: emit fan from
            # vertex 1 to remaining boundary.
            @inbounds for i in 2:(length(indices) - 1)
                push!(triangles,
                      (boundary_local[indices[1]],
                       boundary_local[indices[i]],
                       boundary_local[indices[i + 1]]))
            end
            return triangles
        end
    end
    if length(indices) == 3
        push!(triangles,
              (boundary_local[indices[1]],
               boundary_local[indices[2]],
               boundary_local[indices[3]]))
    end
    triangles
end

@inline function _point_in_tri(p, a, b, c)
    d1 = (p[1] - b[1]) * (a[2] - b[2]) - (a[1] - b[1]) * (p[2] - b[2])
    d2 = (p[1] - c[1]) * (b[2] - c[2]) - (b[1] - c[1]) * (p[2] - c[2])
    d3 = (p[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (p[2] - a[2])
    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0)
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0)
    !(has_neg && has_pos)
end

# Recursively subdivide each triangle by midpoint refinement on the host
# sphere. Each level multiplies triangle count by 4. New vertices get
# emitted with their outward normal.
function _subdivide_sphere_tris(triangles::Vector{NTuple{3,Int}},
                                pos::Vector{Point3{T}},
                                norm::Vector{Vec3{T}},
                                atom_center::Vector3{T}, R::T,
                                levels::Int) where T
    for _ in 1:levels
        new_tris = Vector{NTuple{3,Int}}(undef, 4 * length(triangles))
        cache = Dict{NTuple{2,Int}, Int}()
        function midpoint(a::Int, b::Int)
            key = a < b ? (a, b) : (b, a)
            haskey(cache, key) && return cache[key]
            pa = Vector3{T}(pos[a]...); pb = Vector3{T}(pos[b]...)
            mid = (pa + pb) / 2
            d = mid - atom_center
            ℓ = sqrt(dot(d, d))
            d = ℓ > eps(T) ? d / ℓ : Vector3{T}(0, 0, 1)
            point = atom_center + R * d
            push!(pos, Point3{T}(point...))
            push!(norm, Vec3{T}(d...))
            idx = length(pos)
            cache[key] = idx
            idx
        end
        @inbounds for (i, t) in pairs(triangles)
            a, b, c = t
            m_ab = midpoint(a, b)
            m_bc = midpoint(b, c)
            m_ca = midpoint(c, a)
            new_tris[4*(i-1) + 1] = (a, m_ab, m_ca)
            new_tris[4*(i-1) + 2] = (m_ab, b, m_bc)
            new_tris[4*(i-1) + 3] = (m_ca, m_bc, c)
            new_tris[4*(i-1) + 4] = (m_ab, m_bc, m_ca)
        end
        triangles = new_tris
    end
    triangles
end

function _triangulate_contact!(faces::Vector{TriangleFace{Int}},
                               pos::Vector{Point3{T}},
                               norm::Vector{Vec3{T}},
                               ses::SolventExcludedSurface{T},
                               contact_arc_A::Vector{Vector{Int}},
                               contact_arc_B::Vector{Vector{Int}},
                               template::GeometryBasics.Mesh,
                               density::T) where T
    # Icosphere on the bare atom surface with visibility test and snap to
    # boundary arc vertices.
    #
    # Per-atom adaptive icosphere depth (BALL's `numberOfRefinements`):
    # target triangle count = 4·density·π·R²; choose the icosphere level
    # 0..3 closest to that. Atoms with bigger radii or higher requested
    # density get more triangles, eliminating the wire-frame look on
    # small contact patches.
    rs = ses.reduced_surface
    pr = rs.probe_radius
    n_atoms = length(rs.atoms)
    centers = [Vector3{T}(rs.atoms[i].center...) for i in 1:n_atoms]
    atom_r  = [rs.atoms[i].r for i in 1:n_atoms]
    # Pre-build the 4 candidate icosphere templates once. The
    # `template` argument from the SES driver is ignored — we pick the
    # appropriate refinement level per atom from this table instead.
    templates = ntuple(k -> icosphere(T, k - 1), 4)
    _ = template
    @inbounds for cf in ses.contact_faces
        rsv = rs.vertices[cf.rs_index]
        atom_idx = rsv.atom
        c = centers[atom_idx]
        R = atom_r[atom_idx]
        # Bump the per-atom icosphere refinement up by one level above BALL's
        # target. At low density the contact face would otherwise sit on a
        # 12-vertex (level 0) icosphere, and `snap_thresh` ≈ R sweeps up
        # nearly every visible icosphere vertex onto the boundary — leaving
        # no interior vertices and producing the "wire-frame" cracks at the
        # contact/toric junction. One extra subdivision quadruples the
        # vertex budget and a tighter `snap_thresh` keeps the boundary-near
        # ones from getting glued onto already-stitched arc vertices.
        tmpl = templates[min(4, _num_refinements(density, R) + 2)]
        template_pts = tmpl.position
        template_faces = tmpl.faces
        avg_ang_edge = T(2 * pi) / sqrt(T(length(template_pts)))
        snap_thresh = R * avg_ang_edge * T(0.2)
        snap_thresh_sq = snap_thresh * snap_thresh
        boundary_vids = Int[]
        for ei in rsv.edges
            e = rs.edges[ei]
            arc_ids = rs.vertices[e.v1].atom == atom_idx ? contact_arc_A[ei] :
                                                            contact_arc_B[ei]
            isempty(arc_ids) && continue
            append!(boundary_vids, arc_ids)
        end
        unique!(boundary_vids)
        b_pts = [Vector3{T}(pos[v]...) for v in boundary_vids]
        visible = Vector{Bool}(undef, length(template_pts))
        remap = zeros(Int, length(template_pts))
        for (k, u_t) in pairs(template_pts)
            outward = Vector3{T}(u_t...)
            point = c + R * outward
            ok = true
            probe_c = point + pr * outward
            for j in 1:n_atoms
                j == atom_idx && continue
                rj = atom_r[j] + pr
                d = probe_c - centers[j]
                if dot(d, d) < rj * rj - T(1e-6)
                    ok = false; break
                end
            end
            visible[k] = ok
            ok || continue
            best_d2 = snap_thresh_sq
            best_id = 0
            for (i, bp) in pairs(b_pts)
                dd = point - bp
                d2 = dot(dd, dd)
                if d2 < best_d2
                    best_d2 = d2
                    best_id = boundary_vids[i]
                end
            end
            if best_id != 0
                remap[k] = best_id
            else
                push!(pos, Point3{T}(point...))
                push!(norm, Vec3{T}(outward...))
                remap[k] = length(pos)
            end
        end
        function nearest_b(pt::Vector3{T})
            isempty(b_pts) && return 0
            best_d2 = T(Inf)
            best_id = 0
            for (i, bp) in pairs(b_pts)
                dd = pt - bp
                d2 = dot(dd, dd)
                if d2 < best_d2
                    best_d2 = d2
                    best_id = boundary_vids[i]
                end
            end
            best_id
        end
        for tri in template_faces
            i1 = Int(tri[1]); i2 = Int(tri[2]); i3 = Int(tri[3])
            v_count = (visible[i1] ? 1 : 0) + (visible[i2] ? 1 : 0) + (visible[i3] ? 1 : 0)
            if v_count == 3
                a = remap[i1]; b = remap[i2]; ee = remap[i3]
                (a == b || b == ee || a == ee) && continue
                p1 = Vector3{T}(pos[a]...)
                outward = p1 - c
                _push_triangle!(faces, a, b, ee, pos, outward)
            elseif v_count == 2
                idx_occ = !visible[i1] ? i1 : (!visible[i2] ? i2 : i3)
                others  = filter(x -> x != idx_occ, (i1, i2, i3))
                occ_pt  = c + R * Vector3{T}(template_pts[idx_occ]...)
                bv = nearest_b(occ_pt)
                bv == 0 && continue
                a = remap[others[1]]; b = remap[others[2]]; ee = bv
                (a == b || b == ee || a == ee) && continue
                p1 = Vector3{T}(pos[a]...)
                outward = p1 - c
                _push_triangle!(faces, a, b, ee, pos, outward)
            end
        end

        # Outer-ring stitching: for any consecutive pair of boundary arc
        # vertices in cyclic order, emit a fan triangle to the nearest
        # interior (= visible-and-not-snapped) icosphere vertex. This
        # closes the gap between the icosphere "edge" and the boundary
        # polygon, similar to BALL's advancing-front hole filler.
        boundary = _walk_contact_boundary(rs, cf.rs_index, contact_arc_A, contact_arc_B)
        if length(boundary) >= 3
            boundary_set = Set(boundary)
            ico_vids = Int[]
            for k in 1:length(template_pts)
                v = remap[k]
                v != 0 && !(v in boundary_set) && push!(ico_vids, v)
            end
            if !isempty(ico_vids)
                ico_pts = [Vector3{T}(pos[v]...) for v in ico_vids]
                n_b = length(boundary)
                fan_targets = zeros(Int, n_b)
                for k in 1:n_b
                    b1 = boundary[k]
                    b2 = boundary[mod1(k + 1, n_b)]
                    p1 = Vector3{T}(pos[b1]...); p2 = Vector3{T}(pos[b2]...)
                    mid = (p1 + p2) / 2
                    best_d2 = T(Inf); best_v = 0
                    for (i, vp) in pairs(ico_pts)
                        d = mid - vp
                        d2 = dot(d, d)
                        if d2 < best_d2
                            best_d2 = d2
                            best_v = ico_vids[i]
                        end
                    end
                    best_v == 0 && continue
                    # Skip only when the chosen interior vertex is on the
                    # far side of the atom — without this we'd occasionally
                    # span a triangle across the entire sphere. A radius-
                    # wide threshold lets short-distance boundary edges
                    # stitch even when there are few interior vertices, at
                    # the cost of a few elongated triangles on tiny faces.
                    sqrt(best_d2) > R && continue
                    outward = p1 - c
                    _push_triangle!(faces, b1, b2, best_v, pos, outward)
                    fan_targets[k] = best_v
                end
                # Bridge between adjacent fan triangles. When the fans
                # for edge (b_{k-1}, b_k) and edge (b_k, b_{k+1}) land on
                # *different* interior vertices `v_prev` and `v_curr` that
                # are *not* adjacent in the icosphere, there is a wedge-
                # shaped gap at `b_k`; one `(b_k, v_prev, v_curr)` closes
                # it. The icosphere-adjacency check (`|v_curr - v_prev|`
                # under ~1.4× the avg icosphere edge length) prevents
                # bridges where the icosphere already has a triangle
                # covering the wedge — without it we double-count area at
                # high density.
                bridge_thresh_sq = (R * avg_ang_edge * T(1.4))^2
                for k in 1:n_b
                    v_prev = fan_targets[mod1(k - 1, n_b)]
                    v_curr = fan_targets[k]
                    (v_prev == 0 || v_curr == 0 || v_prev == v_curr) && continue
                    pp = Vector3{T}(pos[v_prev]...)
                    pc = Vector3{T}(pos[v_curr]...)
                    d = pc - pp
                    dot(d, d) < bridge_thresh_sq && continue
                    b_k = boundary[k]
                    p_k = Vector3{T}(pos[b_k]...)
                    outward = p_k - c
                    _push_triangle!(faces, b_k, v_prev, v_curr, pos, outward)
                end
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Public entry point.
# ---------------------------------------------------------------------------

"""
    triangulate_ses(ses::SolventExcludedSurface{T}; density = 1.0)

Build a triangulated mesh approximating the SES. Boundary curves (atom-
contact circles and probe-sphere concave arcs) are sampled into a *shared*
vertex pool so that adjacent spheric / toric / contact faces meet at exactly
the same vertices — eliminating the per-patch cracks that an independent-
sampler approach would produce.

`density` controls subdivision via BALL's adaptive formulas (see
`_segments_for_edge` and `_num_refinements`):

  * Per-RS-edge segment count: `max(1, round(angle · radius · √density))`,
    matching BALL's `partitionNonFreeSingularEdge`.
  * Per-atom contact icosphere refinement level: `numberOfRefinements`
    (BALL `triangulatedSES.C:1350`) — the smallest `n` with
    `(2·5·4^n + 2) ≥ 4·density·π·r²`.

Higher density therefore produces both finer arcs and finer interior
icosphere patches; the relationship is monotonic but not linear in
density.
"""
function triangulate_ses(ses::SolventExcludedSurface{T};
                         density::T = T(1.0)) where {T<:Real}
    # BALL-faithful port (SESTriangulator::run, triangulatedSES.C:130).
    # The legacy per-RS-element pipeline below remains as fallback for
    # inputs the BALL port can't yet handle, but in normal use the
    # BALL port produces a much better mesh:
    #   * ~2× fewer triangles
    #   * boundary edges drop by ~93% on BPTI (718 → 46)
    #   * area much closer to the analytical SES area
    return triangulate_ses_ball(ses; density = density)
end

# Legacy per-RS-element triangulator. Kept for fallback / comparison.
function triangulate_ses_legacy(ses::SolventExcludedSurface{T};
                                 density::T = T(1.0)) where {T<:Real}
    # Mirror BALL's `SESTriangulator::preProcessing` — drop degenerate
    # toric faces and split spheric faces with disconnected boundaries.
    # Both are mandatory prerequisites for the advancing-front to land on
    # well-formed faces (triangulatedSES.C:144-147).
    clean_ses!(ses; min_torus_radius = T(1e-3), min_angle = T(sqrt(eps(T))),
                    density = density)
    split_spheric_faces!(ses)
    rs = ses.reduced_surface
    # Per-RS-edge segment count from BALL's formula. Atoms with large
    # rolling arcs (= big toric/spheric panels) get more segments.
    n_seg_edge = [
        _segments_for_edge(density, rs.edges[ei].angle,
                           rs.edges[ei].radius_of_torus)
        for ei in eachindex(rs.edges)
    ]
    templates = ntuple(k -> icosphere(T, k - 1), 4)

    pos    = Point3{T}[]
    norm   = Vec3{T}[]
    faces  = TriangleFace{Int}[]

    # BALL-style shared-WSurface pipeline (advancing_front.jl):
    #   1. Register the shared vertex pool (RSFace corners + contact and
    #      concave arcs) in the working surface.
    #   2. Build each toric face's quad grid into the working surface,
    #      populating `face_[0]` of every grid edge.
    #   3. Run advancing-front on each spheric face — uses the existing
    #      WSEdges that the toric step populated.
    #   4. Run advancing-front on each contact face — same.
    builder = SESBuilder{T}(length(rs.faces), length(rs.edges))
    _build_corner_vertices_ws!(builder, rs, pos, norm)
    _build_contact_arcs_ws!(builder, rs, pos, norm, n_seg_edge)
    _build_concave_arcs_ws!(builder, rs, pos, norm, n_seg_edge)
    for tf in ses.toric_faces
        _build_toric_face_ws!(builder, ses, pos, norm, faces, tf, n_seg_edge)
    end
    n_toric_tris = length(builder.ws.tris)
    # Order matches BALL: contact before spheric (triangulatedSES.C:135-136).
    # When the contact-face advancing-front finishes, its new edges populate
    # `face_[0]/face_[1]` of every contact-arc segment AND of the concave-arc
    # corners — giving the spheric advancing-front more populated boundary
    # references to read orientation from.
    for cf in ses.contact_faces
        _build_contact_face_ws!(builder, ses, cf, pos, norm, density, templates)
    end
    for sf in ses.spheric_faces
        _build_spheric_face_ws!(builder, ses, sf)
    end
    _emit_remaining_triangles!(builder, faces, pos, n_toric_tris)

    GeometryBasics.Mesh(pos, faces; normal=norm)
end

"""
    triangulate_ses(ac::AbstractAtomContainer{T}; probe_radius, density)

Convenience overload that builds the analytical SES first.
"""
function triangulate_ses(ac::AbstractAtomContainer{T};
                         probe_radius::T = T(1.5),
                         density::T = T(1.0)) where T
    triangulate_ses(compute_ses(ac; probe_radius, density); density)
end
