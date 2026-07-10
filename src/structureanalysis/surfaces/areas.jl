export
    sas_area,
    sas_volume,
    ses_area

# ---------------------------------------------------------------------------
# Surface area / volume for the analytical molecular surfaces.
#
# For SAS we delegate to the Eisenhaber/Argos algorithm (`NumericalSAS`),
# which is the same routine BALL invokes for `calculateSASArea` — closed-form
# per atom modulo icosphere sampling, and exact in the high-density limit.
#
# For SES we sum closed-form contributions per face: spherical triangles on
# the probe sphere (spheric faces), saddle regions via the surface-of-
# revolution formula (toric faces), and Eisenhaber-style sampling on each
# atom's bare sphere clipped to its contact face (contact faces). The contact
# integration reuses the same NumericalSAS infrastructure but with atom radii
# instead of inflated radii, and uses the SES-style visibility test (probe
# placed tangent to the atom must clear all other atoms).
# ---------------------------------------------------------------------------

# Spherical-triangle area on a sphere of radius `R` given three corners and
# the sphere center. Uses interior-angle sum minus π via tangent projections.
function _spherical_triangle_area(R::T,
                                  a::AbstractVector{T}, b::AbstractVector{T},
                                  c::AbstractVector{T},
                                  center::AbstractVector{T}) where {T<:Real}
    ua = (a - center) / R
    ub = (b - center) / R
    uc = (c - center) / R
    @inline tangent(p, q) = begin
        v = q - p * dot(p, q)
        n = sqrt(dot(v, v))
        iszero(n) ? v : v / n
    end
    αa = acos(clamp(dot(tangent(ua, ub), tangent(ua, uc)), -one(T), one(T)))
    αb = acos(clamp(dot(tangent(ub, ua), tangent(ub, uc)), -one(T), one(T)))
    αc = acos(clamp(dot(tangent(uc, ua), tangent(uc, ub)), -one(T), one(T)))
    R^2 * (αa + αb + αc - T(π))
end

"""
    sas_area(ac::AbstractAtomContainer{T};
             probe_radius = 1.5,
             number_of_points = 1600) -> T

Total solvent-accessible surface area (Å²). Uses the Eisenhaber/Argos
icosphere-sampling algorithm (see [`compute_numerical_sas`](@ref)); the
result converges to the analytical SAS area as `number_of_points` grows.
"""
function sas_area(ac::AbstractAtomContainer{T};
                  probe_radius::T = T(1.5),
                  number_of_points::Int = 1600) where T
    r = compute_numerical_sas(ac; probe_radius, number_of_points,
                              compute_area=true, compute_volume=false)
    r.total_area
end

"""
    sas_area(sas::SolventAccessibleSurface{T}; number_of_points = 1600) -> T

Total area of an already-constructed SAS. Re-evaluates the underlying atom
geometry numerically; the topology stored in `sas` is not used.
"""
sas_area(sas::SolventAccessibleSurface{T}; number_of_points::Int = 1600) where T =
    sas_area_from_spheres(sas.reduced_surface.atoms,
                          sas.reduced_surface.probe_radius;
                          number_of_points)

"""
    sas_volume(ac::AbstractAtomContainer{T};
               probe_radius = 1.5,
               number_of_points = 1600) -> T

Total volume enclosed by the SAS (Å³).
"""
function sas_volume(ac::AbstractAtomContainer{T};
                    probe_radius::T = T(1.5),
                    number_of_points::Int = 1600) where T
    r = compute_numerical_sas(ac; probe_radius, number_of_points,
                              compute_area=false, compute_volume=true)
    r.total_volume
end

function sas_area_from_spheres(atoms::Vector{Sphere{T}},
                               probe_radius::T;
                               number_of_points::Int = 1600) where T
    # Synthesize a fake atom container by sampling directly off the spheres.
    template = _tessellation_for(T, number_of_points)
    template_pts = template.position
    np = length(template_pts)
    unit_area = 4 * T(π) / np
    total = zero(T)
    centers = [Vector3{T}(s.center...) for s in atoms]
    pr = probe_radius
    inflated = [s.r + pr for s in atoms]
    inflated2 = inflated .^ 2
    @inbounds for i in eachindex(atoms)
        s = atoms[i]
        s.r > T(0.001) || continue
        R = inflated[i]
        c = centers[i]
        n_visible = 0
        for u in template_pts
            p = c + R * Vector3{T}(u...)
            ok = true
            for j in eachindex(atoms)
                j == i && continue
                d = p - centers[j]
                if dot(d, d) < inflated2[j] - T(1e-6)
                    ok = false
                    break
                end
            end
            ok && (n_visible += 1)
        end
        total += R^2 * unit_area * n_visible
    end
    total
end

# ---------------------------------------------------------------------------
# SES area: closed-form per-face sum.
# ---------------------------------------------------------------------------

# Spheric face: spherical triangle on the probe sphere.
function _ses_spheric_face_area(ses::SolventExcludedSurface{T},
                                face::SESFace{T}) where T
    length(face.vertices) >= 3 || return zero(T)
    a = ses.vertices[face.vertices[1]].point
    b = ses.vertices[face.vertices[2]].point
    c = ses.vertices[face.vertices[3]].point
    _spherical_triangle_area(face.sphere.r, a, b, c, face.sphere.center)
end

# Toric face: surface of revolution between two contact angles on the probe
# sphere, swept by the rolling angle `α`.
#
#   A_toric = α · r · ( R · Δφ + r · (sin φ₂ − sin φ₁) )
#
# where R is the torus's major radius (axis → probe-center distance), r the
# probe radius, and φ₁,φ₂ are the polar angles of the two contact points on
# the probe sphere, measured from the torus plane.
#
# For *spindle* tori (R < r), the integrand `R + r·cos φ` is negative for
# |φ| > acos(-R/r); only the band |φ| ≤ φ_cut contributes to the real SES
# surface. We clip the integration to that band.
function _ses_toric_face_area(ses::SolventExcludedSurface{T},
                              face::SESFace{T}) where T
    rs = ses.reduced_surface
    # Sentinel toric face for a deleted RS edge — area = 0.
    face.rs_index == 0 && return zero(T)
    e  = rs.edges[face.rs_index]
    pr = rs.probe_radius
    R  = e.radius_of_torus
    α  = e.angle
    (R <= 0 || α <= 0) && return zero(T)
    axis = e.contact_circle1.n
    nax = sqrt(dot(axis, axis))
    nax <= eps(T) && return zero(T)
    axis = axis / nax
    h1 = dot(e.contact_circle1.p - e.center_of_torus, axis)
    h2 = dot(e.contact_circle2.p - e.center_of_torus, axis)
    sφ1 = clamp(h1 / pr, -one(T), one(T))
    sφ2 = clamp(h2 / pr, -one(T), one(T))
    φ1 = asin(sφ1)
    φ2 = asin(sφ2)
    if φ1 > φ2; φ1, φ2 = φ2, φ1; end
    # Spindle-torus clipping: only |φ| ≤ φ_cut survives.
    if R < pr
        φ_cut = acos(clamp(-R / pr, -one(T), one(T)))
        φ1 = max(φ1, -φ_cut)
        φ2 = min(φ2,  φ_cut)
        φ2 <= φ1 && return zero(T)
    end
    Δφ = φ2 - φ1
    A = α * pr * (R * Δφ + pr * (sin(φ2) - sin(φ1)))
    A < 0 ? zero(T) : A
end

# Contact face: numerically integrate the visible patch on each atom's bare
# sphere using a probe-tangent occlusion test (SES variant of the
# Eisenhaber/Argos sampling). Returns the per-atom area on the bare-radius
# sphere.
function _ses_contact_area_total(ses::SolventExcludedSurface{T},
                                 number_of_points::Int) where T
    rs = ses.reduced_surface
    pr = rs.probe_radius
    template = _tessellation_for(T, number_of_points)
    template_pts = template.position
    np = length(template_pts)
    unit_area = 4 * T(π) / np
    total = zero(T)
    centers = [Vector3{T}(rs.atoms[i].center...) for i in 1:length(rs.atoms)]
    atom_r  = [rs.atoms[i].r for i in 1:length(rs.atoms)]
    # An atom contributes only if it carries at least one RS vertex
    has_vertex = falses(length(rs.atoms))
    for v in rs.vertices
        v.atom == 0 && continue   # sentinel = deleted vertex
        has_vertex[v.atom] = true
    end
    @inbounds for i in eachindex(rs.atoms)
        has_vertex[i] || continue
        R = atom_r[i]
        R > T(0.001) || continue
        c = centers[i]
        n_visible = 0
        for u in template_pts
            p = c + R * Vector3{T}(u...)
            # Probe placed at p + pr·outward must not overlap any neighbour
            # inflated by pr.
            outward = Vector3{T}(u...)
            probe_center = p + pr * outward
            ok = true
            for j in eachindex(rs.atoms)
                j == i && continue
                rj = atom_r[j] + pr
                d = probe_center - centers[j]
                if dot(d, d) < rj * rj - T(1e-6)
                    ok = false
                    break
                end
            end
            ok && (n_visible += 1)
        end
        total += R^2 * unit_area * n_visible
    end
    total
end

"""
    ses_area(ses::SolventExcludedSurface{T}; number_of_points = 1600) -> T

Total solvent-excluded surface area (Å²) summed over the three face types:

 - **contact** (atom-sphere patches): Eisenhaber-style numerical integration
   using the SES-specific tangent-probe occlusion test;
 - **toric**: closed-form surface-of-revolution between the two probe
   contact angles, swept by the rolling angle;
 - **spheric**: spherical triangle area on the probe sphere.

`number_of_points` controls only the contact-face sampling; toric and
spheric contributions are computed in closed form. Singular toric faces
(probe touching ≥4 atoms) are included with their non-singular formula.
Spindle-torus regions (probe radius exceeding the torus major radius)
are clipped to the physically valid band by `_ses_toric_face_area`;
the analytical sum may still differ slightly from mesh integration in
densely packed structures due to the closed-form / numerical mixture.
Where exact agreement with mesh integration matters, use
`surface_area(triangulate_ses(ses; density))` instead.
"""
function ses_area(ses::SolventExcludedSurface{T};
                  number_of_points::Int = 1600) where {T<:Real}
    A = _ses_contact_area_total(ses, number_of_points)
    for f in ses.toric_faces
        A += _ses_toric_face_area(ses, f)
    end
    for f in ses.spheric_faces
        A += _ses_spheric_face_area(ses, f)
    end
    A
end

"""
    ses_area(ac::AbstractAtomContainer; probe_radius = 1.5, number_of_points = 1600)

Convenience overload that builds the analytical SES first.
"""
ses_area(ac::AbstractAtomContainer{T}; probe_radius::T = T(1.5),
         number_of_points::Int = 1600) where T =
    ses_area(compute_ses(ac; probe_radius); number_of_points)
