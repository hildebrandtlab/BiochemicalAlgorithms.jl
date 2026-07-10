export
    triangulate,
    triangulate_sas

"""
    triangulate(surface; density=1.0)

Dispatch entry point for surface triangulation. Routes to
[`triangulate_sas`](@ref) for a [`SolventAccessibleSurface`](@ref) and to
[`triangulate_ses`](@ref) for a [`SolventExcludedSurface`](@ref).
"""
triangulate(sas::SolventAccessibleSurface{T}; density::T = T(1.0)) where T =
    triangulate_sas(sas; density)

# ---------------------------------------------------------------------------
# BALL-style SAS triangulator. Mirrors `SASTriangulator::run` /
# `triangulateFace` (BALL/source/STRUCTURE/triangulatedSAS.C:111-186):
#
#   for each SAS face (= one per atom):
#     1. Build cutting planes from the face's SAS edges. The plane
#        passes through the edge's intersection-circle centre, oriented
#        outward via the face's `orientation` flag.
#     2. Pick the icosphere template by `numberOfRefinements(density,
#        radius)` and blow it up + shift to the atom's inflated sphere.
#     3. Tag each template TPoint as "inside" iff it falls behind ANY
#        cutting plane (= occluded by a neighbouring inflated atom).
#     4. Emit every template triangle whose three vertices are NOT
#        all inside (BALL's `removeInsideTriangles` keeps mixed-tag
#        triangles intact — they form the patch boundary).
#
# BALL intentionally leaves the SAS mesh "leaky" at neighbouring-atom
# intersection circles: there is no boundary stitching across atoms,
# and the un-tagged vertices stay in the mesh as orphans (BALL calls
# `deleteIsolatedPoints` after the loop; we just don't emit them).
#
# The earlier port had elaborate arc-sampling + snap-to-boundary
# logic that diverged from BALL and produced inconsistent meshes
# (BPTI/F64/d=2: 8636 tris, NM=168, bdy=2001, dup=35). This rewrite
# reproduces BALL's geometry exactly.
# ---------------------------------------------------------------------------

"""
    triangulate_sas(sas::SolventAccessibleSurface{T}; density = 1.0)

Triangulate the analytical SAS produced by [`compute_sas`](@ref).

Per-atom icosphere mesh of each inflated atom, with template vertices
that fall behind any neighbour's intersection plane removed (BALL's
`removeInsideTriangles`). The resulting mesh is patchwise — adjacent
atoms produce separate triangulated caps with no boundary stitching
between them, matching BALL's `TriangulatedSAS` semantics exactly.

`density` controls the icosphere refinement level via the same
`numberOfRefinements` formula used for the SES (see
[`triangulate_ses`](@ref)).
"""
function triangulate_sas(sas::SolventAccessibleSurface{T};
                         density::T = T(1.0)) where {T<:Real}
    pos    = Point3{T}[]
    norms  = Vec3{T}[]
    tris   = TriangleFace{Int}[]

    # Pre-build templates; SES triangulator's helper picks an integer
    # refinement level from (density, atom_radius).
    templates = ntuple(k -> icosphere(T, k - 1), 4)

    @inbounds for sf in sas.faces
        sf.rs_vertex == 0 && continue                  # sentinel = deleted
        c   = Vector3{T}(sf.sphere.center...)
        R   = sf.sphere.r

        # Step 1: cutting planes. BALL `createPlanes`
        # (triangulatedSAS.C:188): plane.p = edge.circle.p; plane.n =
        # orientation ? edge.circle.n : -edge.circle.n. We pre-shift the
        # plane test by the atom centre (`p - c`) to avoid recomputing
        # per template point.
        # Each plane stored as (normal, d) — d = n·p — since we only
        # ever compare `n · point` against `d`. BALL keeps both p and n
        # in TPlane3; we collapse to the precomputed distance.
        planes = Tuple{Vector3{T}, T}[]
        for (k, ei) in pairs(sf.edges)
            ei == 0 && continue
            e  = sas.edges[ei]
            nv = sf.orientation[k] ? e.circle.n : -e.circle.n
            push!(planes, (nv, dot(nv, e.circle.p)))
        end

        # Step 2: blow-up + shift the icosphere template.
        ref  = min(4, _num_refinements(density, R) + 1)
        tmpl = templates[ref]
        n_tmpl = length(tmpl.position)

        # Step 3: tag & emit. To match BALL's `removeInsideTriangles`
        # behaviour (drop only when ALL three vertices are tagged
        # inside), we still need to allocate TPoints for the kept
        # ones, so build a local remap table first and only push the
        # vertices we'll actually reference.
        inside = Vector{Bool}(undef, n_tmpl)
        for i in 1:n_tmpl
            u_t = tmpl.position[i]
            point = c + R * Vector3{T}(u_t...)
            tag = false
            for (nv, d) in planes
                # BALL: `Maths::isLessOrEqual(plane.n * point, plane.d)`
                # — point is inside iff its signed distance to the
                # plane (n · p) is ≤ the cutoff (n · plane.p). Default
                # Maths::EPSILON = 1e-6.
                if dot(nv, point) <= d + T(1e-6)
                    tag = true
                    break
                end
            end
            inside[i] = tag
        end

        # Remap template vertex index → global pos index (only for
        # vertices used by at least one kept triangle).
        remap = zeros(Int, n_tmpl)
        for t in tmpl.faces
            i1, i2, i3 = Int(t[1]), Int(t[2]), Int(t[3])
            # Drop the triangle iff ALL THREE corners are inside.
            (inside[i1] && inside[i2] && inside[i3]) && continue
            for i in (i1, i2, i3)
                remap[i] != 0 && continue
                u_t = tmpl.position[i]
                outward = Vector3{T}(u_t...)
                point = c + R * outward
                push!(pos,   Point3{T}(point...))
                push!(norms, Vec3{T}(outward...))
                remap[i] = length(pos)
            end
            push!(tris, TriangleFace{Int}(remap[i1], remap[i2], remap[i3]))
        end
    end

    GeometryBasics.Mesh(pos, tris; normal=norms)
end

"""
    triangulate_sas(ac::AbstractAtomContainer{T}; probe_radius, density)

Convenience overload that builds the analytical SAS first.
"""
function triangulate_sas(ac::AbstractAtomContainer{T};
                         probe_radius::T = T(1.5),
                         density::T      = T(1.0)) where T
    triangulate_sas(compute_sas(ac; probe_radius); density)
end
