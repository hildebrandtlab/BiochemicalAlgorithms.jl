export
    SolventAccessibleSurface,
    compute_sas

# The SAS is structurally the dual of the reduced surface:
#   * an SASVertex sits at the center of each RS probe sphere (= each RSFace)
#   * an SASEdge corresponds to each RSEdge (the torus's contact circle)
#   * an SASFace covers the part of an inflated atom not occluded by other
#     atoms; it is the dual of an RSVertex.
# All references are stored by index into the parent ReducedSurface arrays
# (they index `rs.faces`, `rs.edges`, and `rs.vertices` respectively).

"""
    SASVertex{T}

A point on the SAS — specifically, the center of a probe sphere resting on
an atom triple. Indexed by the parent `RSFace`.
"""
struct SASVertex{T<:Real}
    point::Vector3{T}
    rs_face::Int            # source RSFace index in the parent RS
    edges::Vector{Int}      # incident SASEdge indices
    faces::Vector{Int}      # incident SASFace indices (the three atoms)
end

"""
    SASEdge{T}

An arc on the SAS — the path traced by the probe center while rolling between
two adjacent contact triples. Indexed by the parent `RSEdge`.
"""
struct SASEdge{T<:Real}
    rs_edge::Int            # source RSEdge index
    v1::Int                 # SASVertex indices (= source RSEdge's f1/f2)
    v2::Int
    f1::Int                 # SASFace indices (= source RSEdge's two atom vertices)
    f2::Int
    angle::T                # arc angle (radians)
    circle::Circle3{T}      # carrier circle in 3D (center + axis + radius)
end

"""
    SASFace{T}

A connected component of the SAS lying on a single inflated atom. The sphere
`sphere` is the atom sphere expanded by the probe radius; `edges` are the
SAS edges that bound this region. `orientation` records whether each edge is
traversed in its native direction (`true`) or reversed (`false`) when walking
the boundary in face-normal order.
"""
struct SASFace{T<:Real}
    rs_vertex::Int          # source RSVertex index
    sphere::Sphere{T}       # atom + probe inflation
    edges::Vector{Int}      # SASEdge indices on this face's boundary
    orientation::Vector{Bool}
    vertices::Vector{Int}   # SASVertex indices around the boundary
end

"""
    SolventAccessibleSurface{T}

Analytical solvent-accessible surface (SAS) built as the dual of a
[`ReducedSurface`](@ref). Each atom contributes one or more spherical-cap
faces; arcs and points on the surface come from the rolling-probe topology.
"""
mutable struct SolventAccessibleSurface{T<:Real} <: AbstractMolecularSurface{T}
    reduced_surface::ReducedSurface{T}
    vertices::Vector{SASVertex{T}}
    edges::Vector{SASEdge{T}}
    faces::Vector{SASFace{T}}
end

@inline nvertices(sas::SolventAccessibleSurface) =
    count(v -> v.rs_face != 0, sas.vertices)
@inline nedges(sas::SolventAccessibleSurface)    =
    count(e -> e.rs_edge != 0, sas.edges)
@inline nfaces(sas::SolventAccessibleSurface)    =
    count(f -> f.rs_vertex != 0, sas.faces)

"""
    $(TYPEDSIGNATURES)

Compute the analytical solvent-accessible surface (SAS) from a precomputed
[`ReducedSurface`](@ref). The SAS shares its `T` parameter with the source RS
and reuses its topology rather than recomputing it.
"""
function compute_sas(rs::ReducedSurface{T}) where {T<:Real}
    # Sentinel entries (atom == 0 / v1 == 0) may appear ANYWHERE in the RS
    # arrays after `_rs_correct!` / `_rs_face_remove!`; iterate raw indices
    # and skip them so SAS arrays stay index-aligned with the RS arrays.
    sas_vertices = Vector{SASVertex{T}}(undef, length(rs.faces))
    for (j, rsface) in enumerate(rs.faces)
        if rsface.v1 == 0
            sas_vertices[j] = SASVertex{T}(Vector3{T}(0, 0, 0), 0, Int[], Int[])  # sentinel: rs_face=0
            continue
        end
        sas_vertices[j] = SASVertex{T}(
            rsface.center,
            j,
            copy(rsface.edges),
            [rsface.v1, rsface.v2, rsface.v3],
        )
    end

    sas_edges = Vector{SASEdge{T}}(undef, length(rs.edges))
    for j in 1:length(rs.edges)
        rsedge = rs.edges[j]
        if rsedge.v1 == 0
            sas_edges[j] = SASEdge{T}(0, 0, 0, 0, 0, zero(T),
                                       Circle3{T}(Vector3{T}(0,0,0), Vector3{T}(0,0,1), zero(T)))
            continue
        end
        v1 = rsedge.f1
        v2 = rsedge.f2
        f1 = rsedge.v1
        f2 = rsedge.v2
        atom1 = rs.atoms[rs.vertices[f1].atom]
        atom2 = rs.atoms[rs.vertices[f2].atom]
        ax = Vector3{T}((atom1.center - atom2.center)...)
        ℓ = norm(ax)
        nvec = iszero(ℓ) ? rsedge.contact_circle1.n : ax / ℓ
        circle = Circle3{T}(rsedge.center_of_torus, nvec, rsedge.radius_of_torus)
        sas_edges[j] = SASEdge{T}(j, v1, v2, f1, f2, rsedge.angle, circle)
    end

    sas_faces = Vector{SASFace{T}}(undef, length(rs.vertices))
    for j in 1:length(rs.vertices)
        rsvertex = rs.vertices[j]
        if rsvertex.atom == 0
            sas_faces[j] = SASFace{T}(0,
                Sphere{T}(Vector3{T}(0, 0, 0), zero(T)), Int[], Bool[], Int[])
            continue
        end
        atom = rs.atoms[rsvertex.atom]
        sphere = Sphere{T}(atom.center, atom.r + rs.probe_radius)
        # The face's edges are all RS edges incident to this vertex.
        edges = sort!(collect(rsvertex.edges))
        # Orientation: traversed forwards if this RS vertex is the "v1" of the
        # source RS edge, backwards if it's "v2".
        orient = Bool[
            (rs.edges[ei].v1 == j) for ei in edges
        ]
        # Vertices around the boundary = the RS faces incident to this vertex
        # (each face has its dual SAS vertex).
        verts = sort!(collect(rsvertex.faces))
        sas_faces[j] = SASFace{T}(j, sphere, edges, orient, verts)
    end

    SolventAccessibleSurface{T}(rs, sas_vertices, sas_edges, sas_faces)
end

"""
    $(TYPEDSIGNATURES)

Compute the analytical solvent-accessible surface (SAS) directly from an
atom container. Equivalent to
`compute_sas(compute_reduced_surface(ac; probe_radius))`.
"""
function compute_sas(ac::AbstractAtomContainer{T}; probe_radius::T = T(1.5)) where T
    compute_sas(compute_reduced_surface(ac; probe_radius))
end
