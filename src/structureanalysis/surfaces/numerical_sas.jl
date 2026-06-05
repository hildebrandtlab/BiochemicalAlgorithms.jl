export
    NumericalSASResult,
    compute_numerical_sas

"""
    NumericalSASResult{T<:Real}

Result of a numerical SAS computation. Vector fields are indexed by atom
*position in the source `AbstractAtomContainer`* (i.e., `at.r[i]`-aligned),
so atoms with `radius == 0` keep zero contributions in `atom_areas` /
`atom_volumes`.

# Fields
 - `total_area::T`             — total solvent-accessible surface area (Å²)
 - `total_volume::T`           — total volume enclosed by the SAS (Å³)
 - `atom_areas::Vector{T}`     — per-atom area contribution (Å²)
 - `atom_volumes::Vector{T}`   — per-atom volume contribution (Å³)
 - `surface_vertices`          — non-occluded sample points (only if
                                 `compute_surface=true`)
 - `surface_normals`           — area-weighted outward normals at each point
                                 (sum of `‖normal‖` equals `total_area`)
 - `atom_surface_vertices`     — per-atom sample points (only if
                                 `compute_surface_per_atom=true`)
 - `atom_surface_normals`      — per-atom area-weighted normals
"""
Base.@kwdef mutable struct NumericalSASResult{T<:Real}
    total_area::T = zero(T)
    total_volume::T = zero(T)
    atom_areas::Vector{T} = T[]
    atom_volumes::Vector{T} = T[]
    surface_vertices::Vector{Point3{T}} = Point3{T}[]
    surface_normals::Vector{Vec3{T}}    = Vec3{T}[]
    atom_surface_vertices::Vector{Vector{Point3{T}}} = Vector{Point3{T}}[]
    atom_surface_normals::Vector{Vector{Vec3{T}}}    = Vector{Vec3{T}}[]
end

@inline function _tessellation_for(::Type{T}, num_points::Int) where {T<:Real}
    subdivisions = 0
    while 10 * 4^subdivisions + 2 < num_points
        subdivisions += 1
    end
    icosphere(T, subdivisions)
end

"""
    $(TYPEDSIGNATURES)

Compute the solvent-accessible surface (SAS) of `ac` numerically via the
Eisenhaber/Argos double-cubic-lattice method (Eisenhaber, Lijnzaad, Argos,
Sander & Scharf, *J. Comput. Chem.* 1995, **15**, 273-284; Eisenhaber & Argos,
*J. Comput. Chem.* 1993, **14**, 1272-1280).

Every atom whose `radius > 0.001` Å is inflated by `probe_radius` and sampled
with `number_of_points` points distributed via icosahedral subdivision (the
returned tessellation has the smallest `10·4ⁿ + 2 ≥ number_of_points`
vertices). For each sample, occlusion against neighboring inflated atoms is
tested; non-occluded samples contribute to area, volume, and (optionally)
surface point clouds.

# Keyword arguments
 - `probe_radius::T = 1.5`           — radius of the spherical probe (Å)
 - `number_of_points::Int = 400`     — lower bound on sample points per atom
 - `compute_area::Bool = true`       — populate `total_area` / `atom_areas`
 - `compute_volume::Bool = true`     — populate `total_volume` / `atom_volumes`
 - `compute_surface::Bool = false`   — populate aggregate `surface_*` fields
 - `compute_surface_per_atom::Bool = false` — populate `atom_surface_*` fields

The volume is integrated via the divergence theorem
`V = ⅓ ∮ r·n dA`, using the geometric center of the radius-bearing atoms as
the origin (per BALL's convention).
"""
function compute_numerical_sas(
    ac::AbstractAtomContainer{T};
    probe_radius::T = T(1.5),
    number_of_points::Int = 400,
    compute_area::Bool = true,
    compute_volume::Bool = true,
    compute_surface::Bool = false,
    compute_surface_per_atom::Bool = false,
) where {T<:Real}
    result = NumericalSASResult{T}()
    at = atoms(ac)
    n = length(at.r)
    n == 0 && return result

    valid_idx = findall(>(T(0.001)), at.radius)
    if compute_area || compute_volume
        result.atom_areas   = zeros(T, n)
        result.atom_volumes = zeros(T, n)
    end
    if compute_surface_per_atom
        result.atom_surface_vertices = [Point3{T}[] for _ in 1:n]
        result.atom_surface_normals  = [Vec3{T}[]   for _ in 1:n]
    end
    isempty(valid_idx) && return result

    template = _tessellation_for(T, number_of_points)
    template_points = template.position
    num_points = length(template_points)
    unit_area_per_point = 4 * T(π) / num_points
    unit_volume = unit_area_per_point / 3

    positions = [Point3{T}(at.r[i]...) for i in valid_idx]
    radii     = T[at.radius[i] for i in valid_idx]
    m = length(valid_idx)

    r_max  = maximum(radii) + probe_radius
    cutoff = 2 * r_max

    neighbors = [Int[] for _ in 1:m]
    if m > 1
        pairs = neighborlist(; positions=positions, cutoff=cutoff)
        @inbounds for (i, j, _) in pairs
            ri = radii[i] + probe_radius
            rj = radii[j] + probe_radius
            d  = positions[i] - positions[j]
            if dot(d, d) <= (ri + rj)^2
                push!(neighbors[i], j)
                push!(neighbors[j], i)
            end
        end
    end

    cog = sum(positions) / m

    @inbounds for k in 1:m
        atom_idx = valid_idx[k]
        center   = positions[k]
        inflated = radii[k] + probe_radius

        nbrs = neighbors[k]
        num_occluded = 0
        dr = Vec3{T}(zero(T), zero(T), zero(T))

        for p in 1:num_points
            unit = template_points[p]
            world_point = center + inflated * unit

            is_occluded = false
            for j in nbrs
                pr = radii[j] + probe_radius
                d  = world_point - positions[j]
                if dot(d, d) <= pr * pr
                    is_occluded = true
                    break
                end
            end

            if is_occluded
                num_occluded += 1
            else
                if compute_volume
                    dr = Vec3{T}(dr[1] + unit[1], dr[2] + unit[2], dr[3] + unit[3])
                end
                if compute_surface
                    weight = inflated * inflated * unit_area_per_point
                    push!(result.surface_vertices, world_point)
                    push!(result.surface_normals,
                          Vec3{T}(unit[1] * weight, unit[2] * weight, unit[3] * weight))
                end
                if compute_surface_per_atom
                    push!(result.atom_surface_vertices[atom_idx], world_point)
                    push!(result.atom_surface_normals[atom_idx], Vec3{T}(unit[1], unit[2], unit[3]))
                end
            end
        end

        n_exposed = num_points - num_occluded

        if compute_area
            atom_area = inflated * inflated * unit_area_per_point * n_exposed
            result.total_area += atom_area
            result.atom_areas[atom_idx] = atom_area
        end
        if compute_volume
            offset = center - cog
            atom_volume = inflated * inflated * unit_volume *
                          (offset[1] * dr[1] + offset[2] * dr[2] + offset[3] * dr[3] +
                           inflated * n_exposed)
            result.total_volume += atom_volume
            result.atom_volumes[atom_idx] = atom_volume
        end
        if compute_surface_per_atom
            weight = inflated * inflated * unit_area_per_point
            normals = result.atom_surface_normals[atom_idx]
            for q in eachindex(normals)
                normals[q] = normals[q] * weight
            end
        end
    end

    result
end
