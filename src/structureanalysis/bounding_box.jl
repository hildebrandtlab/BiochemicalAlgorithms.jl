export BoundingBox, compute_bounding_box

struct BoundingBox{T<:Real}
    min::Vector{T}
    max::Vector{T}
end


"""
    $(TYPEDSIGNATURES)

    Computes the axis-aligned bounding box of the given `AbstractAtomContainer`.
"""

function compute_bounding_box(ac::AbstractAtomContainer{T}) where {T<:Real}
    min_coords = Vector3{T}(
        minimum(getindex.(atoms(ac).r, 1)),
        minimum(getindex.(atoms(ac).r, 2)),
        minimum(getindex.(atoms(ac).r, 3))
   )

    max_coords = Vector3{T}(
        maximum(getindex.(atoms(ac).r, 1)),
        maximum(getindex.(atoms(ac).r, 2)),
        maximum(getindex.(atoms(ac).r, 3))
    )

    BoundingBox{T}(min_coords, max_coords)
end
   