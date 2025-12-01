export
    BoundingBox,
    compute_bounding_box

struct BoundingBox{T<:Real}
    min::Vector3{T}
    max::Vector3{T}
end

"""
    $(TYPEDSIGNATURES)

Computes the axis-aligned bounding box of the given `AbstractAtomContainer`.
"""
function compute_bounding_box(ac::AbstractAtomContainer{T}) where T
    at = atoms(ac)

    min_coords = Vector3{T}(
        minimum(getindex.(at.r, 1); init = zero(T)),
        minimum(getindex.(at.r, 2); init = zero(T)),
        minimum(getindex.(at.r, 3); init = zero(T))
   )

    max_coords = Vector3{T}(
        maximum(getindex.(at.r, 1); init = zero(T)),
        maximum(getindex.(at.r, 2); init = zero(T)),
        maximum(getindex.(at.r, 3); init = zero(T))
    )

    BoundingBox{T}(min_coords, max_coords)
end
