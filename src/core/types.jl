export
    Flags,
    Matrix3,
    MaybeInt,
    Properties,
    Vector3,
    distance,
    squared_norm

const Vector3{T} = StaticArrays.SVector{3, T}
const Matrix3{T} = StaticArrays.SMatrix{3, 3, T}

const MaybeInt = Union{Nothing, Int}
const Properties = Dict{Symbol, Any}
const Flags = Set{Symbol}

@inline squared_norm(v::AbstractVector) = dot(v, v)
@inline distance(v::AbstractVector{T}, w::AbstractVector{T}) where T = norm(v - w)
