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

@inline squared_norm(v::Vector3{T}) where T = dot(v, v)
@inline distance(v::Vector3{T}, w::Vector3{T}) where T = norm(v - w)

struct _RowProjectionVector{T} <: AbstractArray{T, 1}
    _base::Vector{T}
    _rows::Vector{Int}
end

@inline Base.size(M::_RowProjectionVector) = (length(M._rows),)

@inline Base.getindex(
    M::_RowProjectionVector,
    i::Int
) = getindex(M._base, getindex(M._rows, i))

@inline Base.setindex!(
    M::_RowProjectionVector,
    v,
    i::Int
) = setindex!(M._base, v, getindex(M._rows, i))
