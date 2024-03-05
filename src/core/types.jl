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

squared_norm(v::Vector3{T}) where {T<:Real} = dot(v, v)
distance(v::Vector3{T}, w::Vector3{T}) where {T<:Real} = norm(v - w)

struct _RowProjectionVector{T} <: AbstractArray{T, 1}
    _base::Vector{T}
    _rows::Vector{Int}
end

@inline Base.size(M::_RowProjectionVector) = (length(M._rows),)

@inline Base.getindex(
    M::_RowProjectionVector,
    i::Int
) = getindex(getproperty(M, :_base), getindex(getproperty(M, :_rows), i))

@inline Base.setindex!(
    M::_RowProjectionVector{T},
    v::T,
    i::Int
) where T = setindex!(getproperty(M, :_base), v, getindex(getproperty(M, :_rows), i))
