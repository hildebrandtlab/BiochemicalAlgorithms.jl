using StaticArrays

export Vector3, Matrix3, MaybeInt, Properties, Flags, squared_norm, distance

const Vector3{T} = SVector{3, T}
const Matrix3{T} = SMatrix{3, 3, T}

const MaybeInt = Union{Nothing, Int}
const Properties = Dict{Symbol, Any}
const Flags = Set{Symbol}

squared_norm(v::Vector3{T}) where {T<:Real} = dot(v, v)
distance(v::Vector3{T}, w::Vector3{T}) where {T<:Real} = norm(v - w)

struct RowProjectionVector{T} <: AbstractArray{T, 1}
    base::Vector{T}
    rows::Vector{Int}
end

@inline Base.size(M::RowProjectionVector) = (length(M.rows),)

@inline Base.getindex(
    M::RowProjectionVector,
    i::Int
) = getindex(getproperty(M, :base), getindex(getproperty(M, :rows), i))

@inline Base.setindex!(
    M::RowProjectionVector{T},
    v::T,
    i::Int
) where T = setindex!(getproperty(M, :base), v, getindex(getproperty(M, :rows), i))
