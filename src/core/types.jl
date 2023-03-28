using StaticArrays

export Vector3, Matrix3, MaybeInt, Properties, Flags, squared_norm

const Vector3{T} = SVector{3, T}
const Matrix3{T} = SMatrix{3, 3, T}

const MaybeInt = Union{Nothing, Int}
const Properties = Dict{Symbol, Any}
const Flags = Set{Symbol}

squared_norm(v::Vector3{T}) where {T<:Real} = dot(v, v)
distance(v::Vector3{T}, w::Vector3{T}) where {T<:Real} = norm(v - w)
