using StaticArrays

export Vector3, Matrix3, MaybeInt, Properties, squared_norm

const Vector3{T} = SVector{3, T}
const Matrix3{T} = SMatrix{3, 3, T}

const MaybeInt = Union{Missing, Int}
const Properties = Dict{String, Any}

squared_norm(v::Vector3{T}) where {T<:Real} = dot(v, v)
distance(v::Vector3{T}, w::Vector3{T}) where {T<:Real} = norm(v - w)