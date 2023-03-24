using StaticArrays

export Vector3, Matrix3, MaybeInt, Properties

const Vector3{T} = SVector{3, T}
const Matrix3{T} = SMatrix{3, 3, T}

const MaybeInt = Union{Missing, Int}
const Properties = Dict{String, Any}
