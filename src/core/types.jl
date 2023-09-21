using StaticArrays

export Vector3, Matrix3, MaybeInt, Properties, Flags

const Vector3{T} = SVector{3, T}
const Matrix3{T} = SMatrix{3, 3, T}

const MaybeInt = Union{Nothing, Int}
const Properties = Dict{Symbol, Any}
const Flags = Set{Symbol}
