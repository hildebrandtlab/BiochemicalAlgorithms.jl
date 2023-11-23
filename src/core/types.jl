using StaticArrays

export
    Angstrom,
    Flags,
    Matrix3,
    MaybeInt,
    Position,
    Properties,
    Vector3,
    distance,
    squared_norm

const Vector3{T} = SVector{3, T}
const Matrix3{T} = SMatrix{3, 3, T}

const MaybeInt = Union{Nothing, Int}
const Properties = Dict{Symbol, Any}
const Flags = Set{Symbol}

const Angstrom{T <: Real} = Quantity{
    T,
    Unitful.𝐋,
    Unitful.FreeUnits{
        (Unitful.Unit{:Angstrom,Unitful.𝐋}(0, 1),),
        Unitful.𝐋,
        nothing
    }
}

@inline Angstrom(x::T) where {T <: Real} = x * u"Å"

#Base.show(io::IO, ::Type{Angstrom}) = print(io, "Angstrom")
#Base.show(io::IO, ::Type{Angstrom{T}}) where T = print(io, "Angstrom{$T}")

const Position{T} = Vector3{Angstrom{T}}

@inline Position(r::Vector3{T}) where T = r * u"Å"
@inline Position(r::AbstractVector{T}) where T = Vector3(r) * u"Å"
@inline Position(rx::T, ry::T, rz::T) where T = Vector3(rx, ry, rz) * u"Å"

@inline Base.convert(::Type{Position{T}}, r::Vector3{T}) where T = Position(r)
@inline Base.zeros(::Type{Position{T}}) where T = Vector3(zeros(T, 3)u"Å")

#Base.show(io::IO, ::Type{Position}) = print(io, "Position")
#Base.show(io::IO, ::Type{Position{T}}) where T = print(io, "Position{$T}")

squared_norm(v::Vector3{T}) where {T<:Real} = dot(v, v)
distance(v::Vector3{T}, w::Vector3{T}) where {T<:Real} = norm(v - w)
