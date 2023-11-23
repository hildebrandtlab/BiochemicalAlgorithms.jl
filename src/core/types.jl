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
    M::_RowProjectionVector{T},
    v::T,
    i::Int
) where T = setindex!(M._base, v, getindex(M._rows, i))

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

const Position{T} = Vector3{<:Unitful.Length{T}}

@inline Position(r::Vector3{T}) where T = r * u"Å"
@inline Position(r::AbstractVector{T}) where T = Vector3(r) * u"Å"
@inline Position(rx::T, ry::T, rz::T) where T = Vector3(rx, ry, rz) * u"Å"

@inline Base.convert(::Type{Position{T}}, r::Vector3{T}) where T = Position(r)
@inline Base.zeros(::Type{Position{T}}) where T = Vector3(zeros(T, 3)u"Å")
