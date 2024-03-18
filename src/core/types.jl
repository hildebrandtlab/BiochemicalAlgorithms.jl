export
    AbstractColumnTable,
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
    M::_RowProjectionVector{T},
    v::T,
    i::Int
) where T = setindex!(M._base, v, getindex(M._rows, i))

"""
    $(TYPEDEF)

Abstract base type for all Tables.jl-compatible column tables.
"""
abstract type AbstractColumnTable <: Tables.AbstractColumns end

@inline Tables.istable(::Type{<: AbstractColumnTable}) = true
@inline Tables.columnaccess(::Type{<: AbstractColumnTable}) = true
@inline Tables.columns(at::AbstractColumnTable) = at
@inline Tables.rows(at::AbstractColumnTable) = at
@inline Tables.getcolumn(at::AbstractColumnTable, i::Int) = Tables.getcolumn(at, Tables.columnnames(at)[i])

@inline Base.getproperty(at::AbstractColumnTable, nm::Symbol) = getfield(at, nm)
@inline Base.size(at::AbstractColumnTable, dim) = size(at)[dim]
@inline Base.length(at::AbstractColumnTable) = size(at, 1)
@inline Base.keys(at::AbstractColumnTable) = LinearIndices((length(at),))

@inline function Base.show(io::IO, at::AbstractColumnTable) 
    print(io, "$(typeof(at)) with $(length(at)) rows")
end

@inline function Base.show(io::IO, ::MIME"text/html", at::AbstractColumnTable)
    show(io, at; backend = Val(:html))
end

@inline function Base.show(io::IO, ::MIME"text/plain", at::AbstractColumnTable)
    _show(io, at; reserved_display_lines = 2)
end

@inline function _show(io::IO, at::AbstractColumnTable; kwargs...)
    PrettyTables.pretty_table(io, at;
        alignment = :l,
        header = collect(Tables.columnnames(at)),
        max_num_of_rows = 50,
        row_number_column_title = "#",
        show_row_number = true,
        title  = "$(typeof(at)) with $(length(at)) rows:",
        vcrop_mode = :middle,
        kwargs...
    )
end

@inline function _filter_idx(f::Function, at::AbstractColumnTable)
    collect(Int, _filter_select(f, at, :idx))
end

@inline function _filter_select(f::Function, at::AbstractColumnTable, col::Symbol)
    (getproperty(a, col) for a in TableOperations.filter(f, at))
end
