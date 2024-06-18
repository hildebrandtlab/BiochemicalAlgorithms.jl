export
    AbstractColumnTable

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
    _show(io, at; backend = Val(:html))
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

abstract type _AbstractColumnTableRow <: Tables.AbstractRow end
