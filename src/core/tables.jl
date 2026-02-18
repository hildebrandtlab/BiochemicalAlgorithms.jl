export
    AbstractColumnTable,
    ColumnTableRow

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

@inline function Base.show(io::IO, ::MIME"text/html", at::AbstractColumnTable; kwargs...)
    _show(io, at; backend = :html, vertical_crop_mode = :middle, kwargs...)
end

@inline function Base.show(io::IO, ::MIME"text/markdown", at::AbstractColumnTable; kwargs...)
    _show(io, at; backend = :markdown, title = "", kwargs...)
end

@inline function Base.show(io::IO, ::MIME"text/plain", at::AbstractColumnTable; kwargs...)
    _show(io, at; reserved_display_lines = 2, vertical_crop_mode = :middle, kwargs...)
end

@inline function _show(io::IO, at::AbstractColumnTable; kwargs...)
    PrettyTables.pretty_table(io, at;
        alignment = :l,
        column_labels = collect(Tables.columnnames(at)),
        row_number_column_label = "#",
        show_row_number_column = true,
        title  = "$(typeof(at)) with $(length(at)) rows:",
        title_alignment = :l,
        kwargs...
    )
end

@inline function _filter_idx(f::Function, at::AbstractColumnTable)
    collect(Int, _filter_select(f, at, :idx))
end

@inline function _filter_select(f::Function, at::AbstractColumnTable, col::Symbol)
    (getproperty(a, col) for a in TableOperations.filter(f, at))
end

struct ColumnTableRow{CT <: AbstractColumnTable} <: Tables.AbstractRow
    _row::Int
    _tab::CT
end

@inline Tables.getcolumn(ctr::ColumnTableRow, nm::Symbol) = Tables.getcolumn(getfield(ctr, :_tab), nm)[getfield(ctr, :_row)]
@inline Tables.getcolumn(ctr::ColumnTableRow, i::Int) = getproperty(ctr, Tables.columnnames(ctr)[i])

@inline Tables.columnnames(ctr::ColumnTableRow) = Tables.columnnames(getfield(ctr, :_tab))

@inline function Base.getproperty(atr::ColumnTableRow, nm::Symbol)
    getindex(getfield(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

@inline function Base.setproperty!(atr::ColumnTableRow, nm::Symbol, val)
    setindex!(getproperty(getfield(atr, :_tab), nm), val, getfield(atr, :_row))
end

@inline Base.eltype(::CT) where {CT <: AbstractColumnTable} = ColumnTableRow{CT}
@inline Base.iterate(at::AbstractColumnTable, st=1) = st > length(at) ? nothing : (ColumnTableRow(st, at), st + 1)
