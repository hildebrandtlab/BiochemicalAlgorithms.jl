export
    AbstractSystemComponentTable,
    SystemComponentTable,
    SystemComponentTableCol,
    revalidate_indices!

"""
    $(TYPEDEF)

Abstract base type for all Tables.jl-compatible system component tables.
"""
abstract type AbstractSystemComponentTable{T <: Real} <: AbstractColumnTable end

@inline function Base.getproperty(at::T, nm::Symbol) where {T <: AbstractSystemComponentTable}
    hasfield(T, nm) && return getfield(at, nm)
    Tables.getcolumn(at, nm)
end

@inline function Base.setproperty!(at::T, nm::Symbol, val) where {T <: AbstractSystemComponentTable}
    _hascolumn(at, nm) && error("$T columns cannot be set directly! Did you mean to use broadcast assignment (.=)?")
    hasfield(T, nm) || error("type $T has no field $nm")
    setfield!(at, nm, val)
end

@inline _table(at::AbstractSystemComponentTable) = at
@inline _hascolumn(::AbstractSystemComponentTable, ::Symbol) = false

"""
    struct SystemComponentTable{T, C <: AbstractSystemComponent{T}} <: AbstractSystemComponentTable{T}

Tables.jl-compatible system component table for a specific `System{T}` and system component `C` (e.g.,
`Atom{T}`, `Bond{T}`, etc.).
"""
@auto_hash_equals struct SystemComponentTable{T, C <: AbstractSystemComponent{T}} <: AbstractSystemComponentTable{T}
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _table(ct::SystemComponentTable{T, C}) where {T, C} = _table(ct._sys, C)
@inline _hascolumn(::SystemComponentTable{T, C}, nm::Symbol) where {T, C} = _hascolumn(C, nm)

@inline function _element_by_idx(ct::SystemComponentTable{T, C}, idx::Int) where {T, C}
    C(ct._sys, idx)
end

@inline Tables.columnnames(ct::SystemComponentTable) = Tables.columnnames(_table(ct))
@inline Tables.schema(ct::SystemComponentTable) = Tables.schema(_table(ct))

@inline function Base.filter(f::Function, ct::SystemComponentTable{T, C}) where {T, C}
    SystemComponentTable{T, C}(ct._sys, _filter_idx(f, ct))
end

@inline function Base.iterate(ct::SystemComponentTable, st = 1)
    st > length(ct) ?
        nothing :
        (_element_by_idx(ct, ct._idx[st]), st + 1)
end

@inline Base.eltype(::SystemComponentTable{T, C}) where {T, C} = C
@inline Base.size(ct::SystemComponentTable) = (length(ct._idx), length(Tables.columnnames(ct)))
@inline Base.getindex(ct::SystemComponentTable, i::Int) = _element_by_idx(ct, ct._idx[i])
@inline Base.getindex(ct::SystemComponentTable, ::Colon) = ct

@inline function Base.getindex(ct::SystemComponentTable{T, C}, I::AbstractArray) where {T, C}
    SystemComponentTable{T, C}(ct._sys, collect(Int, map(i -> ct._idx[i], I)))
end

@inline _row_by_idx(ct::SystemComponentTable, idx::Int) = _row_by_idx(_table(ct), idx)

"""
    struct SystemComponentTableCol{T} <: AbstractArray{T, 1}

`Vector`-like representation of a single `SystemComponentTable` column.
"""
struct SystemComponentTableCol{T} <: AbstractArray{T, 1}
    _base::Vector{T}
    _idx::Vector{Int}
    _idx_map::Dict{Int, Int}
end

@inline function Base.eltype(M::SystemComponentTableCol)
    eltype(M._base)
end

@inline function Base.size(M::SystemComponentTableCol)
    (length(M._idx),)
end

@inline function Base.getindex(M::SystemComponentTableCol, i::Int)
    getindex(M._base, getindex(M._idx_map, getindex(M._idx, i)))
end

@inline function Base.setindex!(M::SystemComponentTableCol{T}, v::T, i::Int) where T
    setindex!(M._base, v, getindex(M._idx_map, getindex(M._idx, i)))
end

@inline function Tables.getcolumn(ct::SystemComponentTable, nm::Symbol)
    col = Tables.getcolumn(_table(ct), nm)
    SystemComponentTableCol{eltype(col)}(
        col,
        ct._idx,
        _table(ct)._idx_map
    )
end

"""
    revalidate_indices!(::SystemComponentTable)
    revalidate_indices!(::SystemComponentTableCol)

Removes remnants of previously removed system components from the given table or table column.
"""
function revalidate_indices!(ct::SystemComponentTable)
    common_idx = ct._idx ∩ keys(_table(ct)._idx_map)
    empty!(ct._idx)
    append!(ct._idx, common_idx)
    ct
end

function revalidate_indices!(rv::SystemComponentTableCol)
    common_idx = rv._idx ∩ keys(rv._idx_map)
    empty!(rv._idx)
    append!(rv._idx, common_idx)
    rv
end
