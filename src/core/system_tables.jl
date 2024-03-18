export
    SystemComponentTable

@auto_hash_equals struct SystemComponentTable{T, C <: AbstractSystemComponent{T}} <: AbstractSystemComponentTable{T}
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _table(ct::SystemComponentTable{T, C}) where {T, C} = _table(ct._sys, C)
@inline _hascolumn(::SystemComponentTable{T, C}, nm::Symbol) where {T, C} = _hascolumn(C, nm)

@inline function _element_by_idx(ct::SystemComponentTable{T, C}, idx::Int) where {T, C}
    C(ct._sys, _row_by_idx(_table(ct), idx))
end

@inline function Tables.getcolumn(ct::SystemComponentTable, nm::Symbol)
    col = Tables.getcolumn(_table(ct), nm)
    _RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _table(ct)._idx_map[idx], ct._idx)
    )
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

@inline function Base.getindex(ct::SystemComponentTable{T, C}, I) where {T, C}
    SystemComponentTable{T, C}(ct._sys, collect(Int, map(i -> ct._idx[i], I)))
end
