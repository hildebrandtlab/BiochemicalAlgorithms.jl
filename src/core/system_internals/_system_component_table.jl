abstract type _AbstractSystemComponentTable <: AbstractColumnTable end

const _IdxMap = Dict{Int, Int}

@inline function _rowno_by_idx(at::_AbstractSystemComponentTable, idx::Int)
    getindex(getfield(at, :_idx_map), idx)
end

@inline function _row_by_idx(at::_AbstractSystemComponentTable, idx::Int)
    ColumnTableRow(_rowno_by_idx(at, idx), at)
end

@inline function _rebuild_idx_map!(at::_AbstractSystemComponentTable)
    empty!(at._idx_map)
    merge!(at._idx_map, Dict(v => k for (k, v) in enumerate(at.idx)))
    at
end

function Base.delete!(at::_AbstractSystemComponentTable, idx::Int)
    _delete!(at, _rowno_by_idx(at, idx))
    _rebuild_idx_map!(at)
end

function Base.delete!(at::_AbstractSystemComponentTable, idx::Vector{Int})
    rownos = _rowno_by_idx.(Ref(at), idx)
    unique!(rownos)
    sort!(rownos; rev = true)
    for rowno in rownos
        _delete!(at, rowno)
    end
    _rebuild_idx_map!(at)
end

@inline Base.size(at::_AbstractSystemComponentTable) = (length(at.idx), length(Tables.columnnames(at)))

@inline function Tables.getcolumn(at::_AbstractSystemComponentTable, nm::Symbol)
    @assert _hascolumn(at, nm) "type $(typeof(at)) has no column $nm"
    getfield(at, nm)
end
