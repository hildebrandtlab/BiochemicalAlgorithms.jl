const _chain_table_cols_main = (:idx, :name)
const _chain_table_cols_extra = (:properties, :flags, :molecule_idx)
const _chain_table_cols = (_chain_table_cols_main..., _chain_table_cols_extra...)
const _chain_table_cols_set = Set(_chain_table_cols)
const _chain_table_schema = Tables.Schema(
    _chain_table_cols_main,
    (Int, String)
)

@auto_hash_equals struct _ChainTable <: _AbstractSystemComponentTable
    # public columns
    idx::Vector{Int}
    name::Vector{String}

    # private columns
    properties::Vector{Properties}
    flags::Vector{Flags}
    molecule_idx::Vector{Int}

    # internals
    _idx_map::_IdxMap

    function _ChainTable()
        new(
            Int[],
            String[],
            Properties[],
            Flags[],
            Int[],
            _IdxMap()
        )
    end
end

@inline function _hascolumn(::Union{_ChainTable, Type{_ChainTable}}, nm::Symbol)
    nm in _chain_table_cols_set
end

@inline function Tables.columnnames(::_ChainTable)
    _chain_table_cols_main
end

@inline function Tables.schema(::_ChainTable)
    _chain_table_schema
end

function Base.push!(
    ct::_ChainTable,
    idx::Int,
    molecule_idx::Int;
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    ct._idx_map[idx] = length(ct.idx) + 1
    push!(ct.idx, idx)
    push!(ct.name, name)
    push!(ct.properties, properties)
    push!(ct.flags, flags)
    push!(ct.molecule_idx, molecule_idx)
    ct
end

function _delete!(ct::_ChainTable, rowno::Int)
    deleteat!(ct.idx, rowno)
    deleteat!(ct.name, rowno)
    deleteat!(ct.properties, rowno)
    deleteat!(ct.flags, rowno)
    deleteat!(ct.molecule_idx, rowno)
    nothing
end

function Base.empty!(ct::_ChainTable)
    empty!(ct.idx)
    empty!(ct.name)
    empty!(ct.properties)
    empty!(ct.flags)
    empty!(ct.molecule_idx)
    empty!(ct._idx_map)
    ct
end

function _chain_table(itr)
    ct = _ChainTable()
    for c in itr
        push!(ct, c.idx, c.molecule_idx;
            name = c.name,
            properties = c.properties,
            flags = c.flags
       )
    end
    ct
end
@inline Tables.materializer(::Type{_ChainTable}) = itr -> _chain_table(itr)
