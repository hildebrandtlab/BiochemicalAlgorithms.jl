const _chain_table_schema = Tables.Schema(
    (:idx, :name, :properties, :flags),
    (Int, String, Properties, Flags)
)
const _chain_table_cols = _chain_table_schema.names
const _chain_table_cols_set = Set(_chain_table_cols)
const _chain_table_cols_priv = Set([:molecule_idx])

@auto_hash_equals struct _ChainTable <: AbstractColumnTable
    # public columns
    idx::Vector{Int}
    name::Vector{String}
    properties::Vector{Properties}
    flags::Vector{Flags}

    # private columns
    molecule_idx::Vector{Int}

    # internals
    _idx_map::Dict{Int,Int}

    function _ChainTable()
        new(
            Int[],
            String[],
            Properties[],
            Flags[],
            Int[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.columnnames(::_ChainTable) = _chain_table_cols
@inline Tables.schema(::_ChainTable) = _chain_table_schema

@inline function Tables.getcolumn(ct::_ChainTable, nm::Symbol)
    @assert nm in _chain_table_cols_priv || nm in _chain_table_cols_set "type _ChainTable has no column $nm"
    getfield(ct, nm)
end

@inline Base.size(ct::_ChainTable) = (length(ct.idx), length(_chain_table_cols))

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

@auto_hash_equals struct _ChainTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_ChainTable
end

@inline Tables.getcolumn(ctr::_ChainTableRow, nm::Symbol) = Tables.getcolumn(getfield(ctr, :_tab), nm)[getfield(ctr, :_row)]
@inline Tables.getcolumn(ctr::_ChainTableRow, i::Int) = getfield(ctr, Tables.columnnames(ctr)[i])
@inline Tables.columnnames(::_ChainTableRow) = _chain_table_cols

@inline _row_by_idx(ct::_ChainTable, idx::Int) = _ChainTableRow(getfield(ct, :_idx_map)[idx], ct)

@inline function Base.getproperty(ctr::_ChainTableRow, nm::Symbol)
    nm === :idx         && return _getproperty(ctr, :idx)::Int
    nm === :name        && return _getproperty(ctr, :name)::String
    nm === :properties  && return _getproperty(ctr, :properties)::Properties
    nm === :flags       && return _getproperty(ctr, :flags)::Flags
    nm === :molecule_idx && return _getproperty(ctr, :molecule_idx)::Int
    getindex(getfield(getfield(ctr, :_tab), nm), getfield(ctr, :_row))
end

@inline function _getproperty(ctr::_ChainTableRow, nm::Symbol)
    getindex(getproperty(getfield(ctr, :_tab), nm), getfield(ctr, :_row))
end

@inline Base.setproperty!(ctr::_ChainTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(ctr, :_tab), nm), val, getfield(ctr, :_row))

@inline Base.eltype(::_ChainTable) = _ChainTableRow
@inline Base.iterate(ct::_ChainTable, st=1) = st > length(ct) ? nothing : (_ChainTableRow(st, ct), st + 1)
