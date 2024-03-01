@auto_hash_equals struct _Chain
    idx::Int
    name::String
    properties::Properties
    flags::Flags

    function _Chain(;
        idx::Int = 0,
        name::String = "",
        properties::Properties = Properties(),
        flags::Flags = Flags()
    )
        new(idx, name, properties, flags)
    end
end

const _chain_table_cols = fieldnames(_Chain)
const _chain_table_cols_set = Set(_chain_table_cols)
const _chain_table_cols_priv = Set([:molecule_idx])

@auto_hash_equals struct _ChainTable <: Tables.AbstractColumns
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

@inline Tables.istable(::Type{<: _ChainTable}) = true
@inline Tables.columnaccess(::Type{<: _ChainTable}) = true
@inline Tables.columns(ct::_ChainTable) = ct

@inline function Tables.getcolumn(ct::_ChainTable, nm::Symbol)
    @assert nm in _chain_table_cols_priv || nm in _chain_table_cols_set "type _ChainTable has no column $nm"
    getfield(ct, nm)
end
@inline Base.getproperty(ct::_ChainTable, nm::Symbol) = getfield(ct, nm)

@inline Tables.getcolumn(ct::_ChainTable, i::Int) = getfield(ct, Tables.columnnames(ct)[i])
@inline Tables.columnnames(::_ChainTable) = _chain_table_cols
@inline Tables.schema(::_ChainTable) = Tables.Schema(fieldnames(_Chain), fieldtypes(_Chain))

@inline Base.size(ct::_ChainTable) = (length(ct.idx), length(_chain_table_cols))
@inline Base.size(ct::_ChainTable, dim) = size(ct)[dim]
@inline Base.length(ct::_ChainTable) = size(ct, 1)

function Base.push!(ct::_ChainTable, t::_Chain, molecule_idx::Int)
    getfield(ct, :_idx_map)[t.idx] = length(ct.idx) + 1
    for fn in _chain_table_cols
        push!(getfield(ct, Symbol(fn)), getfield(t, Symbol(fn)))
    end
    push!(getfield(ct, :molecule_idx), molecule_idx)
    ct
end

function _chain_table(itr)
    ct = _ChainTable()
    for c in itr
        push!(ct, _Chain(;
                idx = c.idx,
                name = c.name,
                properties = c.properties,
                flags = c.flags
            ),
            Tables.getcolumn(c, :molecule_idx)
       )
    end
    ct
end
@inline Tables.materializer(::Type{_ChainTable}) = itr -> _chain_table(itr)

@auto_hash_equals struct _ChainTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_ChainTable
end

@inline Tables.rowaccess(::Type{<: _ChainTable}) = true
@inline Tables.rows(ct::_ChainTable) = ct
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
