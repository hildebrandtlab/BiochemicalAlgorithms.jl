const _bond_table_schema = Tables.Schema(
    (:idx, :a1, :a2, :order, :properties, :flags),
    (Int, Int, Int, BondOrderType, Properties, Flags)
)
const _bond_table_cols = _bond_table_schema.names
const _bond_table_cols_set = Set(_bond_table_cols)

@auto_hash_equals struct _BondTable <: Tables.AbstractColumns
    idx::Vector{Int}
    a1::Vector{Int}
    a2::Vector{Int}
    order::Vector{BondOrderType}
    properties::Vector{Properties}
    flags::Vector{Flags}

    _idx_map::Dict{Int,Int}

    function _BondTable()
        new(
            Int[],
            Int[],
            Int[],
            BondOrderType[],
            Properties[],
            Flags[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.istable(::Type{<: _BondTable}) = true
@inline Tables.columnaccess(::Type{<: _BondTable}) = true
@inline Tables.columns(bt::_BondTable) = bt

@inline function Tables.getcolumn(bt::_BondTable, nm::Symbol)
    @assert nm in _bond_table_cols "type _BondTable has no column $nm"
    getfield(bt, nm)
end
@inline Base.getproperty(bt::_BondTable, nm::Symbol) = getfield(bt, nm)

@inline Tables.getcolumn(bt::_BondTable, i::Int) = getfield(bt, Tables.columnnames(bt)[i])
@inline Tables.columnnames(::_BondTable) = _bond_table_cols
@inline Tables.schema(::_BondTable) = _bond_table_schema

@inline Base.size(bt::_BondTable) = (length(bt.idx), length(_bond_table_cols))
@inline Base.size(bt::_BondTable, dim) = size(bt)[dim]
@inline Base.length(bt::_BondTable) = size(bt, 1)

function Base.push!(
    bt::_BondTable,
    idx::Int,
    a1::Int,
    a2::Int,
    order::BondOrderType;
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    bt._idx_map[idx] = length(bt.idx) + 1
    push!(bt.idx, idx)
    push!(bt.a1, a1)
    push!(bt.a2, a2)
    push!(bt.order, order)
    push!(bt.properties, properties)
    push!(bt.flags, flags)
    bt
end

function _bond_table(itr)
    bt = _BondTable()
    for b in itr
        push!(bt, b.idx, b.a1, b.a2, b.order;
            properties = b.properties,
            flags = b.flags
        )
    end
    bt
end
Tables.materializer(::Type{_BondTable}) = _bond_table

@auto_hash_equals struct _BondTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_BondTable
end

@inline Tables.rowaccess(::Type{<: _BondTable}) = true
@inline Tables.rows(bt::_BondTable) = bt
@inline Tables.getcolumn(btr::_BondTableRow, nm::Symbol) = Tables.getcolumn(getfield(btr, :_tab), nm)[getfield(btr, :_row)]
@inline Tables.getcolumn(btr::_BondTableRow, i::Int) = getfield(btr, Tables.columnnames(btr)[i])
@inline Tables.columnnames(::_BondTableRow) = _bond_table_cols

@inline _row_by_idx(bt::_BondTable, idx::Int) = _BondTableRow(getfield(bt, :_idx_map)[idx], bt)

@inline function Base.getproperty(btr::_BondTableRow, nm::Symbol)
    nm === :idx        && return _getproperty(btr, :idx)::Int
    nm === :a1         && return _getproperty(btr, :a1)::Int
    nm === :a2         && return _getproperty(btr, :a2)::Int
    nm === :order      && return _getproperty(btr, :order)::BondOrderType
    nm === :properties && return _getproperty(btr, :properties)::Properties
    nm === :flags      && return _getproperty(btr, :flags)::Flags
    getindex(getfield(getfield(btr, :_tab), nm), getfield(btr, :_row))
end

@inline function _getproperty(btr::_BondTableRow, nm::Symbol)
    getindex(getproperty(getfield(btr, :_tab), nm), getfield(btr, :_row))
end

@inline Base.setproperty!(btr::_BondTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(btr, :_tab), nm), val, getfield(btr, :_row))

@inline Base.eltype(::_BondTable) = _BondTableRow
@inline Base.iterate(bt::_BondTable, st=1) = st > length(bt) ? nothing : (_BondTableRow(st, bt), st + 1)
