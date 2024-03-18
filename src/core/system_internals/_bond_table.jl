const _bond_table_schema = Tables.Schema(
    (:idx, :a1, :a2, :order),
    (Int, Int, Int, BondOrderType)
)
const _bond_table_cols = _bond_table_schema.names
const _bond_table_cols_set = Set(_bond_table_cols)
const _bond_table_cols_priv = Set([:properties, :flags])

@auto_hash_equals struct _BondTable <: AbstractColumnTable
    idx::Vector{Int}
    a1::Vector{Int}
    a2::Vector{Int}
    order::Vector{BondOrderType}

    # private columns
    properties::Vector{Properties}
    flags::Vector{Flags}

    # internals
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

@inline Tables.columnnames(::_BondTable) = _bond_table_cols
@inline Tables.schema(::_BondTable) = _bond_table_schema

@inline function Tables.getcolumn(bt::_BondTable, nm::Symbol)
    @assert nm in _bond_table_cols_priv || nm in _bond_table_cols "type _BondTable has no column $nm"
    getfield(bt, nm)
end

@inline Base.size(bt::_BondTable) = (length(bt.idx), length(_bond_table_cols))

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

@inline Tables.getcolumn(btr::_BondTableRow, nm::Symbol) = Tables.getcolumn(getfield(btr, :_tab), nm)[getfield(btr, :_row)]
@inline Tables.getcolumn(btr::_BondTableRow, i::Int) = getfield(btr, Tables.columnnames(btr)[i])
@inline Tables.columnnames(::_BondTableRow) = _bond_table_cols

@inline _row_by_idx(bt::_BondTable, idx::Int) = _BondTableRow(getfield(bt, :_idx_map)[idx], bt)

@inline function Base.getproperty(btr::_BondTableRow, nm::Symbol)
    getindex(getfield(getfield(btr, :_tab), nm), getfield(btr, :_row))
end

@inline function Base.setproperty!(btr::_BondTableRow, nm::Symbol, val)
    setindex!(getproperty(getfield(btr, :_tab), nm), val, getfield(btr, :_row))
end

@inline Base.eltype(::_BondTable) = _BondTableRow
@inline Base.iterate(bt::_BondTable, st=1) = st > length(bt) ? nothing : (_BondTableRow(st, bt), st + 1)
