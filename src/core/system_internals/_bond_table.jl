const _bond_table_cols_main = (:idx, :a1, :a2, :order)
const _bond_table_cols_extra = (:properties, :flags)
const _bond_table_cols = (_bond_table_cols_main..., _bond_table_cols_extra...)
const _bond_table_cols_set = Set(_bond_table_cols)
const _bond_table_schema = Tables.Schema(
    _bond_table_cols_main,
    (Int, Int, Int, BondOrderType)
)

@auto_hash_equals struct _BondTable <: _AbstractSystemComponentTable
    idx::Vector{Int}
    a1::Vector{Int}
    a2::Vector{Int}
    order::Vector{BondOrderType}

    # private columns
    properties::Vector{Properties}
    flags::Vector{Flags}

    # internals
    _idx_map::_IdxMap

    function _BondTable()
        new(
            Int[],
            Int[],
            Int[],
            BondOrderType[],
            Properties[],
            Flags[],
            _IdxMap()
        )
    end
end

@inline function _hascolumn(::Union{_BondTable, Type{_BondTable}}, nm::Symbol)
    nm in _bond_table_cols_set
end

@inline function Tables.columnnames(::_BondTable)
    _bond_table_cols_main
end

@inline function Tables.schema(::_BondTable)
    _bond_table_schema
end

@inline function Base.propertynames(::_BondTable)
    _bond_table_cols
end

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

function _delete!(bt::_BondTable, rowno::Int)
    deleteat!(bt.idx, rowno)
    deleteat!(bt.a1, rowno)
    deleteat!(bt.a2, rowno)
    deleteat!(bt.order, rowno)
    deleteat!(bt.properties, rowno)
    deleteat!(bt.flags, rowno)
    nothing
end

function Base.empty!(bt::_BondTable)
    empty!(bt.idx)
    empty!(bt.a1)
    empty!(bt.a2)
    empty!(bt.order)
    empty!(bt.properties)
    empty!(bt.flags)
    empty!(bt._idx_map)
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
