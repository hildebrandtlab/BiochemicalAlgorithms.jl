const _bond_table_cols_main = (:idx, :order)
const _bond_table_cols_extra = (:properties, :flags, :atom1_idx, :atom2_idx)
const _bond_table_cols = (_bond_table_cols_main..., _bond_table_cols_extra...)
const _bond_table_cols_set = Set(_bond_table_cols)
const _bond_table_schema = Tables.Schema(
    _bond_table_cols_main,
    (Int, Int, Int, BondOrderType)
)

@auto_hash_equals struct _BondTable <: _AbstractSystemComponentTable
    idx::Vector{Int}
    order::Vector{BondOrderType}

    # private columns
    properties::Vector{Properties}
    flags::Vector{Flags}
    atom1_idx::Vector{Int}
    atom2_idx::Vector{Int}

    # internals
    _idx_map::_IdxMap

    function _BondTable()
        new(
            Int[],
            BondOrderType[],
            Properties[],
            Flags[],
            Int[],
            Int[],
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
    order::BondOrderType,
    atom1_idx::Int,
    atom2_idx::Int;
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    bt._idx_map[idx] = length(bt.idx) + 1
    push!(bt.idx, idx)
    push!(bt.order, order)
    push!(bt.properties, properties)
    push!(bt.flags, flags)
    push!(bt.atom1_idx, atom1_idx)
    push!(bt.atom2_idx, atom2_idx)
    bt
end

function _bond_table(itr)
    bt = _BondTable()
    for b in itr
        push!(bt, b.idx, b.order, b.atom1_idx, b.atom2_idx;
            properties = b.properties,
            flags = b.flags
        )
    end
    bt
end
Tables.materializer(::Type{_BondTable}) = _bond_table
