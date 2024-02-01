export
    BondTable,
    BondTableRow

const _bond_table_cols = fieldnames(BondTuple)
const _bond_table_cols_set = Set(_bond_table_cols)

@auto_hash_equals struct BondTable <: Tables.AbstractColumns
    idx::Vector{Int}
    a1::Vector{Int}
    a2::Vector{Int}
    order::Vector{BondOrderType}
    properties::Vector{Properties}
    flags::Vector{Flags}

    _idx_map::Dict{Int,Int}

    function BondTable()
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

Tables.istable(::Type{<: BondTable}) = true
Tables.columnaccess(::Type{<: BondTable}) = true
Tables.columns(bt::BondTable) = bt

@inline function Tables.getcolumn(bt::BondTable, nm::Symbol)
    @assert nm in _bond_table_cols "type BondTable has no column $nm"
    getfield(bt, nm)
end
Base.getproperty(bt::BondTable, nm::Symbol) = getfield(bt, nm)

Tables.getcolumn(bt::BondTable, i::Int) = getfield(bt, Tables.columnnames(bt)[i])
Tables.columnnames(::BondTable) = _bond_table_cols
Tables.schema(::BondTable) = Tables.Schema(fieldnames(BondTuple), fieldtypes(BondTuple))

Base.size(bt::BondTable) = (length(bt.idx), length(_bond_table_cols))
Base.size(bt::BondTable, dim) = size(bt)[dim]
Base.length(bt::BondTable) = size(bt, 1)

function Base.push!(bt::BondTable, t::BondTuple)
    getfield(bt, :_idx_map)[t.idx] = length(bt.idx) + 1
    for fn in _bond_table_cols
        push!(getfield(bt, Symbol(fn)), getfield(t, Symbol(fn)))
    end
    bt
end

function _bond_table(itr)
    bt = BondTable()
    for b in itr
        push!(bt, BondTuple(b.a1, b.a2, b.order;
            idx = b.idx,
            properties = b.properties,
            flags = b.flags
        ))
    end
    bt
end
Tables.materializer(::Type{BondTable}) = _bond_table

struct BondTableRow <: Tables.AbstractRow
    _row::Int
    _tab::BondTable
end

Tables.rowaccess(::Type{<: BondTable}) = true
Tables.rows(bt::BondTable) = bt
Tables.getcolumn(btr::BondTableRow, nm::Symbol) = Tables.getcolumn(getfield(btr, :_tab), nm)[getfield(btr, :_row)]
Tables.getcolumn(btr::BondTableRow, i::Int) = getfield(btr, Tables.columnnames(btr)[i])
Tables.columnnames(::BondTableRow) = _bond_table_cols

_row_by_idx(bt::BondTable, idx::Int) = BondTableRow(getfield(bt, :_idx_map)[idx], bt)

function Base.getproperty(btr::BondTableRow, nm::Symbol)
    nm === :idx        && return _getproperty(btr, :idx)::Int
    nm === :a1         && return _getproperty(btr, :a1)::Int
    nm === :a2         && return _getproperty(btr, :a2)::Int
    nm === :order      && return _getproperty(btr, :order)::BondOrderType
    nm === :properties && return _getproperty(btr, :properties)::Properties
    nm === :flags      && return _getproperty(btr, :flags)::Flags
    getindex(getfield(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

@inline function _getproperty(btr::BondTableRow, nm::Symbol)
    getindex(getproperty(getfield(btr, :_tab), nm), getfield(btr, :_row))
end

Base.setproperty!(btr::BondTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(btr, :_tab), nm), val, getfield(btr, :_row))

Base.eltype(::BondTable) = BondTableRow
Base.iterate(bt::BondTable, st=1) = st > length(bt) ? nothing : (BondTableRow(st, bt), st + 1)
