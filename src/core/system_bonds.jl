const _bond_table_cols = fieldnames(BondTuple)
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

Tables.istable(::Type{<: _BondTable}) = true
Tables.columnaccess(::Type{<: _BondTable}) = true
Tables.columns(bt::_BondTable) = bt

@inline function Tables.getcolumn(bt::_BondTable, nm::Symbol)
    @assert nm in _bond_table_cols "type _BondTable has no column $nm"
    getfield(bt, nm)
end
Base.getproperty(bt::_BondTable, nm::Symbol) = getfield(bt, nm)

Tables.getcolumn(bt::_BondTable, i::Int) = getfield(bt, Tables.columnnames(bt)[i])
Tables.columnnames(::_BondTable) = _bond_table_cols
Tables.schema(::_BondTable) = Tables.Schema(fieldnames(BondTuple), fieldtypes(BondTuple))

Base.size(bt::_BondTable) = (length(bt.idx), length(_bond_table_cols))
Base.size(bt::_BondTable, dim) = size(bt)[dim]
Base.length(bt::_BondTable) = size(bt, 1)

function Base.push!(bt::_BondTable, t::BondTuple)
    getfield(bt, :_idx_map)[t.idx] = length(bt.idx) + 1
    for fn in _bond_table_cols
        push!(getfield(bt, Symbol(fn)), getfield(t, Symbol(fn)))
    end
    bt
end

function _bond_table(itr)
    bt = _BondTable()
    for b in itr
        push!(bt, BondTuple(b.a1, b.a2, b.order;
            idx = b.idx,
            properties = b.properties,
            flags = b.flags
        ))
    end
    bt
end
Tables.materializer(::Type{_BondTable}) = _bond_table

struct _BondTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_BondTable
end

Tables.rowaccess(::Type{<: _BondTable}) = true
Tables.rows(bt::_BondTable) = bt
Tables.getcolumn(btr::_BondTableRow, nm::Symbol) = Tables.getcolumn(getfield(btr, :_tab), nm)[getfield(btr, :_row)]
Tables.getcolumn(btr::_BondTableRow, i::Int) = getfield(btr, Tables.columnnames(btr)[i])
Tables.columnnames(::_BondTableRow) = _bond_table_cols

_row_by_idx(bt::_BondTable, idx::Int) = _BondTableRow(getfield(bt, :_idx_map)[idx], bt)

function Base.getproperty(btr::_BondTableRow, nm::Symbol)
    nm === :idx        && return _getproperty(btr, :idx)::Int
    nm === :a1         && return _getproperty(btr, :a1)::Int
    nm === :a2         && return _getproperty(btr, :a2)::Int
    nm === :order      && return _getproperty(btr, :order)::BondOrderType
    nm === :properties && return _getproperty(btr, :properties)::Properties
    nm === :flags      && return _getproperty(btr, :flags)::Flags
    getindex(getfield(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

@inline function _getproperty(btr::_BondTableRow, nm::Symbol)
    getindex(getproperty(getfield(btr, :_tab), nm), getfield(btr, :_row))
end

Base.setproperty!(btr::_BondTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(btr, :_tab), nm), val, getfield(btr, :_row))

Base.eltype(::_BondTable) = _BondTableRow
Base.iterate(bt::_BondTable, st=1) = st > length(bt) ? nothing : (_BondTableRow(st, bt), st + 1)
