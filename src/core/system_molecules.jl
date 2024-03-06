const _molecule_table_schema = Tables.Schema(
    (:idx, :name, :properties, :flags),
    (Int, String, Properties, Flags)
)
const _molecule_table_cols = _molecule_table_schema.names
const _molecule_table_cols_set = Set(_molecule_table_cols)

@auto_hash_equals struct _MoleculeTable <: Tables.AbstractColumns
    idx::Vector{Int}
    name::Vector{String}
    properties::Vector{Properties}
    flags::Vector{Flags}

    _idx_map::Dict{Int,Int}

    function _MoleculeTable()
        new(
            Int[],
            String[],
            Properties[],
            Flags[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.istable(::Type{<: _MoleculeTable}) = true
@inline Tables.columnaccess(::Type{<: _MoleculeTable}) = true
@inline Tables.columns(mt::_MoleculeTable) = mt

@inline function Tables.getcolumn(mt::_MoleculeTable, nm::Symbol)
    @assert nm in _molecule_table_cols "type _MoleculeTable has no column $nm"
    getfield(mt, nm)
end
@inline Base.getproperty(mt::_MoleculeTable, nm::Symbol) = getfield(mt, nm)

@inline Tables.getcolumn(mt::_MoleculeTable, i::Int) = getfield(mt, Tables.columnnames(mt)[i])
@inline Tables.columnnames(::_MoleculeTable) = _molecule_table_cols
@inline Tables.schema(::_MoleculeTable) = _molecule_table_schema

@inline Base.size(mt::_MoleculeTable) = (length(mt.idx), length(_molecule_table_cols))
@inline Base.size(mt::_MoleculeTable, dim) = size(mt)[dim]
@inline Base.length(mt::_MoleculeTable) = size(mt, 1)

function Base.push!(
    mt::_MoleculeTable,
    idx::Int = 0;
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    mt._idx_map[idx] = length(mt.idx) + 1
    push!(mt.idx, idx)
    push!(mt.name, name)
    push!(mt.properties, properties)
    push!(mt.flags, flags)
    mt
end

function _molecule_table(itr)
    mt = _MoleculeTable()
    for m in itr
        push!(mt, m.idx;
            name = m.name,
            properties = m.properties,
            flags = m.flags
        )
    end
    mt
end
Tables.materializer(::Type{_MoleculeTable}) = _molecule_table

@auto_hash_equals struct _MoleculeTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_MoleculeTable
end

@inline Tables.rowaccess(::Type{<: _MoleculeTable}) = true
@inline Tables.rows(mt::_MoleculeTable) = mt
@inline Tables.getcolumn(mtr::_MoleculeTableRow, nm::Symbol) = Tables.getcolumn(getfield(mtr, :_tab), nm)[getfield(mtr, :_row)]
@inline Tables.getcolumn(mtr::_MoleculeTableRow, i::Int) = getfield(mtr, Tables.columnnames(mtr)[i])
@inline Tables.columnnames(::_MoleculeTableRow) = _molecule_table_cols

@inline _row_by_idx(mt::_MoleculeTable, idx::Int) = _MoleculeTableRow(getfield(mt, :_idx_map)[idx], mt)

@inline function Base.getproperty(mtr::_MoleculeTableRow, nm::Symbol)
    nm === :idx        && return _getproperty(mtr, :idx)::Int
    nm === :name       && return _getproperty(mtr, :name)::String
    nm === :properties && return _getproperty(mtr, :properties)::Properties
    nm === :flags      && return _getproperty(mtr, :flags)::Flags
    getindex(getfield(getfield(mtr, :_tab), nm), getfield(mtr, :_row))
end

@inline function _getproperty(mtr::_MoleculeTableRow, nm::Symbol)
    getindex(getproperty(getfield(mtr, :_tab), nm), getfield(mtr, :_row))
end

@inline Base.setproperty!(mtr::_MoleculeTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(mtr, :_tab), nm), val, getfield(mtr, :_row))

@inline Base.eltype(::_MoleculeTable) = _MoleculeTableRow
@inline Base.iterate(mt::_MoleculeTable, st=1) = st > length(mt) ? nothing : (_MoleculeTableRow(st, mt), st + 1)
