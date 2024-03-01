@auto_hash_equals struct _Residue
    idx::Int
    number::Int
    type::AminoAcid
    properties::Properties
    flags::Flags

    function _Residue(
        number::Int,
        type::AminoAcid;
        idx::Int = 0,
        properties::Properties = Properties(),
        flags::Flags = Flags()
    )
        new(idx, number, type, properties, flags)
    end
end

const _residue_table_cols = fieldnames(_Residue)
const _residue_table_cols_set = Set(_residue_table_cols)
const _residue_table_cols_priv = Set([:molecule_idx, :chain_idx])

@auto_hash_equals struct _ResidueTable <: Tables.AbstractColumns
    # public columns
    idx::Vector{Int}
    number::Vector{Int}
    type::Vector{AminoAcid}
    properties::Vector{Properties}
    flags::Vector{Flags}

    # private columns
    molecule_idx::Vector{Int}
    chain_idx::Vector{Int}

    # internals
    _idx_map::Dict{Int,Int}

    function _ResidueTable()
        new(
            Int[],
            Int[],
            AminoAcid[],
            Properties[],
            Flags[],
            Int[],
            Int[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.istable(::Type{<: _ResidueTable}) = true
@inline Tables.columnaccess(::Type{<: _ResidueTable}) = true
@inline Tables.columns(rt::_ResidueTable) = rt

@inline function Tables.getcolumn(rt::_ResidueTable, nm::Symbol)
    @assert nm in _residue_table_cols_priv || nm in _residue_table_cols_set "type _ResidueTable has no column $nm"
    getfield(rt, nm)
end
@inline Base.getproperty(rt::_ResidueTable, nm::Symbol) = getfield(rt, nm)

@inline Tables.getcolumn(rt::_ResidueTable, i::Int) = getfield(rt, Tables.columnnames(rt)[i])
@inline Tables.columnnames(::_ResidueTable) = _residue_table_cols
@inline Tables.schema(::_ResidueTable) = Tables.Schema(fieldnames(_Residue), fieldtypes(_Residue))

@inline Base.size(rt::_ResidueTable) = (length(rt.idx), length(_residue_table_cols))
@inline Base.size(rt::_ResidueTable, dim) = size(rt)[dim]
@inline Base.length(rt::_ResidueTable) = size(rt, 1)

function Base.push!(rt::_ResidueTable, t::_Residue, molecule_idx::Int, chain_idx::Int)
    getfield(rt, :_idx_map)[t.idx] = length(rt.idx) + 1
    for fn in _residue_table_cols
        push!(getfield(rt, Symbol(fn)), getfield(t, Symbol(fn)))
    end
    push!(getfield(rt, :molecule_idx), molecule_idx)
    push!(getfield(rt, :chain_idx), chain_idx)
    rt
end

function _residue_table(itr)
    rt = _ResidueTable()
    for r in itr
        push!(rt, _Residue(r.number, r.type;
                idx = r.idx,
                properties = r.properties,
                flags = r.flags
            ),
            Tables.getcolumn(r, :molecule_idx),
            Tables.getcolumn(r, :chain_idx)
       )
    end
    rt
end
@inline Tables.materializer(::Type{_ResidueTable}) = itr -> _residue_table(itr)

@auto_hash_equals struct _ResidueTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_ResidueTable
end

@inline Tables.rowaccess(::Type{<: _ResidueTable}) = true
@inline Tables.rows(rt::_ResidueTable) = rt
@inline Tables.getcolumn(rtr::_ResidueTableRow, nm::Symbol) = Tables.getcolumn(getfield(rtr, :_tab), nm)[getfield(rtr, :_row)]
@inline Tables.getcolumn(rtr::_ResidueTableRow, i::Int) = getfield(rtr, Tables.columnnames(rtr)[i])
@inline Tables.columnnames(::_ResidueTableRow) = _residue_table_cols

@inline _row_by_idx(rt::_ResidueTable, idx::Int) = _ResidueTableRow(getfield(rt, :_idx_map)[idx], rt)

@inline function Base.getproperty(rtr::_ResidueTableRow, nm::Symbol)
    nm === :idx         && return _getproperty(rtr, :idx)::Int
    nm === :number      && return _getproperty(rtr, :number)::Int
    nm === :type        && return _getproperty(rtr, :type)::AminoAcid
    nm === :properties  && return _getproperty(rtr, :properties)::Properties
    nm === :flags       && return _getproperty(rtr, :flags)::Flags
    nm === :molecule_idx && return _getproperty(rtr, :molecule_idx)::Int
    nm === :chain_idx    && return _getproperty(rtr, :chain_idx)::Int
    getindex(getfield(getfield(rtr, :_tab), nm), getfield(rtr, :_row))
end

@inline function _getproperty(rtr::_ResidueTableRow, nm::Symbol)
    getindex(getproperty(getfield(rtr, :_tab), nm), getfield(rtr, :_row))
end

@inline Base.setproperty!(rtr::_ResidueTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(rtr, :_tab), nm), val, getfield(rtr, :_row))

@inline Base.eltype(::_ResidueTable) = _ResidueTableRow
@inline Base.iterate(rt::_ResidueTable, st=1) = st > length(rt) ? nothing : (_ResidueTableRow(st, rt), st + 1)
