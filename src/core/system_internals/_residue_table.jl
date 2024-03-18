const _residue_table_schema = Tables.Schema(
    (:idx, :number, :type),
    (Int, Int, AminoAcid)
)
const _residue_table_cols = _residue_table_schema.names
const _residue_table_cols_set = Set(_residue_table_cols)
const _residue_table_cols_priv = Set([:properties, :flags, :molecule_idx, :chain_idx])

@auto_hash_equals struct _ResidueTable <: AbstractColumnTable
    # public columns
    idx::Vector{Int}
    number::Vector{Int}
    type::Vector{AminoAcid}

    # private columns
    properties::Vector{Properties}
    flags::Vector{Flags}
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

@inline Tables.columnnames(::_ResidueTable) = _residue_table_cols
@inline Tables.schema(::_ResidueTable) = _residue_table_schema

@inline function Tables.getcolumn(rt::_ResidueTable, nm::Symbol)
    @assert nm in _residue_table_cols_priv || nm in _residue_table_cols_set "type _ResidueTable has no column $nm"
    getfield(rt, nm)
end

@inline Base.size(rt::_ResidueTable) = (length(rt.idx), length(_residue_table_cols))

function Base.push!(
    rt::_ResidueTable,
    idx::Int,
    number::Int,
    type::AminoAcid,
    molecule_idx::Int,
    chain_idx::Int;
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    rt._idx_map[idx] = length(rt.idx) + 1
    push!(rt.idx, idx)
    push!(rt.number, number)
    push!(rt.type, type)
    push!(rt.properties, properties)
    push!(rt.flags, flags)
    push!(rt.molecule_idx, molecule_idx)
    push!(rt.chain_idx, chain_idx)
    rt
end

function _residue_table(itr)
    rt = _ResidueTable()
    for r in itr
        push!(rt, r.idx, r.number, r.type, r.molecule_idx, r.chain_idx;
            properties = r.properties,
            flags = r.flags
        )
    end
    rt
end
@inline Tables.materializer(::Type{_ResidueTable}) = itr -> _residue_table(itr)

@auto_hash_equals struct _ResidueTableRow <: _AbstractColumnTableRow
    _row::Int
    _tab::_ResidueTable
end

@inline Tables.getcolumn(rtr::_ResidueTableRow, nm::Symbol) = Tables.getcolumn(getfield(rtr, :_tab), nm)[getfield(rtr, :_row)]
@inline Tables.getcolumn(rtr::_ResidueTableRow, i::Int) = getfield(rtr, Tables.columnnames(rtr)[i])
@inline Tables.columnnames(::_ResidueTableRow) = _residue_table_cols

@inline _row_by_idx(rt::_ResidueTable, idx::Int) = _ResidueTableRow(getfield(rt, :_idx_map)[idx], rt)

@inline function Base.getproperty(rtr::_ResidueTableRow, nm::Symbol)
    getindex(getfield(getfield(rtr, :_tab), nm), getfield(rtr, :_row))
end

@inline function Base.setproperty!(rtr::_ResidueTableRow, nm::Symbol, val) 
    setindex!(getproperty(getfield(rtr, :_tab), nm), val, getfield(rtr, :_row))
end

@inline Base.eltype(::_ResidueTable) = _ResidueTableRow
@inline Base.iterate(rt::_ResidueTable, st=1) = st > length(rt) ? nothing : (_ResidueTableRow(st, rt), st + 1)
