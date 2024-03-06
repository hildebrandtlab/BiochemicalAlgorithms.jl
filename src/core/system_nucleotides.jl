const _nucleotide_table_schema = Tables.Schema(
    (:idx, :number, :name, :properties, :flags),
    (Int, Int, String, Properties, Flags)
)
const _nucleotide_table_cols = _nucleotide_table_schema.names
const _nucleotide_table_cols_set = Set(_nucleotide_table_cols)
const _nucleotide_table_cols_priv = Set([:molecule_idx, :chain_idx])

@auto_hash_equals struct _NucleotideTable <: AbstractColumnTable
    # public columns
    idx::Vector{Int}
    number::Vector{Int}
    name::Vector{String}
    properties::Vector{Properties}
    flags::Vector{Flags}

    # private columns
    molecule_idx::Vector{Int}
    chain_idx::Vector{Int}

    # internals
    _idx_map::Dict{Int,Int}

    function _NucleotideTable()
        new(
            Int[],
            Int[],
            String[],
            Properties[],
            Flags[],
            Int[],
            Int[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.columnnames(::_NucleotideTable) = _nucleotide_table_cols
@inline Tables.schema(::_NucleotideTable) = _nucleotide_table_schema

@inline function Tables.getcolumn(nt::_NucleotideTable, nm::Symbol)
    @assert nm in _nucleotide_table_cols_priv || nm in _nucleotide_table_cols_set "type _NucleotideTable has no column $nm"
    getfield(nt, nm)
end

@inline Base.size(nt::_NucleotideTable) = (length(nt.idx), length(_nucleotide_table_cols))

function Base.push!(
    nt::_NucleotideTable,
    idx::Int,
    number::Int,
    molecule_idx::Int,
    chain_idx::Int;
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
    nt._idx_map[idx] = length(nt.idx) + 1
    push!(nt.idx, idx)
    push!(nt.number, number)
    push!(nt.name, name)
    push!(nt.properties, properties)
    push!(nt.flags, flags)
    push!(nt.molecule_idx, molecule_idx)
    push!(nt.chain_idx, chain_idx)
    nt
end

function _nucleotide_table(itr)
    nt = _NucleotideTable()
    for n in itr
        push!(nt, n.idx, n.number, n.molecule_idx, n.chain_idx;
            name = n.name,
            properties = n.properties,
            flags = n.flags
       )
    end
    nt
end
@inline Tables.materializer(::Type{_NucleotideTable}) = itr -> _nucleotide_table(itr)

@auto_hash_equals struct _NucleotideTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_NucleotideTable
end

@inline Tables.getcolumn(ntr::_NucleotideTableRow, nm::Symbol) = Tables.getcolumn(getfield(ntr, :_tab), nm)[getfield(ntr, :_row)]
@inline Tables.getcolumn(ntr::_NucleotideTableRow, i::Int) = getfield(ntr, Tables.columnnames(ntr)[i])
@inline Tables.columnnames(::_NucleotideTableRow) = _nucleotide_table_cols

@inline _row_by_idx(nt::_NucleotideTable, idx::Int) = _NucleotideTableRow(getfield(nt, :_idx_map)[idx], nt)

@inline function Base.getproperty(ntr::_NucleotideTableRow, nm::Symbol)
    nm === :idx         && return _getproperty(ntr, :idx)::Int
    nm === :number      && return _getproperty(ntr, :number)::Int
    nm === :name        && return _getproperty(ntr, :name)::String
    nm === :properties  && return _getproperty(ntr, :properties)::Properties
    nm === :flags       && return _getproperty(ntr, :flags)::Flags
    nm === :molecule_idx && return _getproperty(ntr, :molecule_idx)::Int
    nm === :chain_idx    && return _getproperty(ntr, :chain_idx)::Int
    getindex(getfield(getfield(ntr, :_tab), nm), getfield(ntr, :_row))
end

@inline function _getproperty(ntr::_NucleotideTableRow, nm::Symbol)
    getindex(getproperty(getfield(ntr, :_tab), nm), getfield(ntr, :_row))
end

@inline Base.setproperty!(ntr::_NucleotideTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(ntr, :_tab), nm), val, getfield(ntr, :_row))

@inline Base.eltype(::_NucleotideTable) = _NucleotideTableRow
@inline Base.iterate(nt::_NucleotideTable, st=1) = st > length(nt) ? nothing : (_NucleotideTableRow(st, nt), st + 1)
