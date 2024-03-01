@auto_hash_equals struct _Nucleotide
    idx::Int
    number::Int
    name::String
    properties::Properties
    flags::Flags

    function _Nucleotide(
        number::Int;
        idx::Int = 0,
        name::String = "",
        properties::Properties = Properties(),
        flags::Flags = Flags()
    )
        new(idx, number, name, properties, flags)
    end
end

const _nucleotide_table_cols = fieldnames(_Nucleotide)
const _nucleotide_table_cols_set = Set(_nucleotide_table_cols)
const _nucleotide_table_cols_priv = Set([:molecule_id, :chain_id])

@auto_hash_equals struct _NucleotideTable <: Tables.AbstractColumns
    # public columns
    idx::Vector{Int}
    number::Vector{Int}
    name::Vector{String}
    properties::Vector{Properties}
    flags::Vector{Flags}

    # private columns
    molecule_id::Vector{Int}
    chain_id::Vector{Int}

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

@inline Tables.istable(::Type{<: _NucleotideTable}) = true
@inline Tables.columnaccess(::Type{<: _NucleotideTable}) = true
@inline Tables.columns(nt::_NucleotideTable) = nt

@inline function Tables.getcolumn(nt::_NucleotideTable, nm::Symbol)
    @assert nm in _nucleotide_table_cols_priv || nm in _nucleotide_table_cols_set "type _NucleotideTable has no column $nm"
    getfield(nt, nm)
end
@inline Base.getproperty(nt::_NucleotideTable, nm::Symbol) = getfield(nt, nm)

@inline Tables.getcolumn(nt::_NucleotideTable, i::Int) = getfield(nt, Tables.columnnames(nt)[i])
@inline Tables.columnnames(::_NucleotideTable) = _nucleotide_table_cols
@inline Tables.schema(::_NucleotideTable) = Tables.Schema(fieldnames(_Nucleotide), fieldtypes(_Nucleotide))

@inline Base.size(nt::_NucleotideTable) = (length(nt.idx), length(_nucleotide_table_cols))
@inline Base.size(nt::_NucleotideTable, dim) = size(nt)[dim]
@inline Base.length(nt::_NucleotideTable) = size(nt, 1)

function Base.push!(nt::_NucleotideTable, t::_Nucleotide, molecule_id::Int, chain_id::Int)
    getfield(nt, :_idx_map)[t.idx] = length(nt.idx) + 1
    for fn in _nucleotide_table_cols
        push!(getfield(nt, Symbol(fn)), getfield(t, Symbol(fn)))
    end
    push!(getfield(nt, :molecule_id), molecule_id)
    push!(getfield(nt, :chain_id), chain_id)
    nt
end

function _nucleotide_table(itr)
    nt = _NucleotideTable()
    for n in itr
        push!(nt, _Nucleotide(n.number;
                idx = n.idx,
                name = n.name,
                properties = n.properties,
                flags = n.flags
            );
            molecule_id = Tables.getcolumn(n, :molecule_id),
            chain_id = Tables.getcolumn(n, :chain_id)
       )
    end
    nt
end
@inline Tables.materializer(::Type{_NucleotideTable}) = itr -> _nucleotide_table(itr)

@auto_hash_equals struct _NucleotideTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_NucleotideTable
end

@inline Tables.rowaccess(::Type{<: _NucleotideTable}) = true
@inline Tables.rows(nt::_NucleotideTable) = nt
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
    nm === :molecule_id && return _getproperty(ntr, :molecule_id)::Int
    nm === :chain_id    && return _getproperty(ntr, :chain_id)::Int
    getindex(getfield(getfield(ntr, :_tab), nm), getfield(ntr, :_row))
end

@inline function _getproperty(ntr::_NucleotideTableRow, nm::Symbol)
    getindex(getproperty(getfield(ntr, :_tab), nm), getfield(ntr, :_row))
end

@inline Base.setproperty!(ntr::_NucleotideTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(ntr, :_tab), nm), val, getfield(ntr, :_row))

@inline Base.eltype(::_NucleotideTable) = _NucleotideTableRow
@inline Base.iterate(nt::_NucleotideTable, st=1) = st > length(nt) ? nothing : (_NucleotideTableRow(st, nt), st + 1)
