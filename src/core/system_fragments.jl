@auto_hash_equals struct _Fragment
    idx::Int
    number::Int
    name::String
    properties::Properties
    flags::Flags

    function _Fragment(
        number::Int;
        idx::Int = 0,
        name::String = "",
        properties::Properties = Properties(),
        flags::Flags = Flags()
    )
        new(idx, number, name, properties, flags)
    end
end

const _fragment_table_cols = fieldnames(_Fragment)
const _fragment_table_cols_set = Set(_fragment_table_cols)
const _fragment_table_cols_priv = Set([:molecule_id, :chain_id])

@auto_hash_equals struct _FragmentTable <: Tables.AbstractColumns
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

    function _FragmentTable()
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

@inline Tables.istable(::Type{<: _FragmentTable}) = true
@inline Tables.columnaccess(::Type{<: _FragmentTable}) = true
@inline Tables.columns(ft::_FragmentTable) = ft

@inline function Tables.getcolumn(ft::_FragmentTable, nm::Symbol)
    @assert nm in _fragment_table_cols_priv || nm in _fragment_table_cols_set "type _FragmentTable has no column $nm"
    getfield(ft, nm)
end
@inline Base.getproperty(ft::_FragmentTable, nm::Symbol) = getfield(ft, nm)

@inline Tables.getcolumn(ft::_FragmentTable, i::Int) = getfield(ft, Tables.columnnames(ft)[i])
@inline Tables.columnnames(::_FragmentTable) = _fragment_table_cols
@inline Tables.schema(::_FragmentTable) = Tables.Schema(fieldnames(_Fragment), fieldtypes(_Fragment))

@inline Base.size(ft::_FragmentTable) = (length(ft.idx), length(_fragment_table_cols))
@inline Base.size(ft::_FragmentTable, dim) = size(ft)[dim]
@inline Base.length(ft::_FragmentTable) = size(ft, 1)

function Base.push!(ft::_FragmentTable, t::_Fragment, molecule_id::Int, chain_id::Int)
    getfield(ft, :_idx_map)[t.idx] = length(ft.idx) + 1
    for fn in _fragment_table_cols
        push!(getfield(ft, Symbol(fn)), getfield(t, Symbol(fn)))
    end
    push!(getfield(ft, :molecule_id), molecule_id)
    push!(getfield(ft, :chain_id), chain_id)
    ft
end

function _fragment_table(itr)
    ft = _FragmentTable()
    for f in itr
        push!(ft, _Fragment(f.number;
                idx = f.idx,
                name = f.name,
                properties = f.properties,
                flags = f.flags
            );
            molecule_id = Tables.getcolumn(f, :molecule_id),
            chain_id = Tables.getcolumn(f, :chain_id)
       )
    end
    ft
end
@inline Tables.materializer(::Type{_FragmentTable}) = itr -> _fragment_table(itr)

@auto_hash_equals struct _FragmentTableRow <: Tables.AbstractRow
    _row::Int
    _tab::_FragmentTable
end

@inline Tables.rowaccess(::Type{<: _FragmentTable}) = true
@inline Tables.rows(ft::_FragmentTable) = ft
@inline Tables.getcolumn(ftr::_FragmentTableRow, nm::Symbol) = Tables.getcolumn(getfield(ftr, :_tab), nm)[getfield(ftr, :_row)]
@inline Tables.getcolumn(ftr::_FragmentTableRow, i::Int) = getfield(ftr, Tables.columnnames(ftr)[i])
@inline Tables.columnnames(::_FragmentTableRow) = _fragment_table_cols

@inline _row_by_idx(ft::_FragmentTable, idx::Int) = _FragmentTableRow(getfield(ft, :_idx_map)[idx], ft)

@inline function Base.getproperty(ftr::_FragmentTableRow, nm::Symbol)
    nm === :idx         && return _getproperty(ftr, :idx)::Int
    nm === :number      && return _getproperty(ftr, :number)::Int
    nm === :name        && return _getproperty(ftr, :name)::String
    nm === :properties  && return _getproperty(ftr, :properties)::Properties
    nm === :flags       && return _getproperty(ftr, :flags)::Flags
    nm === :molecule_id && return _getproperty(ftr, :molecule_id)::Int
    nm === :chain_id    && return _getproperty(ftr, :chain_id)::Int
    getindex(getfield(getfield(ftr, :_tab), nm), getfield(ftr, :_row))
end

@inline function _getproperty(ftr::_FragmentTableRow, nm::Symbol)
    getindex(getproperty(getfield(ftr, :_tab), nm), getfield(ftr, :_row))
end

@inline Base.setproperty!(ftr::_FragmentTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(ftr, :_tab), nm), val, getfield(ftr, :_row))

@inline Base.eltype(::_FragmentTable) = _FragmentTableRow
@inline Base.iterate(ft::_FragmentTable, st=1) = st > length(ft) ? nothing : (_FragmentTableRow(st, ft), st + 1)
