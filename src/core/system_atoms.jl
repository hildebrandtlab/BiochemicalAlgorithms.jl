const _atom_table_cols = fieldnames(AtomTuple)
const _atom_table_cols_set = Set(_atom_table_cols)
const _atom_table_cols_priv = Set([:frame_id, :molecule_id, :chain_id, :fragment_id, :nucleotide_id, :residue_id])

@auto_hash_equals struct _AtomTable{T <: Real} <: Tables.AbstractColumns
    # public columns
    idx::Vector{Int}
    number::Vector{Int}
    element::Vector{ElementType}
    name::Vector{String}
    atom_type::Vector{String}
    r::Vector{Vector3{T}}
    v::Vector{Vector3{T}}
    F::Vector{Vector3{T}}
    formal_charge::Vector{Int}
    charge::Vector{T}
    radius::Vector{T}
    properties::Vector{Properties}
    flags::Vector{Flags}

    # private columns
    frame_id::Vector{Int}
    molecule_id::Vector{MaybeInt}
    chain_id::Vector{MaybeInt}
    fragment_id::Vector{MaybeInt}
    nucleotide_id::Vector{MaybeInt}
    residue_id::Vector{MaybeInt}

    # internals
    _idx_map::Dict{Int,Int}

    function _AtomTable{T}() where T
        new(
            Int[],
            Int[],
            ElementType[],
            String[],
            String[],
            Vector3{T}[],
            Vector3{T}[],
            Vector3{T}[],
            Int[],
            T[],
            T[],
            Properties[],
            Flags[],
            Int[],
            MaybeInt[],
            MaybeInt[],
            MaybeInt[],
            MaybeInt[],
            MaybeInt[],
            Dict{Int,Int}()
        )
    end
end

@inline Tables.istable(::Type{<: _AtomTable}) = true
@inline Tables.columnaccess(::Type{<: _AtomTable}) = true
@inline Tables.columns(at::_AtomTable) = at

@inline function Tables.getcolumn(at::_AtomTable, nm::Symbol)
    @assert nm in _atom_table_cols_priv || nm in _atom_table_cols_set "type _AtomTable has no column $nm"
    getfield(at, nm)
end
@inline Base.getproperty(at::_AtomTable{T}, nm::Symbol) where T = getfield(at, nm)

@inline Tables.getcolumn(at::_AtomTable, i::Int) = getfield(at, Tables.columnnames(at)[i])
@inline Tables.columnnames(::_AtomTable) = _atom_table_cols
@inline Tables.schema(::_AtomTable{T}) where T = Tables.Schema(fieldnames(AtomTuple{T}), fieldtypes(AtomTuple{T}))

@inline Base.size(at::_AtomTable) = (length(at.idx), length(_atom_table_cols))
@inline Base.size(at::_AtomTable, dim) = size(at)[dim]
@inline Base.length(at::_AtomTable) = size(at, 1)

function Base.push!(at::_AtomTable{T}, t::AtomTuple{T};
    frame_id::Int = 1,
    molecule_id::MaybeInt = nothing,
    chain_id::MaybeInt = nothing,
    fragment_id::MaybeInt = nothing,
    nucleotide_id::MaybeInt = nothing,
    residue_id::MaybeInt = nothing
) where T
    getfield(at, :_idx_map)[t.idx] = length(at.idx) + 1
    for fn in _atom_table_cols
        push!(getfield(at, Symbol(fn)), getfield(t, Symbol(fn)))
    end
    push!(getfield(at, :frame_id), frame_id)
    push!(getfield(at, :molecule_id), molecule_id)
    push!(getfield(at, :chain_id), chain_id)
    push!(getfield(at, :fragment_id), fragment_id)
    push!(getfield(at, :nucleotide_id), nucleotide_id)
    push!(getfield(at, :residue_id), residue_id)
    at
end

function _atom_table(::Type{T}, itr) where T
    at = _AtomTable{T}()
    for a in itr
        push!(at, AtomTuple{T}(a.number, a.element;
                idx = a.idx,
                name = a.name,
                atom_type = a.atom_type,
                r = a.r,
                v = a.v,
                F = a.F,
                formal_charge = a.formal_charge,
                charge = a.charge,
                radius = a.radius,
                properties = a.properties,
                flags = a.flags
            );
            frame_id = Tables.getcolumn(a, :frame_id),
            molecule_id = Tables.getcolumn(a, :molecule_id),
            chain_id = Tables.getcolumn(a, :chain_id),
            fragment_id = Tables.getcolumn(a, :fragment_id),
            nucleotide_id = Tables.getcolumn(a, :nucleotide_id),
            residue_id = Tables.getcolumn(a, :residue_id)
       )
    end
    at
end
@inline Tables.materializer(::Type{_AtomTable{T}}) where T = itr -> _atom_table(T, itr)

@inline function _filter_select(itr, col::Symbol)
    (getproperty(a, col) for a in itr)
end

@auto_hash_equals struct _AtomTableRow{T} <: Tables.AbstractRow
    _row::Int
    _tab::_AtomTable{T}
end

@inline Tables.rowaccess(::Type{<: _AtomTable}) = true
@inline Tables.rows(at::_AtomTable) = at

@inline Tables.getcolumn(atr::_AtomTableRow, nm::Symbol) = Tables.getcolumn(getfield(atr, :_tab), nm)[getfield(atr, :_row)]
@inline Tables.getcolumn(atr::_AtomTableRow, i::Int) = getfield(atr, Tables.columnnames(atr)[i])
@inline Tables.columnnames(::_AtomTableRow) = _atom_table_cols

@inline _row_by_idx(at::_AtomTable{T}, idx::Int) where T = _AtomTableRow{T}(getfield(at, :_idx_map)[idx], at)

@inline function Base.getproperty(atr::_AtomTableRow{T}, nm::Symbol) where T
    nm === :idx           && return _getproperty(atr, :idx)::Int
    nm === :number        && return _getproperty(atr, :number)::Int
    nm === :element       && return _getproperty(atr, :element)::ElementType
    nm === :name          && return _getproperty(atr, :name)::String
    nm === :atom_type     && return _getproperty(atr, :atom_type)::String
    nm === :r             && return _getproperty(atr, :r)::Vector3{T}
    nm === :v             && return _getproperty(atr, :v)::Vector3{T}
    nm === :F             && return _getproperty(atr, :F)::Vector3{T}
    nm === :properties    && return _getproperty(atr, :properties)::Properties
    nm === :flags         && return _getproperty(atr, :flags)::Flags
    nm === :frame_id      && return _getproperty(atr, :frame_id)::Int
    nm === :molecule_id   && return _getproperty(atr, :molecule_id)::MaybeInt
    nm === :chain_id      && return _getproperty(atr, :chain_id)::MaybeInt
    nm === :fragment_id   && return _getproperty(atr, :fragment_id)::MaybeInt
    nm === :nucleotide_id && return _getproperty(atr, :nucleotide_id)::MaybeInt
    nm === :residue_id    && return _getproperty(atr, :residue_id)::MaybeInt
    getindex(getfield(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

@inline function _getproperty(atr::_AtomTableRow{T}, nm::Symbol) where T
    getindex(getproperty(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

@inline Base.setproperty!(atr::_AtomTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(atr, :_tab), nm), val, getfield(atr, :_row))

@inline Base.eltype(::_AtomTable{T}) where T = _AtomTableRow{T}
@inline Base.iterate(at::_AtomTable, st=1) = st > length(at) ? nothing : (_AtomTableRow(st, at), st + 1)
