export
    AtomTable,
    AtomTableRow,
    AtomTableView

const _atom_table_cols = fieldnames(AtomTuple)
const _atom_table_cols_set = Set(_atom_table_cols)
const _atom_table_cols_priv = Set([:frame_id, :molecule_id, :chain_id, :fragment_id, :nucleotide_id, :residue_id])

@auto_hash_equals struct AtomTable{T <: Real} <: Tables.AbstractColumns
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

    function AtomTable{T}() where T
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

Tables.istable(::Type{<: AtomTable}) = true
Tables.columnaccess(::Type{<: AtomTable}) = true
Tables.columns(at::AtomTable) = at

@inline function Tables.getcolumn(at::AtomTable, nm::Symbol)
    @assert nm in _atom_table_cols_priv || nm in _atom_table_cols_set "type AtomTable has no column $nm"
    getfield(at, nm)
end
Base.getproperty(at::AtomTable, nm::Symbol) = getfield(at, nm)

Tables.getcolumn(at::AtomTable, i::Int) = getfield(at, Tables.columnnames(at)[i])
Tables.columnnames(::AtomTable) = _atom_table_cols
Tables.schema(::AtomTable{T}) where T = Tables.Schema(fieldnames(AtomTuple{T}), fieldtypes(AtomTuple{T}))

Base.size(at::AtomTable) = (length(at.idx), length(_atom_table_cols))
Base.size(at::AtomTable, dim) = size(at)[dim]
Base.length(at::AtomTable) = size(at, 1)

function Base.push!(at::AtomTable{T}, t::AtomTuple{T};
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
    at = AtomTable{T}()
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
Tables.materializer(::Type{AtomTable{T}}) where T = itr -> _atom_table(T, itr)

function _filter_select(itr, col::Symbol)
    (getproperty(a, col) for a in itr)
end

@auto_hash_equals struct AtomTableView{T} <: Tables.AbstractColumns
    _tab::AtomTable{T}
    _idx::Vector{Int}
end

Tables.istable(::Type{<: AtomTableView}) = true
Tables.columnaccess(::Type{<: AtomTableView}) = true
Tables.columns(at::AtomTableView) = at

function Tables.getcolumn(at::AtomTableView, nm::Symbol)
    RowProjectionVector(
        Tables.getcolumn(getfield(at, :_tab), nm),
        map(idx -> getfield(at, :_tab)._idx_map[idx], getfield(at, :_idx))
    )
end
Base.getproperty(at::AtomTableView, nm::Symbol) = Tables.getcolumn(at, nm)

Tables.getcolumn(at::AtomTableView, i::Int) = Tables.getcolumn(at, Tables.columnnames(at)[i])
Tables.columnnames(::AtomTableView) = _atom_table_cols
Tables.schema(at::AtomTableView) = Tables.schema(getfield(at, :_tab))

Base.size(at::AtomTableView) = (length(getfield(at, :_idx)), length(_atom_table_cols))
Base.size(at::AtomTableView, dim) = size(at)[dim]
Base.length(at::AtomTableView) = size(at, 1)

function Base.filter(f, at::AtomTable)
    AtomTableView(at, collect(Int, _filter_select(
        TableOperations.filter(f, at),
        :idx
    )))
end

function Base.filter(f, at::AtomTableView)
    AtomTableView(getfield(at, :_tab), collect(Int, _filter_select(
        TableOperations.filter(f, at),
        :idx
    )))
end

struct AtomTableRow{T} <: Tables.AbstractRow
    _row::Int
    _tab::AtomTable{T}
end

Tables.rowaccess(::Type{<: AtomTable}) = true
Tables.rows(at::AtomTable) = at

Tables.rowaccess(::Type{<: AtomTableView}) = true
Tables.rows(at::AtomTableView) = at

Tables.getcolumn(atr::AtomTableRow, nm::Symbol) = Tables.getcolumn(getfield(atr, :_tab), nm)[getfield(atr, :_row)]
Tables.getcolumn(atr::AtomTableRow, i::Int) = getfield(atr, Tables.columnnames(atr)[i])
Tables.columnnames(::AtomTableRow) = _atom_table_cols

_row_by_idx(at::AtomTable{T}, idx::Int) where T = AtomTableRow{T}(getfield(at, :_idx_map)[idx], at)
_row_by_idx(at::AtomTableView, idx::Int) = _row_by_idx(getfield(at, :_tab), idx)

function Base.getproperty(atr::AtomTableRow{T}, nm::Symbol) where T
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

@inline function _getproperty(atr::AtomTableRow{T}, nm::Symbol) where T
    getindex(getproperty(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

Base.setproperty!(atr::AtomTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(atr, :_tab), nm), val, getfield(atr, :_row))

Base.eltype(::AtomTable{T}) where T = AtomTableRow{T}
Base.iterate(at::AtomTable, st=1) = st > length(at) ? nothing : (AtomTableRow(st, at), st + 1)

Base.eltype(::AtomTableView{T}) where T = AtomTableRow{T}
Base.iterate(at::AtomTableView, st=1) = st > length(at) ? nothing : (_row_by_idx(at, getfield(at, :_idx)[st]), st + 1)
