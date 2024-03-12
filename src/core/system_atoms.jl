const _atom_table_cols = (:idx, :number, :element, :name, :atom_type, :r, :v, :F, :formal_charge, :charge, :radius)
const _atom_table_cols_set = Set(_atom_table_cols)
const _atom_table_cols_priv = Set([:properties, :flags, :frame_id, :molecule_idx, :chain_idx, :fragment_idx, :nucleotide_idx, :residue_idx])

@auto_hash_equals struct _AtomTable{T <: Real} <: AbstractColumnTable
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

    # private columns
    properties::Vector{Properties}
    flags::Vector{Flags}
    frame_id::Vector{Int}
    molecule_idx::Vector{MaybeInt}
    chain_idx::Vector{MaybeInt}
    fragment_idx::Vector{MaybeInt}
    nucleotide_idx::Vector{MaybeInt}
    residue_idx::Vector{MaybeInt}

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

@inline Tables.columnnames(::_AtomTable) = _atom_table_cols

@inline function Tables.schema(::_AtomTable{T}) where T
    Tables.Schema(
        _atom_table_cols,
        (Int, Int, ElementType, String, String, Vector3{T}, Vector3{T}, Vector3{T}, Int, T, T)
    )
end

@inline function Tables.getcolumn(at::_AtomTable, nm::Symbol)
    @assert nm in _atom_table_cols_priv || nm in _atom_table_cols_set "type _AtomTable has no column $nm"
    getfield(at, nm)
end

@inline Base.size(at::_AtomTable) = (length(at.idx), length(_atom_table_cols))

function Base.push!(
    at::_AtomTable{T},
    idx::Int,
    number::Int,
    element::ElementType;
    name::String = "",
    atom_type::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    properties::Properties = Properties(),
    flags::Flags = Flags(),
    frame_id::Int = 1,
    molecule_idx::MaybeInt = nothing,
    chain_idx::MaybeInt = nothing,
    fragment_idx::MaybeInt = nothing,
    nucleotide_idx::MaybeInt = nothing,
    residue_idx::MaybeInt = nothing
) where T
    at._idx_map[idx] = length(at.idx) + 1
    push!(at.idx, idx)
    push!(at.number, number)
    push!(at.element, element)
    push!(at.name, name)
    push!(at.atom_type, atom_type)
    push!(at.r, r)
    push!(at.v, v)
    push!(at.F, F)
    push!(at.formal_charge, formal_charge)
    push!(at.charge, charge)
    push!(at.radius, radius)
    push!(at.properties, properties)
    push!(at.flags, flags)
    push!(at.frame_id, frame_id)
    push!(at.molecule_idx, molecule_idx)
    push!(at.chain_idx, chain_idx)
    push!(at.fragment_idx, fragment_idx)
    push!(at.nucleotide_idx, nucleotide_idx)
    push!(at.residue_idx, residue_idx)
    at
end

function _atom_table(::Type{T}, itr) where T
    at = _AtomTable{T}()
    for a in itr
        push!(at, a.idx, a.number, a.element;
            name = a.name,
            atom_type = a.atom_type,
            r = a.r,
            v = a.v,
            F = a.F,
            formal_charge = a.formal_charge,
            charge = a.charge,
            radius = a.radius,
            properties = a.properties,
            flags = a.flags,
            frame_id = a.frame_id,
            molecule_idx = a.molecule_idx,
            chain_idx = a.chain_idx,
            fragment_idx = a.fragment_idx,
            nucleotide_idx = a.nucleotide_idx,
            residue_idx = a.residue_idx
        )
    end
    at
end
@inline Tables.materializer(::Type{_AtomTable{T}}) where T = itr -> _atom_table(T, itr)

@auto_hash_equals struct _AtomTableRow{T} <: Tables.AbstractRow
    _row::Int
    _tab::_AtomTable{T}
end

@inline Tables.getcolumn(atr::_AtomTableRow, nm::Symbol) = Tables.getcolumn(getfield(atr, :_tab), nm)[getfield(atr, :_row)]
@inline Tables.getcolumn(atr::_AtomTableRow, i::Int) = getfield(atr, Tables.columnnames(atr)[i])
@inline Tables.columnnames(::_AtomTableRow) = _atom_table_cols

@inline _row_by_idx(at::_AtomTable{T}, idx::Int) where T = _AtomTableRow{T}(getfield(at, :_idx_map)[idx], at)

@inline function Base.getproperty(atr::_AtomTableRow{T}, nm::Symbol) where T
    nm === :idx            && return _getproperty(atr, :idx)::Int
    nm === :number         && return _getproperty(atr, :number)::Int
    nm === :element        && return _getproperty(atr, :element)::ElementType
    nm === :name           && return _getproperty(atr, :name)::String
    nm === :atom_type      && return _getproperty(atr, :atom_type)::String
    nm === :r              && return _getproperty(atr, :r)::Vector3{T}
    nm === :v              && return _getproperty(atr, :v)::Vector3{T}
    nm === :F              && return _getproperty(atr, :F)::Vector3{T}
    nm === :properties     && return _getproperty(atr, :properties)::Properties
    nm === :flags          && return _getproperty(atr, :flags)::Flags
    nm === :frame_id       && return _getproperty(atr, :frame_id)::Int
    nm === :molecule_idx   && return _getproperty(atr, :molecule_idx)::MaybeInt
    nm === :chain_idx      && return _getproperty(atr, :chain_idx)::MaybeInt
    nm === :fragment_idx   && return _getproperty(atr, :fragment_idx)::MaybeInt
    nm === :nucleotide_idx && return _getproperty(atr, :nucleotide_idx)::MaybeInt
    nm === :residue_idx    && return _getproperty(atr, :residue_idx)::MaybeInt
    getindex(getfield(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

@inline function _getproperty(atr::_AtomTableRow{T}, nm::Symbol) where T
    getindex(getproperty(getfield(atr, :_tab), nm), getfield(atr, :_row))
end

@inline Base.setproperty!(atr::_AtomTableRow, nm::Symbol, val) = setindex!(getproperty(getfield(atr, :_tab), nm), val, getfield(atr, :_row))

@inline Base.eltype(::_AtomTable{T}) where T = _AtomTableRow{T}
@inline Base.iterate(at::_AtomTable, st=1) = st > length(at) ? nothing : (_AtomTableRow(st, at), st + 1)
