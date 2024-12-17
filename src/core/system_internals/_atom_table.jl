const _atom_table_cols_main = (:idx, :number, :element, :name, :atom_type, :r, :v, :F, :formal_charge, :charge, :radius)
const _atom_table_cols_extra = (:properties, :flags, :frame_id, :molecule_idx, :chain_idx, :fragment_idx)
const _atom_table_cols = (_atom_table_cols_main..., _atom_table_cols_extra...)
const _atom_table_cols_set = Set(_atom_table_cols)

@auto_hash_equals struct _AtomTable{T <: Real} <: _AbstractSystemComponentTable
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

    # internals
    _idx_map::_IdxMap

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
            _IdxMap()
        )
    end
end

@inline function _hascolumn(::Union{_AtomTable, Type{<: _AtomTable}}, nm::Symbol)
    nm in _atom_table_cols_set
end

@inline function Tables.columnnames(::_AtomTable)
    _atom_table_cols_main
end

@inline function Tables.schema(::_AtomTable{T}) where T
    Tables.Schema(
        _atom_table_cols_main,
        (Int, Int, ElementType, String, String, Vector3{T}, Vector3{T}, Vector3{T}, Int, T, T)
    )
end

@inline function Base.propertynames(::_AtomTable)
    _atom_table_cols
end

function Base.push!(
    at::_AtomTable{T},
    idx::Int,
    number::Int,
    element::ElementType;
    name::AbstractString = "",
    atom_type::AbstractString = "",
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
    fragment_idx::MaybeInt = nothing
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
            fragment_idx = a.fragment_idx
        )
    end
    at
end
@inline Tables.materializer(::Type{_AtomTable{T}}) where T = itr -> _atom_table(T, itr)
