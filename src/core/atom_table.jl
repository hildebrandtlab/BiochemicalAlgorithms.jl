export
    AtomTable

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

Tables.getcolumn(at::AtomTable, nm::Symbol) = nm in _atom_table_cols_set ? getfield(at, nm) : error("type AtomTable has no column $nm")
Tables.getcolumn(at::AtomTable, i::Int) = getfield(at, Tables.columnnames(at)[i])
Tables.columnnames(::AtomTable) = _atom_table_cols
Tables.schema(::AtomTable{T}) where T = Tables.Schema(fieldnames(AtomTuple{T}), fieldtypes(AtomTuple{T}))

_getcolumn(at::AtomTable, nm::Symbol) = nm in _atom_table_cols_priv || nm in _atom_table_cols_set ? getfield(at, nm) : error("type AtomTable has no column $nm")

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

