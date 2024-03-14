export
    Residue,
    ResidueTable,
    residue_by_idx,
    residues,
    nresidues,
    parent_residue

"""
    $(TYPEDEF)

Tables.jl-compatible representation of system residues (or a subset thereof). Residue tables can be
generated using [`residues`](@ref) or filtered from other residue tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `number::AbstractVector{Int}`
 - `type::AbstractVector{AminoAcid}`

# Private columns
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
 - `molecule_idx::AbstractVector{Int}`
 - `chain_idx::AbstractVector{Int}`
"""
@auto_hash_equals struct ResidueTable{T} <: AbstractSystemComponentTable{T}
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _residues(rt::ResidueTable) = getfield(getfield(rt, :_sys), :_residues)
@inline _hascolumn(::ResidueTable, nm::Symbol) = nm in _residue_table_cols_priv || nm in _residue_table_cols_set

@inline function Tables.getcolumn(rt::ResidueTable, nm::Symbol)
    col = Tables.getcolumn(_residues(rt), nm)
    _RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _residues(rt)._idx_map[idx], rt._idx)
    )
end

@inline Tables.columnnames(rt::ResidueTable) = Tables.columnnames(_residues(rt))
@inline Tables.schema(rt::ResidueTable) = Tables.schema(_residues(rt))

@inline function _filter_residues(f::Function, sys::System)
    ResidueTable(sys, _filter_idx(f, sys._residues))
end

@inline function Base.filter(f::Function, rt::ResidueTable)
    ResidueTable(rt._sys, _filter_idx(f, rt))
end

@inline function Base.iterate(rt::ResidueTable, st = 1)
    st > length(rt) ?
        nothing :
        (residue_by_idx(rt._sys, rt._idx[st]), st + 1)
end

@inline Base.eltype(::ResidueTable{T}) where T = Residue{T}
@inline Base.size(rt::ResidueTable) = (length(rt._idx), length(Tables.columnnames(rt)))
@inline Base.getindex(rt::ResidueTable, i::Int) = residue_by_idx(rt._sys, rt._idx[i])
@inline Base.getindex(rt::ResidueTable, ::Colon) = rt

@inline function Base.getindex(rt::ResidueTable, I)
    ResidueTable(rt._sys, collect(Int, map(i -> rt._idx[i], I)))
end

"""
    $(TYPEDEF)

Mutable representation of an individual residue in a system.

# Public fields
 - `idx::Int`
 - `number::Int`
 - `type::AminoAcid`

# Private fields
 - `properties::Properties`
 - `flags::Flags`
 - `molecule_idx::Int`
 - `chain_idx::Int`

# Constructors
```julia
Residue(
    chain::Chain{T},
    number::Int,
    type::AminoAcid;
    # keyword arguments
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Residue{T}` in the given chain.
"""
@auto_hash_equals struct Residue{T} <: AbstractAtomContainer{T}
    _sys::System{T}
    _row::_ResidueTableRow
end

@inline function Residue(
    chain::Chain{T},
    number::Int,
    type::AminoAcid;
    kwargs...
) where T
    sys = parent(chain)
    idx = _next_idx(sys)
    push!(sys._residues, idx, number, type, chain.molecule_idx, chain.idx; kwargs...)
    residue_by_idx(sys, idx)
end

@inline function Base.getproperty(res::Residue, name::Symbol)
    hasfield(typeof(res), name) && return getfield(res, name)
    getproperty(getfield(res, :_row), name)
end

@inline function Base.setproperty!(res::Residue, name::Symbol, val)
    hasfield(typeof(res), name) && return setfield!(res, name, val)
    setproperty!(getfield(res, :_row), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", res::Residue) = show(io, res)
@inline function Base.show(io::IO, res::Residue)
    print(io, "$(typeof(res)): ")
    show(io, NamedTuple(res._row))
end

@inline Base.parent(res::Residue) = res._sys
@inline parent_system(res::Residue) = parent(res)
@inline parent_molecule(res::Residue) = molecule_by_idx(parent(res), res.molecule_idx)
@inline parent_chain(res::Residue) = chain_by_idx(parent(res), res.chain_idx)

@doc raw"""
    parent_residue(::Atom)

Returns the `Residue{T}` containing the given atom. Returns `nothing` if no such residue exists.
""" parent_residue

"""
    $(TYPEDSIGNATURES)

Returns the `Residue{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
residue exists.
"""
@inline function residue_by_idx(sys::System{T}, idx::Int) where T
    Residue{T}(sys, _row_by_idx(sys._residues, idx))
end

"""
    residues(::Chain)
    residues(::Molecule)
    residues(::System)

Returns a `ResidueTable{T}` containing all residues of the given atom container.

# Supported keyword arguments
 - `molecule_idx::MaybeInt = nothing`
 - `chain_idx::MaybeInt = nothing`
All keyword arguments limit the results to residues matching the given IDs. Keyword arguments set to
`nothing` are ignored.
"""
function residues(sys::System{T};
    molecule_idx::MaybeInt = nothing,
    chain_idx::MaybeInt = nothing
) where T
    isnothing(molecule_idx) && isnothing(chain_idx) && return ResidueTable{T}(sys, sys._residues.idx)
    _filter_residues(res ->
        (isnothing(molecule_idx) || res.molecule_idx == something(molecule_idx)) &&
        (isnothing(chain_idx)    || res.chain_idx    == something(chain_idx)),
        sys
    )
end

"""
    nresidues(::Chain)
    nresidues(::Molecule)
    nresidues(::System)

Returns the number of residues in the given atom container.

# Supported keyword arguments
See [`residues`](@ref)
"""
@inline function nresidues(sys::System; kwargs...)
    length(residues(sys; kwargs...))
end

#=
    Molecule residues
=#
@inline residues(mol::Molecule; kwargs...) = residues(parent(mol); molecule_idx = mol.idx, kwargs...)
@inline nresidues(mol::Molecule; kwargs...) = nresidues(parent(mol); molecule_idx = mol.idx, kwargs...)

#=
    Chain residues
=#
@inline residues(chain::Chain; kwargs...) = residues(parent(chain); chain_idx = chain.idx, kwargs...)
@inline nresidues(chain::Chain; kwargs...) = nresidues(parent(chain); chain_idx = chain.idx, kwargs...)

"""
    push!(::Chain{T}, ::Residue{T})

Creates a copy of the given residue in the given chain. The new residue is automatically assigned a
new `idx`.
"""
@inline function Base.push!(chain::Chain{T}, res::Residue{T}) where T
    Residue(chain, res.number, res.type;
        properties = res.properties,
        flags = res.flags
    )
    chain
end

#=
    Residue atoms
=#
@inline atoms(res::Residue; kwargs...) = atoms(parent(res); residue_idx = res.idx, kwargs...)
@inline natoms(res::Residue; kwargs...) = natoms(parent(res); residue_idx = res.idx, kwargs...)

@inline function Atom(res::Residue, number::Int, element::ElementType; kwargs...)
    Atom(parent(res), number, element;
        molecule_idx = res.molecule_idx,
        chain_idx = res.chain_idx,
        residue_idx = res.idx,
        kwargs...
    )
end

@inline function Base.push!(res::Residue{T}, atom::Atom{T}; kwargs...) where T
    push!(parent(res), atom;
        molecule_idx = res.molecule_idx,
        chain_idx = res.chain_idx,
        residue_idx = res.idx,
        kwargs...
    )
    res
end

#=
    Residue bonds
=#
@inline bonds(res::Residue; kwargs...) = bonds(parent(res); residue_idx = res.idx, kwargs...)
@inline nbonds(res::Residue; kwargs...) = nbonds(parent(res); residue_idx = res.idx, kwargs...)
