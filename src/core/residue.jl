export
    Residue,
    ResidueTable,
    residue_by_idx,
    residues,
    nresidues,
    parent_residue

@auto_hash_equals struct ResidueTable{T} <: Tables.AbstractColumns
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _residues(rt::ResidueTable) = getproperty(getfield(rt, :_sys), :_residues)

@inline Tables.istable(::Type{<: ResidueTable}) = true
@inline Tables.columnaccess(::Type{<: ResidueTable}) = true
@inline Tables.columns(rt::ResidueTable) = rt

@inline function Tables.getcolumn(rt::ResidueTable, nm::Symbol)
    col = Tables.getcolumn(_residues(rt), nm)
    RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _residues(rt)._idx_map[idx], getfield(rt, :_idx))
    )
end

@inline function Base.getproperty(rt::ResidueTable, nm::Symbol)
    hasfield(typeof(rt), nm) && return getfield(rt, nm)
    Tables.getcolumn(rt, nm)
end

@inline Tables.getcolumn(rt::ResidueTable, i::Int) = Tables.getcolumn(rt, Tables.columnnames(rt)[i])
@inline Tables.columnnames(rt::ResidueTable) = Tables.columnnames(_residues(rt))
@inline Tables.schema(rt::ResidueTable) = Tables.schema(_residues(rt))

@inline Base.size(rt::ResidueTable) = (length(getfield(rt, :_idx)), length(_residue_table_cols))
@inline Base.size(rt::ResidueTable, dim) = size(rt)[dim]
@inline Base.length(rt::ResidueTable) = size(rt, 1)

@inline function _filter_residues(f::Function, sys::System{T}) where T
    ResidueTable{T}(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._residues),
        :idx
    )))
end

@inline function Base.filter(f::Function, rt::ResidueTable)
    ResidueTable(getfield(rt, :_sys), collect(Int, _filter_select(
        TableOperations.filter(f, rt),
        :idx
    )))
end

@inline function Base.iterate(rt::ResidueTable, st = 1)
    st > length(rt) ?
        nothing :
        (residue_by_idx(getfield(rt, :_sys), getfield(rt, :_idx)[st]), st + 1)
end
@inline Base.eltype(::ResidueTable{T}) where T = Residue{T}
@inline Base.getindex(rt::ResidueTable{T}, i::Int) where T = residue_by_idx(getfield(rt, :_sys), getfield(rt, :_idx)[i])
@inline Base.keys(rt::ResidueTable) = LinearIndices((length(rt),))

"""
    $(TYPEDEF)

Mutable representation of an individual residue in a system.

# Fields
 - `idx::Int`
 - `number::Int`
 - `type::AminoAcid`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Residue(
    chain::Chain{T},
    number::Int,
    type::AminoAcid,
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
    push!(
        sys._residues,
        _Residue(number, type; idx = idx, kwargs...),
        chain.molecule_id,
        chain.idx
    )
    residue_by_idx(sys, idx)
end

@inline Tables.rows(rt::ResidueTable) = rt
@inline Tables.getcolumn(res::Residue, name::Symbol) = Tables.getcolumn(getfield(res, :_row), name)

@inline function Base.getproperty(res::Residue, name::Symbol)
    hasfield(typeof(res), name) && return getfield(res, name)
    getproperty(getfield(res, :_row), name)
end

@inline function Base.setproperty!(res::Residue, name::Symbol, val)
    hasfield(typeof(res), name) && return setfield!(res, name, val)
    setproperty!(getfield(res, :_row), name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", res::Residue) = show(io, getfield(res, :_row))
@inline Base.show(io::IO, res::Residue) = show(io, getfield(res, :_row))

@inline Base.parent(res::Residue) = res._sys
@inline parent_system(res::Residue) = parent(res)
@inline parent_molecule(res::Residue) = molecule_by_idx(parent(res), res.molecule_id)
@inline parent_chain(res::Residue) = chain_by_idx(parent(res), res.chain_id)

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
    residues(::Protein)
    residues(::System)

Returns a `ResidueTable{T}` containing all residues of the given atom container.

# Supported keyword arguments
 - `molecule_id::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to residues belonging to the molecule with the given ID.
- `chain_id::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to residues belonging to the chain with the given ID.
"""
function residues(sys::System{T};
    molecule_id::MaybeInt = nothing,
    chain_id::MaybeInt = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return ResidueTable{T}(sys, sys._residues.idx)
    _filter_residues(res ->
        (isnothing(molecule_id) || res.molecule_id == something(molecule_id)) &&
        (isnothing(chain_id)    || res.chain_id    == something(chain_id)),
        sys
    )
end

"""
    nresidues(::Chain)
    nresidues(::Molecule)
    nresidues(::Protein)
    nresidues(::System)

Returns the number of residues in the given atom container.
"""
@inline function nresidues(sys::System; kwargs...)
    length(residues(sys; kwargs...))
end

#=
    Molecule residues
=#
@inline residues(mol::Molecule; kwargs...) = residues(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nresidues(mol::Molecule; kwargs...) = nresidues(parent(mol); molecule_id = mol.idx, kwargs...)

#=
    Chain residues
=#
@inline residues(chain::Chain; kwargs...) = residues(parent(chain); chain_id = chain.idx, kwargs...)
@inline nresidues(chain::Chain; kwargs...) = nresidues(parent(chain); chain_id = chain.idx, kwargs...)

"""
    push!(::Chain{T}, res::Residue{T})

Creates a copy of the given residue in the given chain. The new residue is automatically assigned a
new `idx`.
"""
@inline function Base.push!(chain::Chain{T}, res::Residue{T}) where T
    Residue(parent(chain), res.number, res.type;
        properties = res.properties,
        flags = res.flags
    )
    chain
end

#=
    Residue atoms
=#
@inline atoms(res::Residue; kwargs...) = atoms(parent(res); residue_id = res.idx, kwargs...)
@inline natoms(res::Residue; kwargs...) = natoms(parent(res); residue_id = res.idx, kwargs...)

@inline function Atom(res::Residue, number::Int, element::ElementType; kwargs...)
    Atom(parent(res), number, element;
        molecule_id = res.molecule_id,
        chain_id = res.chain_id,
        residue_id = res.idx,
        kwargs...
    )
end

@inline function Base.push!(res::Residue{T}, atom::Atom{T}; kwargs...) where T
    push!(parent(res), atom;
        molecule_id = res.molecule_id,
        chain_id = res.chain_id,
        residue_id = res.idx,
        kwargs...
    )
    res
end

#=
    Residue bonds
=#
@inline bonds(res::Residue; kwargs...) = bonds(parent(res); residue_id = res.idx, kwargs...)
@inline nbonds(res::Residue; kwargs...) = nbonds(parent(res); residue_id = res.idx, kwargs...)
