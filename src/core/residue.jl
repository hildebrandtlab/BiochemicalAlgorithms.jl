export
    Residue,
    ResidueTable,
    residue_by_idx,
    residues,
    residues_df,
    eachresidue,
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

function Base.push!(rt::ResidueTable, t::ResidueTuple, molecule_id::Int, chain_id::Int)
    sys = getfield(rt, :_sys)
    push!(sys._residues, t, molecule_id, chain_id)
    push!(getfield(rt, :_idx), sys._curr_idx)
    rt
end

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

function Residue(
    chain::Chain{T},
    number::Int,
    type::AminoAcid,
    properties::Properties = Properties(),
    flags::Flags = Flags()
) where T
    sys = parent(chain)
    idx = _next_idx(sys)
    push!(sys._residues, ResidueTuple(number, type;
            idx = idx,
            properties = properties,
            flags = flags
        ), chain._row.molecule_id, chain.idx
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
@inline parent_molecule(res::Residue) = molecule_by_idx(parent(res), res._row.molecule_id)
@inline parent_chain(res::Residue) = chain_by_idx(parent(res), res._row.chain_id)

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
    residues_df(::Chain)
    residues_df(::Molecule)
    residues_df(::Protein)
    residues_df(::System)

Returns a `DataFrame{T}` containing all residues of the given atom container.
"""
@inline function residues_df(sys::System; kwargs...)
    DataFrame(residues(sys; kwargs...))
end

"""
    eachresidue(::Chain)
    eachresidue(::Molecule)
    eachresidue(::Protein)
    eachresidue(::System)

Returns a `Residue{T}` generator for all residues of the given atom container.
"""
@inline function eachresidue(sys::System{T}; kwargs...) where T
    (res for res in residues(sys; kwargs...))
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
@inline residues_df(mol::Molecule; kwargs...) = residues_df(parent(mol); molecule_id = mol.idx, kwargs...)
@inline eachresidue(mol::Molecule; kwargs...) = eachresidue(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nresidues(mol::Molecule; kwargs...) = nresidues(parent(mol); molecule_id = mol.idx, kwargs...)

#=
    Chain residues
=#
@inline residues(chain::Chain; kwargs...) = residues(parent(chain); chain_id = chain.idx, kwargs...)
@inline residues_df(chain::Chain; kwargs...) = residues_df(parent(chain); chain_id = chain.idx, kwargs...)
@inline eachresidue(chain::Chain; kwargs...) = eachresidue(parent(chain); chain_id = chain.idx, kwargs...)
@inline nresidues(chain::Chain; kwargs...) = nresidues(parent(chain); chain_id = chain.idx, kwargs...)

"""
    push!(::Chain, frag::FragmentTuple)

Creates a new fragment in the given chain, based on the given tuple. The new residue is automatically
assigned a new `idx`.
"""
@inline function Base.push!(chain::Chain, res::ResidueTuple)
    sys = parent(chain)
    push!(sys._residues,
        (; res..., idx = _next_idx(sys)),
        chain._row.molecule_id,
        chain.idx
    )
    chain
end

#=
    Residue atoms
=#
@inline atoms(res::Residue; kwargs...) = atoms(parent(res); residue_id = res.idx, kwargs...)
@inline atoms_df(res::Residue; kwargs...) = atoms_df(parent(res); residue_id = res.idx, kwargs...)
@inline eachatom(res::Residue; kwargs...) = eachatom(parent(res); residue_id = res.idx, kwargs...)
@inline natoms(res::Residue; kwargs...) = natoms(parent(res); residue_id = res.idx, kwargs...)

@inline function Base.push!(res::Residue{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(parent(res), atom; molecule_id = res._row.molecule_id, chain_id = res._row.chain_id,
        residue_id = res.idx, kwargs...)
    res
end

#=
    Residue bonds
=#
@inline bonds(res::Residue; kwargs...) = bonds(parent(res); residue_id = res.idx, kwargs...)
@inline bonds_df(res::Residue; kwargs...) = bonds_df(parent(res); residue_id = res.idx, kwargs...)
@inline eachbond(res::Residue; kwargs...) = eachbond(parent(res); residue_id = res.idx, kwargs...)
@inline nbonds(res::Residue; kwargs...) = nbonds(parent(res); residue_id = res.idx, kwargs...)
