using AutoHashEquals
export Residue, residues, residues_df, eachresidue, nresidues, parent_residue

"""
    $(TYPEDEF)

Mutable representation of an individual residue in a system.

# Fields
- `idx::Int`
- `number::Int`
- `type::AminoAcid`
- `properties::Properties`

# Constructors
```julia
Residue(
    chain::Chain{T},
    number::Int,
    type::AminoAcid,
    properties::Properties = Properties()
)
```
Creates a new `Residue{T}` in the given chain.
"""
@auto_hash_equals struct Residue{T}
    sys::System{T}
    row::DataFrameRow
end

function Residue(
    chain::Chain{T},
    number::Int,
    type::AminoAcid,
    properties::Properties = Properties()
) where T
    sys = chain.sys
    idx = _next_idx(sys)
    push!(sys.residues, (idx, number, type, properties, chain.row.molecule_id, chain.idx))
    _residue_by_idx(sys, idx)
end

function Base.getproperty(res::Residue, name::Symbol)
    in(name, fieldnames(ResidueTuple)) && return getproperty(getfield(res, :row), name)
    getfield(res, name)
end

function Base.setproperty!(res::Residue, name::Symbol, val)
    in(name, fieldnames(ResidueTuple)) && return setproperty!(getfield(res, :row), name, val)
    setfield!(res, name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", res::Residue) = show(io, getfield(res, :row))
@inline Base.show(io::IO, res::Residue) = show(io, getfield(res, :row))

@inline Base.parent(res::Residue) = res.sys
@inline parent_system(res::Residue) = parent(res)
@inline parent_molecule(res::Residue) = _molecule_by_idx(res.sys, res.row.molecule_id)
@inline parent_chain(res::Residue) = _chain_by_idx(res.sys, res.row.chain_id)

@doc raw"""
    parent_residue(::Atom)

Returns the `Residue{T}` containing the given atom. Returns `nothing` if no such residue exists.
""" parent_residue

"""
    $(TYPEDSIGNATURES)

Returns the `Residue{T}` associated with the given `idx` in `sys`.
"""
@inline function _residue_by_idx(sys::System{T}, idx::Int) where T
    Residue{T}(sys, DataFrameRow(sys.residues, findfirst(sys.residues.idx .== idx), :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's residues matching the given criteria. Fields given
as `nothing` are ignored. The returned `DataFrame` contains all public and private residue fields.
"""
function _residues(sys::System{T};
    molecule_id::Union{Nothing, Int} = nothing,
    chain_id::Union{Nothing, Int} = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return sys.residues

    cols = Tuple{Symbol, Int}[]
    isnothing(molecule_id) || push!(cols, (:molecule_id, molecule_id))
    isnothing(chain_id)    || push!(cols, (:chain_id, chain_id))

    get(
        groupby(sys.residues, getindex.(cols, 1)),
        ntuple(i -> cols[i][2], length(cols)),
        DataFrame(_SystemResidueTuple[])
    )
end

"""
    residues(::Chain)
    residues(::Molecule)
    residues(::Protein)
    residues(::System)

Returns a `Vector{Residue{T}}` containing all residues of the given atom container.
"""
@inline function residues(sys::System; kwargs...)
    collect(eachresidue(sys; kwargs...))
end

"""
    residues_df(::Chain)
    residues_df(::Molecule)
    residues_df(::Protein)
    residues_df(::System)

Returns a `SystemDataFrame{T}` containing all residues of the given atom container.
"""
@inline function residues_df(sys::System; kwargs...)
    SystemDataFrame(sys, view(_residues(sys; kwargs...), :, 1:length(fieldnames(ResidueTuple))))
end

"""
    eachresidue(::Chain)
    eachresidue(::Molecule)
    eachresidue(::Protein)
    eachresidue(::System)

Returns a `Residue{T}` generator for all residues of the given atom container.
"""
@inline function eachresidue(sys::System{T}; kwargs...) where T
    (Residue{T}(sys, row) for row in eachrow(_residues(sys; kwargs...)))
end

"""
    nresidues(::Chain)
    nresidues(::Molecule)
    nresidues(::Protein)
    nresidues(::System)

Returns the number of residues in the given atom container.
"""
@inline function nresidues(sys::System; kwargs...)
    nrow(_residues(sys; kwargs...))
end

#=
    Molecule residues
=#
@inline _residues(mol::Molecule; kwargs...) = _residues(mol.sys; molecule_id = mol.idx, kwargs...)
@inline residues(mol::Molecule; kwargs...) = residues(mol.sys; molecule_id = mol.idx, kwargs...)
@inline residues_df(mol::Molecule; kwargs...) = residues_df(mol.sys; molecule_id = mol.idx, kwargs...)
@inline eachresidue(mol::Molecule; kwargs...) = eachresidue(mol.sys; molecule_id = mol.idx, kwargs...)
@inline nresidues(mol::Molecule; kwargs...) = nresidues(mol.sys; molecule_id = mol.idx, kwargs...)

#=
    Chain residues
=#
@inline _residues(chain::Chain; kwargs...) = _residues(chain.sys; chain_id = chain.idx, kwargs...)
@inline residues(chain::Chain; kwargs...) = residues(chain.sys; chain_id = chain.idx, kwargs...)
@inline residues_df(chain::Chain; kwargs...) = residues_df(chain.sys; chain_id = chain.idx, kwargs...)
@inline eachresidue(chain::Chain; kwargs...) = eachresidue(chain.sys; chain_id = chain.idx, kwargs...)
@inline nresidues(chain::Chain; kwargs...) = nresidues(chain.sys; chain_id = chain.idx, kwargs...)

"""
    push!(::Chain, frag::FragmentTuple)

Creates a new fragment in the given chain, based on the given tuple. The new residue is automatically
assigned a new `idx`.
"""
@inline function Base.push!(chain::Chain, res::ResidueTuple)
    push!(chain.sys.residues, (_with_idx(res, _next_idx(chain.sys))..., chain.idx))
    chain
end

#=
    Residue atoms
=#
@inline _atoms(res::Residue; kwargs...) = _atoms(res.sys; residue_id = res.idx, kwargs...)
@inline atoms(res::Residue; kwargs...) = atoms(res.sys; residue_id = res.idx, kwargs...)
@inline atoms_df(res::Residue; kwargs...) = atoms_df(res.sys; residue_id = res.idx, kwargs...)
@inline eachatom(res::Residue; kwargs...) = eachatom(res.sys; residue_id = res.idx, kwargs...)
@inline natoms(res::Residue; kwargs...) = natoms(res.sys; residue_id = res.idx, kwargs...)

@inline function Base.push!(res::Residue{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(res.sys, atom; molecule_id = res.row.molecule_id, chain_id = res.row.chain_id,
        residue_id = res.idx, kwargs...)
    res
end

#=
    Residue bonds
=#
@inline _bonds(res::Residue; kwargs...) = _bonds(res.sys; residue_id = res.idx, kwargs...)
@inline bonds(res::Residue; kwargs...) = bonds(res.sys; residue_id = res.idx, kwargs...)
@inline bonds_df(res::Residue; kwargs...) = bonds_df(res.sys; residue_id = res.idx, kwargs...)
@inline eachbond(res::Residue; kwargs...) = eachbond(res.sys; residue_id = res.idx, kwargs...)
@inline nbonds(res::Residue; kwargs...) = nbonds(res.sys; residue_id = res.idx, kwargs...)

@inline function Base.push!(res::Residue, bond::Bond)
    push!(res.sys, bond)
    res
end
