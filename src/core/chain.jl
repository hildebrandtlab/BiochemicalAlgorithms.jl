using AutoHashEquals
export Chain, chain_by_idx, chains, chains_df, eachchain, nchains, parent_chain

"""
    $(TYPEDEF)

Mutable representation of an individual chain in a system.

# Fields
- `idx::Int`
- `name::String`
- `properties::Properties`

# Constructors
```julia
Chain(mol::Molecule{T}, name::String = "", properties::Properties = Properties())
```
Creates a new `Chain{T}` in the given molecule.
"""
@auto_hash_equals struct Chain{T} <: AbstractAtomContainer{T}
    _sys::System{T}
    _row::DataFrameRow
end

function Chain(
    mol::Molecule{T},
    name::String = "",
    properties::Properties = Properties()
) where T
    sys = parent(mol)
    idx = _next_idx(sys)
    push!(sys._chains, (idx, name, properties, mol.idx))
    chain_by_idx(sys, idx)
end

function Base.getproperty(chain::Chain, name::Symbol)
    in(name, fieldnames(ChainTuple)) && return getproperty(getfield(chain, :_row), name)
    getfield(chain, name)
end

function Base.setproperty!(chain::Chain, name::Symbol, val)
    in(name, fieldnames(ChainTuple)) && return setproperty!(getfield(chain, :_row), name, val)
    setfield!(chain, name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", chain::Chain) = show(io, getfield(chain, :_row))
@inline Base.show(io::IO, chain::Chain) = show(io, getfield(chain, :_row))

@inline Base.parent(chain::Chain) = chain._sys
@inline parent_system(chain::Chain) = parent(chain)
@inline parent_molecule(chain::Chain) = _molecule_by_idx(parent(chain), chain._row.molecule_id)

@doc raw"""
    parent_chain(::Atom)
    parent_chain(::Fragment)
    parent_chain(::Nucleotide)
    parent_chain(::Residue)

Returns the `Chain{T}` containing the given object. Returns `nothing` if no such chain exists.
""" parent_chain

"""
    $(TYPEDSIGNATURES)

Returns the `Chain{T}` associated with the given `idx` in `sys`. Returns `nothing` if no such chain
exists.
"""
@inline function chain_by_idx(sys::System{T}, idx::Int) where T
    rn = _row_by_idx(sys._chains, idx)
    isnothing(rn) ? nothing : Chain{T}(sys, DataFrameRow(sys._chains, rn, :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's chains matching the given criteria. Fields
given as `nothing` are ignored. The returned `DataFrame` contains all public and private chain fields.
"""
function _chains(sys::System; molecule_id::Union{Nothing, Int} = nothing)
    isnothing(molecule_id) && return sys._chains

    get(
        groupby(sys._chains, :molecule_id),
        (molecule_id = molecule_id,),
        DataFrame(_SystemChainTuple[])
    )
end

"""
    chains(::Molecule)
    chains(::Protein)
    chains(::System)

Returns a `Vector{Chain{T}}` containing all chains of the given atom container.
"""
@inline function chains(sys::System; kwargs...)
    collect(eachchain(sys; kwargs...))
end

"""
    chains_df(::Molecule)
    chains_df(::Protein)
    chains_df(::System)

Returns a `SubDataFrame` containing all chains of the given atom container.
"""
@inline function chains_df(sys::System; kwargs...)
    view(_chains(sys; kwargs...), :, 1:length(fieldnames(ChainTuple)))
end

"""
    eachchain(::Molecule)
    eachchain(::Protein)
    eachchain(::System)

Returns a `Chain{T}` generator for all chains of the given atom container.
"""
@inline function eachchain(sys::System{T}; kwargs...) where T
    (Chain{T}(sys, row) for row in eachrow(_chains(sys; kwargs...)))
end

"""
    nchains(::Molecule)
    nchains(::Protein)
    nchains(::System)

Returns the number of chains in the given atom container.
"""
@inline function nchains(sys::System; kwargs...)
    nrow(_chains(sys; kwargs...))
end

#=
    Molecule chains
=#
@inline _chains(mol::Molecule) = _chains(parent(mol); molecule_id = mol.idx)
@inline chains(mol::Molecule) = chains(parent(mol), molecule_id = mol.idx)
@inline chains_df(mol::Molecule) = chains_df(parent(mol), molecule_id = mol.idx)
@inline eachchain(mol::Molecule) = eachchain(parent(mol), molecule_id = mol.idx)
@inline nchains(mol::Molecule) = nchains(parent(mol), molecule_id = mol.idx)

"""
    push!(::Molecule, chain::ChainTuple)
    push!(::Protein, chain::ChainTuple)

Creates a new chain in the given molecule, based on the given tuple. The new chain is automatically
assigned a new `idx`.
"""
@inline function Base.push!(mol::Molecule, chain::ChainTuple)
    push!(sys._chains, (_with_idx(chain, _next_idx(sys))..., mol.idx))
    mol
end

#=
    Chain atoms
=#
@inline _atoms(chain::Chain; kwargs...) = _atoms(parent(chain); chain_id = chain.idx, kwargs...)
@inline atoms(chain::Chain; kwargs...) = atoms(parent(chain); chain_id = chain.idx, kwargs...)
@inline atoms_df(chain::Chain; kwargs...) = atoms_df(parent(chain); chain_id = chain.idx, kwargs...)
@inline eachatom(chain::Chain; kwargs...) = eachatom(parent(chain); chain_id = chain.idx, kwargs...)
@inline natoms(chain::Chain; kwargs...) = natoms(parent(chain); chain_id = chain.idx, kwargs...)

#=
    Chain bonds
=#
@inline _bonds(chain::Chain; kwargs...) = _bonds(parent(chain); chain_id = chain.idx, kwargs...)
@inline bonds(chain::Chain; kwargs...) = bonds(parent(chain); chain_id = chain.idx, kwargs...)
@inline bonds_df(chain::Chain; kwargs...) = bonds_df(parent(chain); chain_id = chain.idx, kwargs...)
@inline eachbond(chain::Chain; kwargs...) = eachbond(parent(chain); chain_id = chain.idx, kwargs...)
@inline nbonds(chain::Chain; kwargs...) = nbonds(parent(chain); chain_id = chain.idx, kwargs...)
