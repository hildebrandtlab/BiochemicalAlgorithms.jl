export
    Chain,
    ChainTable,
    chain_by_idx,
    chains,
    nchains,
    parent_chain

"""
    $(TYPEDEF)

Tables.jl-compatible representation of system chains (or a subset thereof). Chain tables can be
generated using [`chains`](@ref) or filtered from other chain tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `name::AbstractVector{String}`
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`

# Private columns
 - `molecule_idx::AbstractVector{Int}`
"""
@auto_hash_equals struct ChainTable{T} <: AbstractSystemComponentTable{T}
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _chains(ct::ChainTable) = getfield(getfield(ct, :_sys), :_chains)

@inline function Tables.getcolumn(ct::ChainTable, nm::Symbol)
    col = Tables.getcolumn(_chains(ct), nm)
    _RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _chains(ct)._idx_map[idx], ct._idx)
    )
end

@inline Tables.columnnames(ct::ChainTable) = Tables.columnnames(_chains(ct))
@inline Tables.schema(ct::ChainTable) = Tables.schema(_chains(ct))

@inline function Base.getproperty(ct::ChainTable, nm::Symbol)
    hasfield(typeof(ct), nm) && return getfield(ct, nm)
    Tables.getcolumn(ct, nm)
end

@inline function Base.setproperty!(ct::ChainTable, nm::Symbol, val)
    if nm in _chain_table_cols_priv || nm in _chain_table_cols_set
        error("ChainTable columns cannot be set directly! Did you mean to use broadcast assignment (.=)?")
    end
    if !hasfield(typeof(ct), nm)
        error("type ChainTable has no field $nm")
    end
    setfield!(ct, nm, val)
end

@inline function _filter_chains(f::Function, sys::System)
    ChainTable(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._chains),
        :idx
    )))
end

@inline function Base.filter(f::Function, ct::ChainTable)
    ChainTable(ct._sys, collect(Int, _filter_select(
        TableOperations.filter(f, ct),
        :idx
    )))
end

@inline function Base.iterate(ct::ChainTable, st = 1)
    st > length(ct) ?
        nothing :
        (chain_by_idx(ct._sys, ct._idx[st]), st + 1)
end

@inline Base.eltype(::ChainTable{T}) where T = Chain{T}
@inline Base.size(ct::ChainTable) = (length(ct._idx), length(Tables.columnnames(ct)))
@inline Base.getindex(ct::ChainTable, i::Int) = chain_by_idx(ct._sys, ct._idx[i])
@inline Base.getindex(ct::ChainTable, ::Colon) = ct

@inline function Base.getindex(ct::ChainTable, I)
    ChainTable(ct._sys, collect(Int, map(i -> ct._idx[i], I)))
end

"""
    $(TYPEDEF)

Mutable representation of an individual chain in a system.

# Public fields
 - `idx::Int`
 - `name::String`
 - `properties::Properties`
 - `flags::Flags`

# Private fields
 - `molecule_idx::Int`

# Constructors
```julia
Chain(
    mol::Molecule{T};
    # keyword arguments
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Chain{T}` in the given molecule.
"""
@auto_hash_equals struct Chain{T} <: AbstractAtomContainer{T}
    _sys::System{T}
    _row::_ChainTableRow
end

@inline function Chain(
    mol::Molecule;
    kwargs...
)
    sys = parent(mol)
    idx = _next_idx(sys)
    push!(sys._chains, idx, mol.idx; kwargs...)
    chain_by_idx(sys, idx)
end

@inline function Base.getproperty(chain::Chain, name::Symbol)
    hasfield(typeof(chain), name) && return getfield(chain, name)
    getproperty(getfield(chain, :_row), name)
end

@inline function Base.setproperty!(chain::Chain, name::Symbol, val)
    hasfield(typeof(chain), name) && return setfield!(chain, name, val)
    setproperty!(getfield(chain, :_row), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", chain::Chain) = show(io, getfield(chain, :_row))
@inline Base.show(io::IO, chain::Chain) = show(io, getfield(chain, :_row))

@inline Base.parent(chain::Chain) = chain._sys
@inline parent_system(chain::Chain) = parent(chain)
@inline parent_molecule(chain::Chain) = molecule_by_idx(parent(chain), chain.molecule_idx)

@doc raw"""
    parent_chain(::Atom)
    parent_chain(::Fragment)
    parent_chain(::Nucleotide)
    parent_chain(::Residue)

Returns the `Chain{T}` containing the given object. Returns `nothing` if no such chain exists.
""" parent_chain

"""
    $(TYPEDSIGNATURES)

Returns the `Chain{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
chain exists.
"""
@inline function chain_by_idx(sys::System{T}, idx::Int) where T
    Chain{T}(sys, _row_by_idx(sys._chains, idx))
end

"""
    chains(::Molecule)
    chains(::System; kwargs...)

Returns a `ChainTable{T}` containing all chains of the given atom container.

# Supported keyword arguments
 - `molecule_idx::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to chains belonging to the molecule with the given ID.
"""
@inline function chains(sys::System{T}; molecule_idx::MaybeInt = nothing) where T
    isnothing(molecule_idx) && return ChainTable{T}(sys, sys._chains.idx)
    _filter_chains(chain -> chain.molecule_idx == molecule_idx, sys)
end

"""
    nchains(::Molecule)
    nchains(::System; kwargs...)

Returns the number of chains in the given atom container.

# Supported keyword arguments
See [`chains`](@ref)
"""
@inline function nchains(sys::System; kwargs...)
    length(chains(sys; kwargs...))
end

#=
    Molecule chains
=#
@inline chains(mol::Molecule) = chains(parent(mol), molecule_idx = mol.idx)
@inline nchains(mol::Molecule) = nchains(parent(mol), molecule_idx = mol.idx)

"""
    push!(::Molecule{T}, ::Chain{T})

Creates a copy of the given chain in the given molecule. The new chain is automatically assigned a
new `idx`.
"""
@inline function Base.push!(mol::Molecule{T}, chain::Chain{T}) where T
    Chain(mol;
        name = chain.name,
        properties = chain.properties,
        flags = chain.flags
    )
    mol
end

#=
    Chain atoms
=#
@inline atoms(chain::Chain; kwargs...) = atoms(parent(chain); chain_idx = chain.idx, kwargs...)
@inline natoms(chain::Chain; kwargs...) = natoms(parent(chain); chain_idx = chain.idx, kwargs...)

@inline function Atom(chain::Chain, number::Int, element::ElementType; kwargs...)
    Atom(parent(chain), number, element;
        molecule_idx = chain.molecule_idx,
        chain_idx = chain.idx,
        kwargs...
    )
end

@inline function Base.push!(chain::Chain{T}, atom::Atom{T}; kwargs...) where T
    push!(parent(chain), atom;
        molecule_idx = chain.molecule_idx,
        chain_idx = chain.idx,
        kwargs...
    )
    chain
end

#=
    Chain bonds
=#
@inline bonds(chain::Chain; kwargs...) = bonds(parent(chain); chain_idx = chain.idx, kwargs...)
@inline nbonds(chain::Chain; kwargs...) = nbonds(parent(chain); chain_idx = chain.idx, kwargs...)
