export
    Chain,
    ChainTable,
    chain_by_idx,
    chains,
    chains_df,
    eachchain,
    nchains,
    parent_chain

@auto_hash_equals struct ChainTable{T} <: Tables.AbstractColumns
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _chains(ct::ChainTable) = getproperty(getfield(ct, :_sys), :_chains)

@inline Tables.istable(::Type{<: ChainTable}) = true
@inline Tables.columnaccess(::Type{<: ChainTable}) = true
@inline Tables.columns(ct::ChainTable) = ct

@inline function Tables.getcolumn(ct::ChainTable, nm::Symbol)
    col = Tables.getcolumn(_chains(ct), nm)
    RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _chains(ct)._idx_map[idx], getfield(ct, :_idx))
    )
end

@inline function Base.getproperty(ct::ChainTable, nm::Symbol)
    hasfield(typeof(ct), nm) && return getfield(ct, nm)
    Tables.getcolumn(ct, nm)
end

@inline Tables.getcolumn(ct::ChainTable, i::Int) = Tables.getcolumn(ct, Tables.columnnames(ct)[i])
@inline Tables.columnnames(ct::ChainTable) = Tables.columnnames(_chains(ct))
@inline Tables.schema(ct::ChainTable) = Tables.schema(_chains(ct))

@inline Base.size(ct::ChainTable) = (length(getfield(ct, :_idx)), length(_chain_table_cols))
@inline Base.size(ct::ChainTable, dim) = size(ct)[dim]
@inline Base.length(ct::ChainTable) = size(ct, 1)

function Base.push!(ct::ChainTable, t::ChainTuple, molecule_id::Int)
    sys = getfield(ct, :_sys)
    push!(sys._chains, t, molecule_id)
    push!(getfield(ct, :_idx), sys._curr_idx)
    ct
end

@inline function _filter_chains(f::Function, sys::System{T}) where T
    ChainTable{T}(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._chains),
        :idx
    )))
end

@inline function Base.filter(f::Function, ct::ChainTable)
    ChainTable(getfield(ct, :_sys), collect(Int, _filter_select(
        TableOperations.filter(f, ct),
        :idx
    )))
end

@inline function Base.iterate(ct::ChainTable, st = 1)
    st > length(ct) ?
        nothing :
        (chain_by_idx(getfield(ct, :_sys), getfield(ct, :_idx)[st]), st + 1)
end
@inline Base.eltype(::ChainTable{T}) where T = Chain{T}
@inline Base.getindex(ct::ChainTable{T}, i::Int) where T = chain_by_idx(getfield(ct, :_sys), getfield(ct, :_idx)[i])
@inline Base.keys(ct::ChainTable) = LinearIndices((length(ct),))

"""
    $(TYPEDEF)

Mutable representation of an individual chain in a system.

# Fields
 - `idx::Int`
 - `name::String`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Chain(
    mol::Molecule{T},
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

function Chain(
    mol::Molecule{T},
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
) where T
    sys = parent(mol)
    idx = _next_idx(sys)
    push!(sys._chains, ChainTuple(
            idx = idx,
            name = name,
            properties = properties,
            flags = flags
        ), mol.idx
    )
    chain_by_idx(sys, idx)
end

@inline Tables.rows(ct::ChainTable) = ct
@inline Tables.getcolumn(chain::Chain, name::Symbol) = Tables.getcolumn(getfield(chain, :_row), name)

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
@inline parent_molecule(chain::Chain) = molecule_by_idx(parent(chain), chain._row.molecule_id)

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
    $(TYPEDSIGNATURES)

Returns a raw `ChainTable` for all of the given system's chains matching the given criteria. Fields
given as `nothing` are ignored. The returned table contains all public and private chain fields.
"""
@inline function _chains(sys::System{T}; molecule_id::MaybeInt = nothing) where T
    isnothing(molecule_id) && return ChainTable{T}(sys, sys._chains.idx)
    _filter_chains(chain -> chain.molecule_id == molecule_id, sys)
end

"""
    chains(::Molecule)
    chains(::Protein)
    chains(::System)

Returns a `ChainTable{T}` containing all chains of the given atom container.

# Supported keyword arguments
 - `molecule_id::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to chains belonging to the molecule with the given ID.
"""
@inline function chains(sys::System; kwargs...)
    _chains(sys; kwargs...)
end

"""
    chains_df(::Molecule)
    chains_df(::Protein)
    chains_df(::System)

Returns a `DataFrame` containing all chains of the given atom container.
"""
@inline function chains_df(sys::System; kwargs...)
    DataFrame(_chains(sys; kwargs...))
end

"""
    eachchain(::Molecule)
    eachchain(::Protein)
    eachchain(::System)

Returns a `Chain{T}` generator for all chains of the given atom container.
"""
@inline function eachchain(sys::System{T}; kwargs...) where T
    (chain for chain in _chains(sys; kwargs...))
end

"""
    nchains(::Molecule)
    nchains(::Protein)
    nchains(::System)

Returns the number of chains in the given atom container.
"""
@inline function nchains(sys::System; kwargs...)
    length(_chains(sys; kwargs...))
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
    sys = parent(mol)
    push!(sys._chains, (; chain..., idx = _next_idx(sys)), mol.idx)
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
