export
    Chain,
    ChainTable,
    chain_by_idx,
    chains,
    nchains,
    parent_chain

"""
    Chain{T} <: AbstractAtomContainer{T}

Mutable representation of an individual chain in a system.

# Public fields
 - `idx::Int`
 - `name::String`

# Private fields
 - `properties::Properties`
 - `flags::Flags`
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
const Chain{T} = AtomContainer{T, :Chain}

@inline function Chain(
    mol::Molecule;
    kwargs...
)
    sys = parent(mol)
    idx = _next_idx(sys)
    push!(sys._chains, idx, mol.idx; kwargs...)
    chain_by_idx(sys, idx)
end

"""
    ChainTable{T} <: AbstractSystemComponentTable{T}

Tables.jl-compatible representation of system chains (or a subset thereof). Chain tables can be
generated using [`chains`](@ref) or filtered from other chain tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `name::AbstractVector{String}`

# Private columns
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
 - `molecule_idx::AbstractVector{Int}`
"""
const ChainTable{T} = SystemComponentTable{T, Chain{T}}

@inline function _filter_chains(f::Function, sys::System{T}) where T
    ChainTable{T}(sys, _filter_idx(f, sys._chains))
end

@inline _table(sys::System{T}, ::Type{Chain{T}}) where T = sys._chains

@inline function _hascolumn(::Type{<: Chain}, nm::Symbol)
    nm in _chain_table_cols_set || nm in _chain_table_cols_priv
end

@inline parent_molecule(chain::Chain) = molecule_by_idx(parent(chain), chain.molecule_idx)

@doc raw"""
    parent_chain(::Atom)
    parent_chain(::Fragment)

Returns the `Chain{T}` containing the given object. Returns `nothing` if no such chain exists.
""" parent_chain

"""
    chain_by_idx(
        sys::System{T} = default_system(),
        idx::Int
    ) -> Chain{T}

Returns the `Chain{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
chain exists.
"""
@inline function chain_by_idx(sys::System{T}, idx::Int) where T
    _rowno_by_idx(_table(sys, Chain{T}), idx) # check idx
    Chain{T}(sys, idx)
end

@inline function chain_by_idx(idx::Int)
    chain_by_idx(default_system(), idx)
end

"""
    chains(::Molecule)
    chains(::System = default_system(); kwargs...)

Returns a `ChainTable{T}` containing all chains of the given atom container.

# Supported keyword arguments
 - `molecule_idx::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to chains belonging to the molecule with the given ID.
"""
@inline function chains(
    sys::System{T} = default_system();
    molecule_idx::MaybeInt = nothing
) where T
    isnothing(molecule_idx) && return ChainTable{T}(sys, sys._chains.idx)
    _filter_chains(chain -> chain.molecule_idx == molecule_idx, sys)
end

"""
    nchains(::Molecule)
    nchains(::System = default_system(); kwargs...)

Returns the number of chains in the given atom container.

# Supported keyword arguments
See [`chains`](@ref)
"""
@inline function nchains(sys::System = default_system(); kwargs...)
    length(chains(sys; kwargs...))
end

#=
    Molecule chains
=#
@inline chains(mol::Molecule) = chains(parent(mol), molecule_idx = mol.idx)
@inline nchains(mol::Molecule) = nchains(parent(mol), molecule_idx = mol.idx)

"""
    chains(::MoleculeTable)

Returns a `ChainTable{T}` containing all chains of the given molecule table.
"""
@inline function chains(mt::MoleculeTable)
    idx = Set(mt.idx)
    _filter_chains(chain -> chain.molecule_idx in idx, mt._sys)
end

"""
    nchains(::ChainTable)
    nchains(::MolculeTable)

Returns the number of chains belonging to the given molecule table.
"""
@inline function nchains(ct::ChainTable)
    length(ct)
end

@inline function nchains(mt::MoleculeTable)
    length(chains(mt))
end

"""
    delete!(::Chain)
    delete!(::ChainTable)

Removes the given chain(s) and all associated fragments from the associated system.

# Supported keyword arguments
 - `keep_atoms::Bool = false`
   Determines whether associated atoms (and their bonds) are removed as well
"""
function Base.delete!(chain::Chain; keep_atoms::Bool = false)
    keep_atoms ? atoms(chain).chain_idx .= Ref(nothing) : delete!(atoms(chain))
    delete!(fragments(chain); keep_atoms = keep_atoms)
    delete!(parent(chain)._chains, chain.idx)
    nothing
end

function Base.delete!(ct::ChainTable; keep_atoms::Bool = false)
    keep_atoms ? atoms(ct).chain_idx .= Ref(nothing) : delete!(atoms(ct))
    delete!(fragments(ct); keep_atoms = keep_atoms)
    delete!(_table(ct), ct._idx)
    empty!(ct._idx)
    ct
end

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

@inline function atoms(ct::ChainTable)
    idx = Set(ct._idx)
    _filter_atoms(atom -> atom.chain_idx in idx, ct._sys)
end

@inline natoms(ct::ChainTable) = length(atoms(ct))

#=
    Chain bonds
=#
@inline bonds(chain::Chain; kwargs...) = bonds(parent(chain); chain_idx = chain.idx, kwargs...)
@inline nbonds(chain::Chain; kwargs...) = nbonds(parent(chain); chain_idx = chain.idx, kwargs...)

@inline bonds(ct::ChainTable) = bonds(atoms(ct))
@inline nbonds(ct::ChainTable) = nbonds(atoms(ct))
