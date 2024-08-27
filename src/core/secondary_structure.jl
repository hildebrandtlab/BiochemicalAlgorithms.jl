export
    SecondaryStructure,
    SecondaryStructureTable,
    secondary_structure_by_idx,
    secondary_structures,
    nsecondary_structures,
    parent_secondary_structure

"""
    SecondaryStructure{T} <: AbstractAtomContainer{T}

Mutable representation of an individual secondary structure element in a Chain.

# Public fields
 - `idx::Int`
 - `number::Int`
 - `type::SecondaryStructureType`
 - `name::String`

# Private fields
 - `properties::Properties`
 - `flags::Flags`
 - `molecule_idx::Int`
 - `chain_idx::Int`

# Constructors
```julia
SecondaryStructure(
    chain::Chain{T};
    number::Int,
    type::SecondaryStructureType;
    # keyword arguments
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `SecondaryStructure{T}` in the given chain.
"""
const SecondaryStructure{T} = AtomContainer{T, :SecondaryStructure}

@inline function SecondaryStructure(
    chain::Chain,
    number::Int,
    type::SecondaryStructureType;
    name::String="",
    kwargs...
)
    sys = parent(chain)
    idx = _next_idx!(sys)
    push!(sys._secondary_structures, idx, number, type, chain.molecule_idx, chain.idx; name=name, kwargs...)
    secondary_structure_by_idx(sys, idx)
end

"""
    SecondaryStructureTable{T} <: AbstractSystemComponentTable{T}

Tables.jl-compatible representation of system secondary structures (or a subset thereof). Secondary structure tables can be
generated using [`secondary_structures`](@ref) or filtered from other secondary structure tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `number::Int`
 - `type::SecondaryStructureType`
 - `name::AbstractVector{String}`

# Private columns
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
 - `molecule_idx::AbstractVector{Int}`
 - `chain_idx::AbstractVector{Int}`
"""
const SecondaryStructureTable{T} = SystemComponentTable{T, SecondaryStructure{T}}

@inline function _filter_secondary_structures(f::Function, sys::System{T}) where T
    SecondaryStructureTable{T}(sys, _filter_idx(f, sys._secondary_structures))
end

@inline _table(sys::System{T}, ::Type{SecondaryStructure{T}}) where T = sys._secondary_structures

@inline function _hascolumn(::Type{<: SecondaryStructure}, nm::Symbol)
    _hascolumn(_SecondaryStructureTable, nm)
end

@inline parent_molecule(secondary_structure::SecondaryStructure) = molecule_by_idx(parent(secondary_structure), secondary_structure.molecule_idx)
@inline parent_chain(secondary_structure::SecondaryStructure) = chain_by_idx(parent(secondary_structure), secondary_structure.chain_idx)

@doc raw"""
    parent_secondary_structure(::Atom)
    parent_secondary_structure(::Fragment)
    parent_secondary_structure(::Nucleotide)
    parent_secondary_structure(::Residue)

Returns the `SecondaryStructure{T}` containing the given object. Returns `nothing` if no such secondary structure exists.
""" parent_secondary_structure

"""
    $(TYPEDSIGNATURES)

Returns the `SecondaryStructure{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
secondary structure exists.
"""
@inline function secondary_structure_by_idx(sys::System{T}, idx::Int) where T
    _rowno_by_idx(_table(sys, SecondaryStructure{T}), idx) # check idx
    SecondaryStructure{T}(sys, idx)
end

@inline function secondary_structure_by_idx(idx::Int)
    secondary_structure_by_idx(default_system(), idx)
end


"""
    secondary_structures(::Molecule)
    secondary_structures(::Chain)
    secondary_structures(::System; kwargs...)

Returns a `SecondaryStructureTable{T}` containing all secondary structures of the given atom container.

# Supported keyword arguments
 - `molecule_idx::MaybeInt = nothing`
 - `chain_idx::MaybeInt = nothing`
All keyword arguments limit the results to secondary structures matching the given IDs. Keyword arguments set to
`nothing` are ignored.
"""
@inline function secondary_structures(sys::System{T}; molecule_idx::MaybeInt = nothing, chain_idx::MaybeInt = nothing) where T
    isnothing(molecule_idx) &&
        isnothing(chain_idx) &&
        return SecondaryStructureTable{T}(sys, copy(sys._secondary_structures.idx))
    _filter_secondary_structures(ss ->
        (isnothing(molecule_idx) || ss.molecule_idx == something(molecule_idx)) &&
        (isnothing(chain_idx)    || ss.chain_idx    == something(chain_idx)),
        sys
    )
end

"""
    nsecondary_structures(::Chain)
    nsecondary_structures(::Molecule)
    nsecondary_structures(::System; kwargs...)

Returns the number of secondary structures in the given atom container.

# Supported keyword arguments
See [`secondary_structures`](@ref)
"""
@inline function nsecondary_structures(sys::System; kwargs...)
    length(secondary_structures(sys; kwargs...))
end

#=
    Molecule secondary structures
=#
@inline secondary_structures(mol::Molecule; kwargs...) = secondary_structures(parent(mol); molecule_idx = mol.idx, kwargs...)
@inline nsecondary_structures(mol::Molecule; kwargs...) = nsecondary_structures(parent(mol); molecule_idx = mol.idx, kwargs...)

"""
    secondary_structures(::ChainTable)
    secondary_structures(::MoleculeTable)

Returns a `SecondaryStructureTable{T}` containing all secondary structures of the given table.

# Supported keyword arguments
See [`secondary_structures`](@ref secondary_structures)
"""
@inline function secondary_structures(mt::MoleculeTable)
    idx = Set(mt.idx)

    _filter_secondary_structures(ss -> ss.molecule_idx in idx, mt._sys)
end

"""
    nsecondary_structures(::ChainTable)
    nsecondary_structures(::SecondaryStructureTable)
    nsecondary_structures(::MoleculeTable)

Returns the number of secondary_structures belonging to the given table.

# Supported keyword arguments
See [`secondary_structures`](@ref secondary_structures)
"""
@inline function nsecondary_structures(st::SecondaryStructureTable)
    length(st)
end

@inline function nsecondary_structures(mt::MoleculeTable; kwargs...)
    length(secondary_structures(mt; kwargs...))
end

#=
    Chain secondary structures
=#
@inline secondary_structures(chain::Chain; kwargs...) = secondary_structures(parent(chain); chain_idx = chain.idx, kwargs...)
@inline nsecondary_structures(chain::Chain; kwargs...) = nsecondary_structures(parent(chain); chain_idx = chain.idx, kwargs...)

@inline function secondary_structures(ct::ChainTable)
    idx = Set(ct.idx)
    _filter_secondary_structures(ss -> ss.chain_idx in idx, ct._sys)
end

@inline nsecondary_structures(ct::ChainTable; kwargs...) = length(secondary_structures(ct; kwargs...))

"""
    delete!(::SecondaryStructure; keep_fragments::Bool = false)
    delete!(::SecondaryStructureTable; keep_fragments::Bool = false)
    delete!(::SecondaryStructureTable, idx::Int; keep_fragments::Bool = false)

Removes the given secondary_structure(s) from the associated system.

# Supported keyword arguments
- `keep_fragments::Bool = false`
   Determines whether the fragments contained in this secondary structure are removed as well. All atoms contained in those
   fragments are deleted as well.
"""
function Base.delete!(ss::SecondaryStructure; keep_fragments::Bool = false)
    if keep_fragments
        fragments(ss).secondary_structure_idx .= Ref(nothing)
    else
        delete!(fragments(ss); keep_atoms = false)
    end

    delete!(parent(ss)._secondary_structures, ss.idx)
    nothing
end

function Base.delete!(st::SecondaryStructureTable; keep_fragments::Bool = false)
    if keep_fragments
        fragments(st).secondary_structure_idx .= Ref(nothing)
    else
        delete!(fragments(st); keep_atoms = false)
    end

    delete!(_table(st), st._idx)
    empty!(st._idx)
    st
end

function Base.delete!(st::SecondaryStructureTable, idx::Int; kwargs...)
    idx in st._idx || throw(KeyError(idx))
    delete!(secondary_structure_by_idx(st._sys, idx); kwargs...)
    deleteat!(st._idx, findall(i -> i == idx, st._idx))
    st
end

"""
    push!(::Chain{T}, ::SecondaryStructure{T})

Creates a copy of the given secondary structure in the given chain. The new secondary structure is automatically assigned a
new `idx`.
"""
@inline function Base.push!(chain::Chain{T}, ss::SecondaryStructure{T}) where T
    SecondaryStructure(chain, ss.number, ss.type;
        name = ss.name,
        properties = ss.properties,
        flags = ss.flags
    )
    chain
end

#=
    Secondary structure atoms
=#
@inline atoms(ss::SecondaryStructure; kwargs...) = atoms(fragments(ss); kwargs...)
@inline natoms(ss::SecondaryStructure; kwargs...) = natoms(fragments(ss); kwargs...)

@inline atoms(st::SecondaryStructureTable; kwargs...) = atoms(fragments(st); kwargs...)
@inline natoms(st::SecondaryStructureTable; kwargs...) = natoms(fragments(st); kwargs...)

#=
    Secondary structure bonds
=#
@inline bonds(ss::SecondaryStructure; kwargs...) = bonds(fragments(ss); kwargs...)
@inline nbonds(ss::SecondaryStructure; kwargs...) = nbonds(fragments(ss); kwargs...)

@inline bonds(st::SecondaryStructureTable; kwargs...) = bonds(fragments(st); kwargs...)
@inline nbonds(st::SecondaryStructureTable; kwargs...) = nbonds(fragments(st); kwargs...)
