export
    SecondaryStructure,
    SecondaryStructureTable,
    secondary_structure_by_idx,
    secondary_structures,
    nsecondary_structures

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
 - `first_fragment_idx::Int`
 - `last_fragment_idx::Int`

# Constructors
```julia
SecondaryStructure(
    first_fragment::Fragment{T},
    last_fragment::Fragment{T};
    number::Int,
    type::SecondaryStructureType;
    # keyword arguments
    name::AbstractString = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `SecondaryStructure{T}` in the given chain.
"""
const SecondaryStructure{T} = AtomContainer{T, :SecondaryStructure}

@inline function SecondaryStructure(
    first_fragment::Fragment{T},
    last_fragment::Fragment{T},
    number::Int,
    type::SecondaryStructureType;
    name::AbstractString="",
    kwargs...
) where T
    @assert first_fragment.chain_idx == last_fragment.chain_idx "Invalid secondary structure: fragments on different chains!"
    sys = parent(first_fragment)
    chain = parent_chain(first_fragment)
    idx = _next_idx!(sys)
    push!(sys._secondary_structures, idx, number, type, chain.molecule_idx, chain.idx, first_fragment.idx, last_fragment.idx; name=name, kwargs...)
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
 - `first_fragment_idx::AbstractVector{Int}`
 - `last_fragment_idx::AbstractVector{Int}`
"""
const SecondaryStructureTable{T} = SystemComponentTable{T, SecondaryStructure{T}}

@inline function _wrap_secondary_structures(sys::System{T}) where T
    SecondaryStructureTable{T}(sys, getfield(getfield(sys, :_secondary_structures), :idx))
end

@inline function _filter_secondary_structures(f::Function, sys::System{T}) where T
    SecondaryStructureTable{T}(sys, _filter_idx(f, sys._secondary_structures))
end

@inline _table(sys::System{T}, ::Type{SecondaryStructure{T}}) where T = sys._secondary_structures

@inline function _hascolumn(::Type{<: SecondaryStructure}, nm::Symbol)
    _hascolumn(_SecondaryStructureTable, nm)
end

@inline parent_molecule(secondary_structure::SecondaryStructure) = molecule_by_idx(parent(secondary_structure), secondary_structure.molecule_idx)
@inline parent_chain(secondary_structure::SecondaryStructure) = chain_by_idx(parent(secondary_structure), secondary_structure.chain_idx)

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
    secondary_structures(::Fragment)
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
    nsecondary_structures(::Fragment)
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
    secondary_structures(::FragmentTable)
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
    nsecondary_structures(::FragmentTable)
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

#=
    Fragment secondary structures
=#
@inline function secondary_structures(frag::Fragment; kwargs...)
    filter(ss -> frag.idx in fragments(ss).idx, secondary_structures(parent(frag); kwargs...))
end

@inline function secondary_structures(ft::FragmentTable; kwargs...)
    idx = Set(ft.idx)
    filter(ss -> !isempty(idx ∩ fragments(ss).idx), secondary_structures(ft._sys; kwargs...))
end

@inline function nsecondary_structures(frag::Union{Fragment, FragmentTable}; kwargs...)
    length(secondary_structures(frag; kwargs...))
end

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
    if !keep_fragments
        delete!(fragments(ss); keep_atoms = false)
    end

    delete!(parent(ss)._secondary_structures, ss.idx)
    nothing
end

function Base.delete!(st::SecondaryStructureTable; keep_fragments::Bool = false)
    if !keep_fragments
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

#=
    SecondaryStructure fragments
=#

@inline function fragments(ss::SecondaryStructure; kwargs...)
    sys = parent(ss)
    first_frag = fragment_by_idx(sys, ss.first_fragment_idx)
    last_frag = fragment_by_idx(sys, ss.last_fragment_idx)
    number_range = first_frag.number:last_frag.number
    filter(frag -> frag.number in number_range, fragments(parent_chain(ss); kwargs...))
end

@inline function nfragments(ss::SecondaryStructure; kwargs...)
    length(fragments(ss; kwargs...))
end

@inline function fragments(st::SecondaryStructureTable; kwargs...)
    idx = Set(Iterators.flatten(getproperty.(fragments.(secondary_structures(st._sys); kwargs...), :idx)))
    _filter_fragments(frag -> frag.idx in idx, st._sys)
end

@inline function nfragments(st::SecondaryStructureTable; kwargs...)
    length(fragments(st; kwargs...))
end

@inline function nucleotides(ac::Union{SecondaryStructure, SecondaryStructureTable}; kwargs...)
    fragments(ac; variant = FragmentVariant.Nucleotide, kwargs...)
end

@inline function nnucleotides(ac::Union{SecondaryStructure, SecondaryStructureTable}; kwargs...)
    nfragments(ac; variant = FragmentVariant.Nucleotide, kwargs...)
end

@inline function residues(ac::Union{SecondaryStructure, SecondaryStructureTable}; kwargs...)
    fragments(ac; variant = FragmentVariant.Residue, kwargs...)
end

@inline function nresidues(ac::Union{SecondaryStructure, SecondaryStructureTable}; kwargs...)
    nfragments(ac; variant = FragmentVariant.Residue, kwargs...)
end
