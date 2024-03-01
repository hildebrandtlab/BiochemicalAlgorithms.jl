export
    Nucleotide,
    NucleotideTable,
    nucleotide_by_idx,
    nucleotides,
    nnucleotides,
    parent_nucleotide

@auto_hash_equals struct NucleotideTable{T} <: Tables.AbstractColumns
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _nucleotides(nt::NucleotideTable) = getproperty(getfield(nt, :_sys), :_nucleotides)

@inline Tables.istable(::Type{<: NucleotideTable}) = true
@inline Tables.columnaccess(::Type{<: NucleotideTable}) = true
@inline Tables.columns(nt::NucleotideTable) = nt

@inline function Tables.getcolumn(nt::NucleotideTable, nm::Symbol)
    col = Tables.getcolumn(_nucleotides(nt), nm)
    RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _nucleotides(nt)._idx_map[idx], getfield(nt, :_idx))
    )
end

@inline function Base.getproperty(nt::NucleotideTable, nm::Symbol)
    hasfield(typeof(nt), nm) && return getfield(nt, nm)
    Tables.getcolumn(nt, nm)
end

@inline Tables.getcolumn(nt::NucleotideTable, i::Int) = Tables.getcolumn(nt, Tables.columnnames(nt)[i])
@inline Tables.columnnames(nt::NucleotideTable) = Tables.columnnames(_nucleotides(nt))
@inline Tables.schema(nt::NucleotideTable) = Tables.schema(_nucleotides(nt))

@inline Base.size(nt::NucleotideTable) = (length(getfield(nt, :_idx)), length(_nucleotide_table_cols))
@inline Base.size(nt::NucleotideTable, dim) = size(nt)[dim]
@inline Base.length(nt::NucleotideTable) = size(nt, 1)

@inline function _filter_nucleotides(f::Function, sys::System{T}) where T
    NucleotideTable{T}(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._nucleotides),
        :idx
    )))
end

@inline function Base.filter(f::Function, nt::NucleotideTable)
    NucleotideTable(getfield(nt, :_sys), collect(Int, _filter_select(
        TableOperations.filter(f, nt),
        :idx
    )))
end

@inline function Base.iterate(nt::NucleotideTable, st = 1)
    st > length(nt) ?
        nothing :
        (nucleotide_by_idx(getfield(nt, :_sys), getfield(nt, :_idx)[st]), st + 1)
end
@inline Base.eltype(::NucleotideTable{T}) where T = Nucleotide{T}
@inline Base.getindex(nt::NucleotideTable{T}, i::Int) where T = nucleotide_by_idx(getfield(nt, :_sys), getfield(nt, :_idx)[i])
@inline Base.keys(nt::NucleotideTable) = LinearIndices((length(nt),))

"""
    $(TYPEDEF)

Mutable representation of an individual nucleotide in a system.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Nucleotide(
    chain::Chain{T},
    number::Int,
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Nucleotide{T}` in the given chain.
"""
@auto_hash_equals struct Nucleotide{T} <: AbstractAtomContainer{T}
    _sys::System{T}
    _row::_NucleotideTableRow
end

@inline function Nucleotide(
    chain::Chain{T},
    number::Int;
    kwargs...
) where T
    sys = parent(chain)
    idx = _next_idx(sys)
    push!(
        sys._nucleotides,
        _Nucleotide(number; idx = idx, kwargs...),
        chain.molecule_id,
        chain.idx
    )
    nucleotide_by_idx(sys, idx)
end

@inline Tables.rows(nt::NucleotideTable) = nt
@inline Tables.getcolumn(nuc::Nucleotide, name::Symbol) = Tables.getcolumn(getfield(nuc, :_row), name)

@inline function Base.getproperty(nuc::Nucleotide, name::Symbol)
    hasfield(typeof(nuc), name) && return getfield(nuc, name)
    getproperty(getfield(nuc, :_row), name)
end

@inline function Base.setproperty!(nuc::Nucleotide, name::Symbol, val)
    hasfield(typeof(nuc), name) && return setfield!(nuc, name, val)
    setproperty!(getfield(nuc, :_row), name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", nuc::Nucleotide) = show(io, getfield(nuc, :_row))
@inline Base.show(io::IO, nuc::Nucleotide) = show(io, getfield(nuc, :_row))

@inline Base.parent(nuc::Nucleotide) = nuc._sys
@inline parent_system(nuc::Nucleotide) = parent(nuc)
@inline parent_molecule(nuc::Nucleotide) = molecule_by_idx(parent(nuc), nuc.molecule_id)
@inline parent_chain(nuc::Nucleotide) = chain_by_idx(parent(nuc), nuc.chain_id)

@doc raw"""
    parent_nucleotide(::Atom)

Returns the `Nucleotide{T}` containing the given atom. Returns `nothing` if no such nucleotide exists.
""" parent_nucleotide

"""
    $(TYPEDSIGNATURES)

Returns the `Nucleotide{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
nucleotide exists.
"""
@inline function nucleotide_by_idx(sys::System{T}, idx::Int) where T
    Nucleotide{T}(sys, _row_by_idx(sys._nucleotides, idx))
end

"""
    nucleotides(::Chain)
    nucleotides(::Molecule)
    nucleotides(::Protein)
    nucleotides(::System)

Returns a `NucleotideTable{T}` containing all nucleotides of the given atom container.

# Supported keyword arguments
 - `molecule_id::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to nucleotides belonging to the molecule with the given ID.
- `chain_id::MaybeInt = nothing`: \
Any value other than `nothing` limits the result to nucleotides belonging to the chain with the given ID.
"""
function nucleotides(sys::System{T};
    molecule_id::MaybeInt = nothing,
    chain_id::MaybeInt = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return NucleotideTable{T}(sys, sys._nucleotides.idx)
    _filter_nucleotides(nuc ->
        (isnothing(molecule_id) || nuc.molecule_id == something(molecule_id)) &&
        (isnothing(chain_id)    || nuc.chain_id    == something(chain_id)),
        sys
    )
end

"""
    nnucleotides(::Chain)
    nnucleotides(::Molecule)
    nnucleotides(::Protein)
    nnucleotides(::System)

Returns the number of nucleotides in the given atom container.
"""
@inline function nnucleotides(sys::System; kwargs...)
    length(nucleotides(sys; kwargs...))
end

#=
    Nucleotides
=#
@inline nucleotides(mol::Molecule; kwargs...) = nucleotides(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nnucleotides(mol::Molecule; kwargs...) = nnucleotides(parent(mol); molecule_id = mol.idx, kwargs...)

#=
    Chain nucleotides
=#
@inline nucleotides(chain::Chain; kwargs...) = nucleotides(parent(chain); chain_id = chain.idx, kwargs...)
@inline nnucleotides(chain::Chain; kwargs...) = nnucleotides(parent(chain); chain_id = chain.idx, kwargs...)

"""
    push!(::Chain{T}, nuc::Nucleotide{T})

Creates a copy of the given nucleotide in the given chain. The new nucleotide is automatically assigned a
new `idx`.
"""
@inline function Base.push!(chain::Chain{T}, nuc::Nucleotide{T}) where T
    Nucleotide(parent(chain), nuc.number;
        name = nuc.name,
        properties = nuc.properties,
        flags = nuc.flags
    )
    chain
end

#=
    Nucleotide atoms
=#
@inline atoms(nuc::Nucleotide; kwargs...) = atoms(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline natoms(nuc::Nucleotide; kwargs...) = natoms(parent(nuc); nucleotide_id = nuc.idx, kwargs...)

@inline function Atom(nuc::Nucleotide, number::Int, element::ElementType; kwargs...)
    Atom(parent(nuc), number, element;
        molecule_id = nuc.molecule_id,
        chain_id = nuc.chain_id,
        nucleotide_id = nuc.idx,
        kwargs...
    )
end

@inline function Base.push!(nuc::Nucleotide{T}, atom::Atom{T}; kwargs...) where T
    push!(parent(nuc), atom;
        molecule_id = nuc.molecule_id,
        chain_id = nuc.chain_id,
        nucleotide_id = nuc.idx,
        kwargs...
    )
    nuc
end

#=
    Nucleotide bonds
=#
@inline bonds(nuc::Nucleotide; kwargs...) = bonds(parent(nuc); nucleotide_id = nuc.idx, kwargs...)
@inline nbonds(nuc::Nucleotide; kwargs...) = nbonds(parent(nuc); nucleotide_id = nuc.idx, kwargs...)

# TODO: we should come up with a better test than just checking the name
is_nucleotide(name::String) = name ∈ ["A", "C", "G", "T", "U", "I", "DA", "DC", "DG", "DT", "DU", "DI"]
