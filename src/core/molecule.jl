export
    AbstractMolecule,
    Molecule,
    MoleculeTable,
    molecule_by_idx,
    molecules,
    nmolecules,
    parent_molecule

"""
    $(TYPEDEF)

Tables.jl-compatible representation of system molecules (or a subset thereof). Molecule tables can be
generated using [`molecules`](@ref) or filtered from other molecule tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `name::AbstractVector{String}`

# Private columns
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
"""
@auto_hash_equals struct MoleculeTable{T} <: AbstractSystemComponentTable{T}
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _molecules(mt::MoleculeTable) = getfield(getfield(mt, :_sys), :_molecules)

@inline function Tables.getcolumn(mt::MoleculeTable, nm::Symbol)
    col = Tables.getcolumn(_molecules(mt), nm)
    _RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _molecules(mt)._idx_map[idx], mt._idx)
    )
end

@inline Tables.columnnames(mt::MoleculeTable) = Tables.columnnames(_molecules(mt))
@inline Tables.schema(mt::MoleculeTable) = Tables.schema(_molecules(mt))

@inline function Base.getproperty(mt::MoleculeTable, nm::Symbol)
    hasfield(typeof(mt), nm) && return getfield(mt, nm)
    Tables.getcolumn(mt, nm)
end

@inline function Base.setproperty!(mt::MoleculeTable, nm::Symbol, val)
    if nm in _molecule_table_cols_priv || nm in _molecule_table_cols_set
        error("MoleculeTable columns cannot be set directly! Did you mean to use broadcast assignment (.=)?")
    end
    if !hasfield(typeof(mt), nm)
        error("type MoleculeTable has no field $nm")
    end
    setfield!(mt, nm, val)
end

@inline function _filter_molecules(f::Function, sys::System)
    MoleculeTable(sys, _filter_idx(f, sys._molecules))
end

@inline function Base.filter(f::Function, mt::MoleculeTable)
    MoleculeTable(mt._sys, _filter_idx(f, mt))
end

@inline function Base.iterate(mt::MoleculeTable, st = 1)
    st > length(mt) ?
        nothing :
        (molecule_by_idx(mt._sys, mt._idx[st]), st + 1)
end

@inline Base.eltype(::MoleculeTable{T}) where T = Molecule{T}
@inline Base.size(mt::MoleculeTable) = (length(mt._idx), length(Tables.columnnames(mt)))
@inline Base.getindex(mt::MoleculeTable, i::Int) = molecule_by_idx(mt._sys, mt._idx[i])
@inline Base.getindex(mt::MoleculeTable, ::Colon) = mt

@inline function Base.getindex(mt::MoleculeTable, I)
    MoleculeTable(mt._sys, collect(Int, map(i -> mt._idx[i], I)))
end

    """
    $(TYPEDEF)

Abstract base type for all molecules.
"""
abstract type AbstractMolecule{T} <: AbstractAtomContainer{T} end

"""
    $(TYPEDEF)

Mutable representation of an individual molecule in a system.

# Public fields
 - `idx::Int`
 - `name::String`

# Private fields
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Molecule(
    sys::System{T};
    # keyword arguments
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Molecule{T}` in the given system.

```julia
Molecule(;
    #keyword arguments
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Molecule{Float32}` in the default system.
"""
@auto_hash_equals struct Molecule{T} <: AbstractMolecule{T}
    _sys::System{T}
    _row::_MoleculeTableRow
end

@inline function Molecule(
    sys::System{T};
    kwargs...
) where T
    idx = _next_idx(sys)
    push!(sys._molecules, idx; kwargs...)
    molecule_by_idx(sys, idx)
end

@inline function Molecule(; kwargs...)
    Molecule(default_system(); kwargs...)
end

@inline function Base.getproperty(mol::Molecule, name::Symbol)
    hasfield(typeof(mol), name) && return getfield(mol, name)
    getproperty(getfield(mol, :_row), name)
end

@inline function Base.setproperty!(mol::Molecule, name::Symbol, val)
    hasfield(typeof(mol), name) && return setfield!(mol, name, val)
    setproperty!(getfield(mol, :_row), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", mol::Molecule) = show(io, mol)
@inline function Base.show(io::IO, mol::Molecule)
    print(io, "$(typeof(mol)): ")
    show(io, NamedTuple(mol._row))
end

@inline Base.parent(mol::Molecule) = mol._sys
@inline parent_system(mol::Molecule) = parent(mol)

@doc raw"""
    parent_molecule(::Atom)
    parent_molecule(::Chain)
    parent_molecule(::Fragment)
    parent_molecule(::Nucleotide)
    parent_molecule(::Residue)

Returns the `Molecule{T}` containing the given object. Returns `nothing` if no such molecule exists.
""" parent_molecule

"""
    push!(::System{T}, ::Molecule{T})

Creates a copy of the given molecule in the given system. The new molecule is automatically assigned a
new `idx`.
"""
@inline function Base.push!(sys::System{T}, mol::Molecule{T}) where T
    Molecule(sys;
        name = mol.name,
        properties = mol.properties,
        flags = mol.flags
    )
    sys
end

"""
    $(TYPEDSIGNATURES)

Returns the `Molecule{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
molecule exists.
"""
@inline function molecule_by_idx(sys::System{T}, idx::Int) where T
    Molecule{T}(sys, _row_by_idx(sys._molecules, idx))
end

"""
    $(TYPEDSIGNATURES)

Returns a `MoleculeTable{T}` containing all molecules of the given system.
"""
@inline function molecules(sys::System{T}) where T
    MoleculeTable{T}(sys, sys._molecules.idx)
end

"""
    $(TYPEDSIGNATURES)

Returns the number of molecules in the given system.
"""
function nmolecules(sys::System)
    length(sys._molecules)
end

#=
    Molecule atoms
=#
@inline atoms(mol::Molecule; kwargs...) = atoms(parent(mol); molecule_idx = mol.idx, kwargs...)
@inline natoms(mol::Molecule; kwargs...) = natoms(parent(mol); molecule_idx = mol.idx, kwargs...)

@inline function Atom(mol::Molecule, number::Int, element::ElementType; kwargs...)
    Atom(parent(mol), number, element; molecule_idx = mol.idx, kwargs...)
end

@inline function Base.push!(mol::Molecule{T}, atom::Atom{T}; kwargs...) where T
    push!(parent(mol), atom; molecule_idx = mol.idx, kwargs...)
    mol
end

#=
    Molecule bonds
=#
@inline bonds(mol::Molecule; kwargs...) = bonds(parent(mol); molecule_idx = mol.idx, kwargs...)
@inline nbonds(mol::Molecule; kwargs...) = nbonds(parent(mol); molecule_idx = mol.idx, kwargs...)
