using AutoHashEquals
export AbstractMolecule, Molecule, molecule_by_idx, molecules, molecules_df, eachmolecule, nmolecules,
    parent_molecule, has_property, get_property, set_property

"""
    $(TYPEDEF)

Abstract base type for all molecules.
"""
abstract type AbstractMolecule{T} <: AbstractAtomContainer{T} end

"""
    $(TYPEDEF)

Mutable representation of an individual molecule in a system.

# Fields
 - `idx::Int`
 - `name::String`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
```julia
Molecule(
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Molecule{Float32}` in the default system.

```julia
Molecule(
    sys::System{T},
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
)
```
Creates a new `Molecule{T}` in the given system.
"""
@auto_hash_equals struct Molecule{T} <: AbstractMolecule{T}
    _sys::System{T}
    _row::DataFrameRow
end

function Molecule(
    sys::System{T},
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
) where T
    idx = _next_idx(sys)
    push!(sys._molecules, (idx = idx, name = name, properties = properties, flags = flags))
    molecule_by_idx(sys, idx)
end

function Molecule(name::String = "", properties::Properties = Properties(), flags::Flags = Flags())
    Molecule(default_system(), name, properties, flags)
end

function Base.getproperty(mol::Molecule, name::Symbol)
    in(name, fieldnames(MoleculeTuple)) && return getproperty(getfield(mol, :_row), name)
    getfield(mol, name)
end

function Base.setproperty!(mol::Molecule, name::Symbol, val)
    in(name, fieldnames(MoleculeTuple)) && return setproperty!(getfield(mol, :_row), name, val)
    setfield!(mol, name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", mol::Molecule) = show(io, getfield(mol, :_row))
@inline Base.show(io::IO, mol::Molecule) = show(io, getfield(mol, :_row))

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

# TODO this should also add all related entities
#function Base.push!(sys::System{T}, mol::Molecule{T}) where T
#    ...
#    sys
#end

"""
    $(TYPEDSIGNATURES)

Returns the `Molecule{T}` associated with the given `idx` in `sys`. Returns `nothing` if no such
molecule exists.
"""
@inline function molecule_by_idx(sys::System{T}, idx::Int) where T
    rn = _row_by_idx(sys._molecules, idx)
    isnothing(rn) ? nothing : Molecule{T}(sys, DataFrameRow(sys._molecules, rn, :))
end

"""
    $(TYPEDSIGNATURES)

Returns a raw `DataFrame` for all of the given system's molecules. The returned `DataFrame`
contains all public and private molecule fields.
"""
@inline function _molecules(sys::System)
    sys._molecules
end

"""
    $(TYPEDSIGNATURES)

Returns a `Vector{Molecule{T}}` containing all molecules of the given system.
"""
@inline function molecules(sys::System)
    collect(eachmolecule(sys))
end

"""
    $(TYPEDSIGNATURES)

Returns a `SubDataFrame` containing all molecules of the given system.
"""
@inline function molecules_df(sys::System{T}) where T
    view(_molecules(sys), :, :)
end

"""
    $(TYPEDSIGNATURES)

Returns a `Molecule{T}` generator for all molecules of the given system.
"""
@inline function eachmolecule(sys::System{T}) where T
    (Molecule{T}(sys, row) for row in eachrow(_molecules(sys)))
end

"""
    $(TYPEDSIGNATURES)

Returns the number of molecules in the given system.
"""
function nmolecules(sys::System)
    nrow(_molecules(sys))
end

#=
    Molecule atoms
=#
@inline _atoms(mol::Molecule; kwargs...) = _atoms(parent(mol), molecule_id = mol.idx, kwargs...)
@inline atoms(mol::Molecule; kwargs...) = atoms(parent(mol); molecule_id = mol.idx, kwargs...)
@inline atoms_df(mol::Molecule; kwargs...) = atoms_df(parent(mol); molecule_id = mol.idx, kwargs...)
@inline eachatom(mol::Molecule; kwargs...) = eachatom(parent(mol); molecule_id = mol.idx, kwargs...)
@inline natoms(mol::Molecule; kwargs...) = natoms(parent(mol); molecule_id = mol.idx, kwargs...)

@inline function Base.push!(mol::Molecule{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(parent(mol), atom; molecule_id = mol.idx, kwargs...)
    mol
end

#=
    Molecule bonds
=#
@inline _bonds(mol::Molecule; kwargs...) = _bonds(parent(mol); molecule_id = mol.idx, kwargs...)
@inline bonds(mol::Molecule; kwargs...) = bonds(parent(mol); molecule_id = mol.idx, kwargs...)
@inline bonds_df(mol::Molecule; kwargs...) = bonds_df(parent(mol); molecule_id = mol.idx, kwargs...)
@inline eachbond(mol::Molecule; kwargs...) = eachbond(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nbonds(mol::Molecule; kwargs...) = nbonds(parent(mol); molecule_id = mol.idx, kwargs...)

@inline function Base.push!(mol::Molecule, bond::BondTuple)
    push!(parent(mol), bond)
    mol
end

function has_property(m::AbstractMolecule, key::String)
    haskey(m.properties, key)
end

function get_property(m::AbstractMolecule, key::String)
    m.properties[key]
end

function set_property(m::AbstractMolecule, key::String, value)
    m.properties[key] = value
end
