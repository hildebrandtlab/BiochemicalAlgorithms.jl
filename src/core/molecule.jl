export
    AbstractMolecule,
    Molecule,
    MoleculeTable,
    molecule_by_idx,
    molecules,
    molecules_df,
    eachmolecule,
    nmolecules,
    parent_molecule

@auto_hash_equals struct MoleculeTable{T} <: Tables.AbstractColumns
    _sys::System{T}
    _idx::Vector{Int}
end

@inline _molecules(mt::MoleculeTable) = getproperty(getfield(mt, :_sys), :_molecules)

@inline Tables.istable(::Type{<: MoleculeTable}) = true
@inline Tables.columnaccess(::Type{<: MoleculeTable}) = true
@inline Tables.columns(mt::MoleculeTable) = mt

@inline function Tables.getcolumn(mt::MoleculeTable, nm::Symbol)
    col = Tables.getcolumn(_molecules(mt), nm)
    RowProjectionVector{eltype(col)}(
        col,
        map(idx -> _molecules(mt)._idx_map[idx], getfield(mt, :_idx))
    )
end

@inline function Base.getproperty(mt::MoleculeTable, nm::Symbol)
    hasfield(typeof(mt), nm) && return getfield(mt, nm)
    Tables.getcolumn(mt, nm)
end

@inline Tables.getcolumn(mt::MoleculeTable, i::Int) = Tables.getcolumn(mt, Tables.columnnames(mt)[i])
@inline Tables.columnnames(mt::MoleculeTable) = Tables.columnnames(_molecules(mt))
@inline Tables.schema(mt::MoleculeTable) = Tables.schema(_molecules(mt))

@inline Base.size(mt::MoleculeTable) = (length(getfield(mt, :_idx)), length(_molecule_table_cols))
@inline Base.size(mt::MoleculeTable, dim) = size(mt)[dim]
@inline Base.length(mt::MoleculeTable) = size(mt, 1)

function Base.push!(mt::MoleculeTable, t::MoleculeTuple)
    sys = getfield(mt, :_sys)
    push!(sys._molecules, t)
    push!(getfield(mt, :_idx), sys._curr_idx)
    mt
end

@inline function _filter_molecules(f::Function, sys::System{T}) where T
    MoleculeTable(sys, collect(Int, _filter_select(
        TableOperations.filter(f, sys._molecules),
        :idx
    )))
end

@inline function Base.filter(f::Function, mt::MoleculeTable)
    MoleculeTable(getfield(mt, :_sys), collect(Int, _filter_select(
        TableOperations.filter(f, mt),
        :idx
    )))
end

@inline function Base.iterate(mt::MoleculeTable, st = 1)
    st > length(mt) ?
        nothing :
        (molecule_by_idx(getfield(mt, :_sys), getfield(mt, :_idx)[st]), st + 1)
end
@inline Base.eltype(::MoleculeTable{T}) where T = Molecule{T}
@inline Base.getindex(mt::MoleculeTable{T}, i::Int) where T = molecule_by_idx(getfield(mt, :_sys), getfield(mt, :_idx)[i])
@inline Base.keys(mt::MoleculeTable) = LinearIndices((length(mt),))

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
    _row::_MoleculeTableRow
end

function Molecule(
    sys::System{T},
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
) where T
    idx = _next_idx(sys)
    push!(sys._molecules, MoleculeTuple(
        idx = idx,
        name = name,
        properties = properties,
        flags = flags
    ))
    molecule_by_idx(sys, idx)
end

@inline function Molecule(name::String = "", properties::Properties = Properties(), flags::Flags = Flags())
    Molecule(default_system(), name, properties, flags)
end

@inline Tables.rows(mt::MoleculeTable) = mt
@inline Tables.getcolumn(mol::Molecule, name::Symbol) = Tables.getcolumn(getfield(mol, :_row), name)

@inline function Base.getproperty(mol::Molecule, name::Symbol)
    hasfield(typeof(mol), name) && return getfield(mol, name)
    getproperty(getfield(mol, :_row), name)
end

@inline function Base.setproperty!(mol::Molecule, name::Symbol, val)
    hasfield(typeof(mol), name) && return setfield!(mol, name, val)
    setproperty!(getfield(mol, :_row), name, val)
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

Returns a `DataFrame` containing all molecules of the given system.
"""
@inline function molecules_df(sys::System{T}) where T
    DataFrame(molecules(sys))
end

"""
    $(TYPEDSIGNATURES)

Returns a `Molecule{T}` generator for all molecules of the given system.
"""
@inline function eachmolecule(sys::System{T}) where T
    (mol for mol in molecules(sys))
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
@inline bonds(mol::Molecule; kwargs...) = bonds(parent(mol); molecule_id = mol.idx, kwargs...)
@inline bonds_df(mol::Molecule; kwargs...) = bonds_df(parent(mol); molecule_id = mol.idx, kwargs...)
@inline eachbond(mol::Molecule; kwargs...) = eachbond(parent(mol); molecule_id = mol.idx, kwargs...)
@inline nbonds(mol::Molecule; kwargs...) = nbonds(parent(mol); molecule_id = mol.idx, kwargs...)
