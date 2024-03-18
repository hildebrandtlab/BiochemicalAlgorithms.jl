export
    System,
    default_system,
    parent_system

"""
    $(TYPEDEF)

Mutable representation of a biomolecular system.

# Fields
 - `name::String`
 - `properties::Properties`
 - `flags::Flags`

# Constructors
    System(name::String = "", properties::Properties = Properties(), flags::Flags = Flags())

Creates a new and empty `System{Float32}`.

    System{T}(name::String = "", properties::Properties = Properties(), flags::Flags = Flags())

Creates a new and empty `System{T}`.
"""
@auto_hash_equals mutable struct System{T} <: AbstractAtomContainer{T}
    name::String
    properties::Properties
    flags::Flags

    _atoms::_AtomTable{T}
    _bonds::_BondTable
    _molecules::_MoleculeTable
    _chains::_ChainTable
    _fragments::_FragmentTable
    _nucleotides::_NucleotideTable
    _residues::_ResidueTable
    _curr_idx::Int

    function System{T}(
        name::String = "",
        properties::Properties = Properties(),
        flags::Flags = Flags()
    ) where T
        new(
            name,
            properties,
            flags,
            _AtomTable{T}(),
            _BondTable(),
            _MoleculeTable(),
            _ChainTable(),
            _FragmentTable(),
            _NucleotideTable(),
            _ResidueTable(),
            0
        )
    end
end

System(
    name::String = "",
    properties::Properties = Properties(),
    flags::Flags = Flags()
) = System{Float32}(name, properties, flags)

"""
    const _default_system

Global default system.
"""
const _default_system = System("default")

"""
    $(TYPEDSIGNATURES)

Returns the global default system.
"""
@inline function default_system()
    _default_system
end

"""
    $(TYPEDSIGNATURES)

Returns the next available `idx` for the given system.
"""
@inline function _next_idx(sys::System{T}) where T
    sys._curr_idx += 1
end

Base.show(io::IO, ::MIME"text/plain", sys::System) = show(io, sys)
Base.show(io::IO, sys::System) = print(io, 
    "System with ", natoms(sys), " atoms", isempty(sys.name) ? "" : " ($(sys.name))")

@doc raw"""
    parent(::Atom)
    parent(::Bond)
    parent(::Chain)
    parent(::Fragment)
    parent(::Molecule)
    parent(::Nucleotide)
    parent(::Residue)
    parent(::System)

Returns the `System{T}` containing the given object.
""" Base.parent(::System)
Base.parent(s::System) = s

@doc raw"""
    parent_system(::Atom)
    parent_system(::Bond)
    parent_system(::Chain)
    parent_system(::Fragment)
    parent_system(::Molecule)
    parent_system(::Nucleotide)
    parent_system(::Residue)
    parent_system(::System)

Returns the `System{T}` containing the given object. Alias for 
[`Base.parent`](@ref Base.parent(::System)).
""" parent_system
parent_system(s::System) = s