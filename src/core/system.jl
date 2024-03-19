export
    AbstractAtomContainer,
    AbstractSystemComponent,
    System,
    SystemComponent,
    default_system,
    frame_ids,
    get_property,
    has_flag,
    has_property,
    nframes,
    parent_system,
    set_flag!,
    set_property!,
    unset_flag!

"""
    $(TYPEDEF)

Abstract base type for all components of a system, including the system itself.
"""
abstract type AbstractSystemComponent{T <: Real} end

"""
    $(TYPEDSIGNATURES)

Returns a `Bool` indicating whether the given system component has the given property.
"""
@inline has_property(ac::AbstractSystemComponent, key::Symbol) = haskey(ac.properties, key)

"""
    $(TYPEDSIGNATURES)

Returns the property associated with the given key in `ac`.
"""
@inline get_property(ac::AbstractSystemComponent, key::Symbol) = ac.properties[key]

"""
    $(TYPEDSIGNATURES)

Returns the property associated with the given key in `ac`. If no such property exists, returns `default`.
"""
@inline get_property(ac::AbstractSystemComponent, key::Symbol, default) = get(ac.properties, key, default)

"""
    $(TYPEDSIGNATURES)

Sets the property associated with the given key in `ac` to the given `value`.
"""
@inline set_property!(ac::AbstractSystemComponent, key::Symbol, value) = ac.properties[key] = value

"""
    $(TYPEDSIGNATURES)

Returns a `Bool` indicating whether the given system component has the given flag.
"""
@inline has_flag(ac::AbstractSystemComponent, flag::Symbol) = flag in ac.flags

"""
    $(TYPEDSIGNATURES)

Adds the given flag to `ac`.
"""
@inline set_flag!(ac::AbstractSystemComponent, flag::Symbol) = push!(ac.flags, flag)

"""
    $(TYPEDSIGNATURES)

Removes the given flag from `ac`.
"""
@inline unset_flag!(ac::AbstractSystemComponent, flag::Symbol) = delete!(ac.flags, flag)

"""
    $(TYPEDEF)

Abstract base type for all atom containers.
"""
abstract type AbstractAtomContainer{T} <: AbstractSystemComponent{T} end

@inline frame_ids(ac::AbstractAtomContainer) = unique(atoms(ac).frame_id)
@inline nframes(ac::AbstractAtomContainer) = length(frame_ids(ac))

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

@auto_hash_equals struct SystemComponent{T, R <: _AbstractColumnTableRow} <: AbstractSystemComponent{T}
    _sys::System{T}
    _row::R
end

@inline function Base.getproperty(sc::SystemComponent, name::Symbol)
    (name === :_sys || name === :_row) && return getfield(sc, name)
    getproperty(getfield(sc, :_row), name)
end

@inline function Base.setproperty!(sc::SystemComponent, name::Symbol, val)
    (name === :_sys || name === :_row) && return setfield!(sc, name, val)
    setproperty!(getfield(sc, :_row), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", sc::SystemComponent) = show(io, sc)
@inline function Base.show(io::IO, sc::SystemComponent)
    print(io, "$(typeof(sc)): ")
    show(io, NamedTuple(sc._row))
end

@inline Base.parent(sc::SystemComponent) = sc._sys
@inline parent_system(sc::SystemComponent) = parent(sc)
