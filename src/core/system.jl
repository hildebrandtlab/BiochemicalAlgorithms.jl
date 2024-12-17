export
    AbstractAtomContainer,
    AbstractSystemComponent,
    AtomContainer,
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
    sort_atoms!,
    sort_bonds!,
    sort_chains!,
    sort_fragments!,
    sort_molecules!,
    sort_secondary_structures!,
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
    System(name::AbstractString = "", properties::Properties = Properties(), flags::Flags = Flags())

Creates a new and empty `System{Float32}`.

    System{T}(name::AbstractString = "", properties::Properties = Properties(), flags::Flags = Flags())

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
    _secondary_structures::_SecondaryStructureTable
    _fragments::_FragmentTable
    _curr_idx::Int

    function System{T}(
        name::AbstractString = "",
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
            _SecondaryStructureTable(),
            _FragmentTable(),
            0
        )
    end
end

System(
    name::AbstractString = "",
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
@inline function _next_idx!(sys::System{T}) where T
    sys._curr_idx += 1
end

Base.show(io::IO, ::MIME"text/plain", sys::System) = show(io, sys)
Base.show(io::IO, sys::System) = print(io,
    "$(typeof(sys)) with ", natoms(sys), " atoms", isempty(sys.name) ? "" : " ($(sys.name))")

"""
    empty!(::System)

Removes all components from the system.
"""
function Base.empty!(sys::System)
    empty!(sys.properties)
    empty!(sys.flags)
    empty!(sys._atoms)
    empty!(sys._bonds)
    empty!(sys._molecules)
    empty!(sys._chains)
    empty!(sys._secondary_structures)
    empty!(sys._fragments)
    sys
end

@doc raw"""
    parent(::Atom)
    parent(::Bond)
    parent(::Chain)
    parent(::SecondaryStructure)
    parent(::Fragment)
    parent(::Molecule)
    parent(::System)

Returns the `System{T}` containing the given object.
""" Base.parent(::System)
Base.parent(s::System) = s

@doc raw"""
    parent_system(::Atom)
    parent_system(::Bond)
    parent_system(::Chain)
    parent_system(::SecondaryStructure)
    parent_system(::Fragment)
    parent_system(::Molecule)
    parent_system(::System)

Returns the `System{T}` containing the given object. Alias for
[`Base.parent`](@ref Base.parent(::System)).
""" parent_system
parent_system(s::System) = s

@inline function _sort_table!(at::_AbstractSystemComponentTable, view; kwargs...)
    permute!(at,
        map(i -> at._idx_map[i], getproperty.(sort(view; by=e -> e.idx, kwargs...), :idx))
    )
end

"""
    sort_atoms!(::System)

Sorts the atoms in the given system by `idx` (default) or according to the given
keyword arguments.

# Supported keyword arguments
Same as `Base.sort`
"""
@inline function sort_atoms!(sys::System; kwargs...)
    _sort_table!(sys._atoms, _wrap_atoms(sys); kwargs...)
    sys
end

"""
    sort_bonds!(::System)

Sorts the bonds in the given system by `idx` (default) or according to the given
keyword arguments.

# Supported keyword arguments
Same as `Base.sort`
"""
@inline function sort_bonds!(sys::System; kwargs...)
    _sort_table!(sys._bonds, _wrap_bonds(sys); kwargs...)
    sys
end

"""
    sort_molecules!(::System)

Sorts the molecules in the given system by `idx` (default) or according to the given
keyword arguments.

# Supported keyword arguments
Same as `Base.sort`
"""
@inline function sort_molecules!(sys::System; kwargs...)
    _sort_table!(sys._molecules, _wrap_molecules(sys); kwargs...)
    sys
end

"""
    sort_chains!(::System)

Sorts the chains in the given system by `idx` (default) or according to the given
keyword arguments.

# Supported keyword arguments
Same as `Base.sort`
"""
@inline function sort_chains!(sys::System; kwargs...)
    _sort_table!(sys._chains, _wrap_chains(sys); kwargs...)
    sys
end

"""
    sort_secondary_structures!(::System)

Sorts the secondary structures in the given system by `idx` (default) or according
to the given keyword arguments.

# Supported keyword arguments
Same as `Base.sort`
"""
@inline function sort_secondary_structures!(sys::System; kwargs...)
    _sort_table!(sys._secondary_structures, _wrap_secondary_structures(sys); kwargs...)
    sys
end

"""
    sort_fragments!(::System)

Sorts the fragments in the given system by `idx` (default) or according to the given
keyword arguments.

# Supported keyword arguments
Same as `Base.sort`
"""
@inline function sort_fragments!(sys::System; kwargs...)
    _sort_table!(sys._fragments, _wrap_fragments(sys); kwargs...)
    sys
end

@auto_hash_equals struct SystemComponent{T, C} <: AbstractSystemComponent{T}
    _sys::System{T}
    _idx::Int
end

@inline _table(sc::SystemComponent{T, :Atom}) where T = getfield(getfield(sc, :_sys), :_atoms)
@inline _table(sc::SystemComponent{T, :Bond}) where T = getfield(getfield(sc, :_sys), :_bonds)
@inline _table(sc::SystemComponent{T, :Chain}) where T = getfield(getfield(sc, :_sys), :_chains)
@inline _table(sc::SystemComponent{T, :SecondaryStructure}) where T = getfield(getfield(sc, :_sys), :_secondary_structures)
@inline _table(sc::SystemComponent{T, :Fragment}) where T = getfield(getfield(sc, :_sys), :_fragments)
@inline _table(sc::SystemComponent{T, :Molecule}) where T = getfield(getfield(sc, :_sys), :_molecules)

@inline function Base.propertynames(sc::SystemComponent)
    propertynames(_table(sc))
end

@inline function Base.getproperty(sc::SystemComponent, name::Symbol)
    (name === :_sys || name === :_idx) && return getfield(sc, name)
    tab = _table(sc)
    getindex(getfield(tab, name), _rowno_by_idx(tab, getfield(sc, :_idx)))
end

@inline function Base.setproperty!(sc::SystemComponent, name::Symbol, val)
    (name === :_sys || name === :_idx) && return setfield!(sc, name, val)
    setindex!(getfield(_table(sc), name), val, _rowno_by_idx(_table(sc), getfield(sc, :_idx)))
end

@inline Base.show(io::IO, ::MIME"text/plain", sc::SystemComponent) = show(io, sc)
@inline function Base.show(io::IO, sc::SystemComponent; display_name::AbstractString=repr(typeof(sc)))
    print(io, "$display_name: ")
    show(io, NamedTuple(_row_by_idx(_table(sc), sc._idx)))
end

@inline Base.parent(sc::SystemComponent) = sc._sys
@inline parent_system(sc::SystemComponent) = parent(sc)

@auto_hash_equals struct AtomContainer{T, C} <: AbstractAtomContainer{T}
    _comp::SystemComponent{T, C}

    function AtomContainer{T, C}(
        sys::System{T},
        idx::Int
    ) where {T, C}
        new(SystemComponent{T, C}(sys, idx))
    end
end

@inline function Base.propertynames(ac::AtomContainer)
    propertynames(ac._comp)
end

@inline function Base.getproperty(ac::AtomContainer, name::Symbol)
    name === :_comp && return getfield(ac, :_comp)
    getproperty(getfield(ac, :_comp), name)
end

@inline function Base.setproperty!(ac::AtomContainer, name::Symbol, val)
    name === :_comp && return setfield!(ac, :_comp, val)
    setproperty!(getfield(ac, :_comp), name, val)
end

@inline Base.show(io::IO, ::MIME"text/plain", ac::AtomContainer) = show(io, ac)
@inline Base.show(io::IO, ac::AtomContainer) = show(io, ac._comp; display_name = repr(typeof(ac)))

@inline Base.parent(at::AtomContainer) = parent(at._comp)
@inline parent_system(at::AtomContainer) = parent_system(at._comp)
