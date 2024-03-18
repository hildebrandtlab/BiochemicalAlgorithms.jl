export
    AbstractSystemComponent,
    get_property,
    has_flag,
    has_property,
    set_flag!,
    set_property!,
    unset_flag!

"""
    $(TYPEDEF)

Abstract base type for all components of a system, including the system itself.
"""
abstract type AbstractSystemComponent{T <: Real} end

#=
    Properties
=#

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

#=
    Flags
=#

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
