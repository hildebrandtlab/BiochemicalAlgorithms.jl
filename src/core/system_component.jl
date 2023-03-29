export AbstractSystemComponent, has_property, get_property, set_property

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
@inline set_property(ac::AbstractSystemComponent, key::Symbol, value) = ac.properties[key] = value
