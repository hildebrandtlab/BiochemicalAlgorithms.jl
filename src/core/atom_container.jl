export AbstractAtomContainer, frame_ids, nframes, 
    has_property, get_property, set_property

"""
    $(TYPEDEF)

Abstract base type for all atom containers.
"""
abstract type AbstractAtomContainer{T} <: AbstractSystemComponent{T} end

frame_ids(ac::AbstractAtomContainer{T}) where {T<:Real} = unique(_atoms(ac).frame_id)
nframes(ac::AbstractAtomContainer{T}) where {T<:Real} = length(frame_ids(ac))


function has_property(ac::AbstractAtomContainer{T}, key::Symbol) where {T<:Real}
    haskey(ac.properties, key)
end

function get_property(ac::AbstractAtomContainer{T}, key::Symbol) where {T<:Real}
    ac.properties[key]
end

function get_property(ac::AbstractAtomContainer{T}, key::Symbol, default) where {T<:Real}
    has_property(ac, key) ? ac.properties[key] : default
end

function set_property(ac::AbstractAtomContainer{T}, key::Symbol, value) where {T<:Real}
    ac.properties[key] = value
end
