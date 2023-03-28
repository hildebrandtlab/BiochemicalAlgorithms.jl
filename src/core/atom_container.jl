export AbstractAtomContainer, frame_ids, nframes, 
    has_property, get_property, set_property

"""
    $(TYPEDEF)

Abstract base type for all atom containers.
"""
abstract type AbstractAtomContainer{T <: Real} end

frame_ids(ac::AbstractAtomContainer{T}) where {T<:Real} = unique(_atoms(ac).frame_id)
nframes(ac::AbstractAtomContainer{T}) where {T<:Real} = length(frame_ids(ac))


function has_property(ac::AbstractAtomContainer{T}, key::String) where {T<:Real}
    haskey(ac.properties, key)
end

function get_property(ac::AbstractAtomContainer{T}, key::String) where {T<:Real}
    ac.properties[key]
end

function get_property(ac::AbstractAtomContainer{T}, key::String, default) where {T<:Real}
    has_property(ac, key) ? ac.properties[key] : default
end

function set_property(ac::AbstractAtomContainer{T}, key::String, value) where {T<:Real}
    ac.properties[key] = value
end
