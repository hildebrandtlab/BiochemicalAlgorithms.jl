export AbstractAtomContainer, frame_ids, nframes

"""
    $(TYPEDEF)

Abstract base type for all atom containers.
"""
abstract type AbstractAtomContainer{T <: Real} end

frame_ids(ac::AbstractAtomContainer{T}) where {T<:Real} = unique(_atoms(ac).frame_id)
nframes(ac::AbstractAtomContainer{T}) where {T<:Real} = length(frame_ids(ac))

