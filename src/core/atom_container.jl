export AbstractAtomContainer, frame_ids, nframes, 
    has_property, get_property, set_property

"""
    $(TYPEDEF)

Abstract base type for all atom containers.
"""
abstract type AbstractAtomContainer{T} <: AbstractSystemComponent{T} end

function Base.push!(ac::AbstractAtomContainer, bond::BondTuple)
    push!(parent(ac), bond)
    ac
end

frame_ids(ac::AbstractAtomContainer{T}) where {T<:Real} = unique(atoms(ac).frame_id)
nframes(ac::AbstractAtomContainer{T}) where {T<:Real} = length(frame_ids(ac))
