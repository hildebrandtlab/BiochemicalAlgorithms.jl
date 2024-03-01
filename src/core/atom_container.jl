export
    AbstractAtomContainer,
    frame_ids,
    nframes

"""
    $(TYPEDEF)

Abstract base type for all atom containers.
"""
abstract type AbstractAtomContainer{T} <: AbstractSystemComponent{T} end

@inline frame_ids(ac::AbstractAtomContainer) = unique(atoms(ac).frame_id)
@inline nframes(ac::AbstractAtomContainer) = length(frame_ids(ac))
