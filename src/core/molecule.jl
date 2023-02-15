using DataFrames

export AbstractMolecule, Molecule

abstract type AbstractMolecule{T} <: AbstractAtomContainer{T} end

const _MoleculeAtom{T} = NamedTuple{
    (fieldnames(Atom{T})..., :frame_id),
    Tuple{fieldtypes(Atom{T})..., Int}
}

mutable struct Molecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame
    properties::Properties

    function Molecule{T}(name = "") where T
        new(name, DataFrame(_MoleculeAtom{T}[]), DataFrame(Bond[]), Properties())
    end
end

Molecule(name = "") = Molecule{Float32}(name)

@inline function atoms(ac::Molecule{T}; frame_id::Int = 1) where T
    df = get(groupby(ac.atoms, :frame_id), (frame_id = frame_id,), DataFrame(_MoleculeAtom{T}[]))
    view(df, :, 1:length(fieldnames(Atom{T})))
end

function Base.push!(ac::Molecule{T}, atom::Atom{T}; frame_id::Int = 1) where T
    push!(ac.atoms, (atom..., frame_id = frame_id))
end
