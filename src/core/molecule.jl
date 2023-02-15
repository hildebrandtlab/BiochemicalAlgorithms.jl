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
    # if the DataFrame is not empty, it requires at least one molecule w/ frame_id=1
    df = nrow(ac.atoms) == 0 ? ac.atoms : groupby(ac.atoms, :frame_id)[(frame_id = frame_id,)]
    view(df, :, 1:length(fieldnames(Atom{T})))
end

function Base.push!(ac::Molecule{T}, atom::Atom{T}; frame_id::Int = 1) where T
    push!(ac.atoms, (atom..., frame_id = frame_id))
end
