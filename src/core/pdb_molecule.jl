using DataFrames

export PDBMolecule

const _PDBMoleculeAtom{T} = NamedTuple{
    (fieldnames(Atom{T})..., :frame_id, :fragment_id),
    Tuple{fieldtypes(Atom{T})..., Int, Int}
}

mutable struct PDBMolecule{T<:Real} <: AbstractMolecule{T}
    name::String

    atoms::DataFrame
    bonds::DataFrame
    fragments::DataFrame
    properties::Properties

    function PDBMolecule{T}(name = "") where T
        new(name, DataFrame(_PDBMoleculeAtom{T}[]), DataFrame(Bond[]), DataFrame(Fragment[]), Properties())
    end
end

PDBMolecule(name = "") = PDBMolecule{Float32}(name)

@inline function fragments(ac::PDBMolecule{T}) where T
    ac.fragments
end

@inline function atoms(ac::PDBMolecule{T}; frame_id::Int = 1) where T
    df = get(groupby(ac.atoms, :frame_id), (frame_id = frame_id,), DataFrame(_PDBMoleculeAtom{T}[]))
    view(df, :, 1:length(fieldnames(Atom{T})))
end

@inline function atoms(ac::PDBMolecule{T}, f::Fragment; frame_id::Int = 1) where T
    df = get(
        groupby(ac.atoms, [:frame_id, :fragment_id]),
        (frame_id = frame_id, fragment_id = f.number),
        DataFrame(_PDBMoleculeAtom{T}[])
    )
    view(df, :, 1:length(fieldnames(Atom{T})))
end

function Base.push!(ac::PDBMolecule{T}, f::Fragment, atoms::Atom{T}...; frame_id::Int = 1) where T
    # check whether the given fragment is known, add it otherwise
    dff = ac.fragments
    if isempty(@view dff[dff.number .== f.number, :])
        push!(dff, f)
    end

    # add all atoms
    for atom in atoms
        push!(ac.atoms, (atom..., frame_id = frame_id, fragment_id = f.number))
    end
end
